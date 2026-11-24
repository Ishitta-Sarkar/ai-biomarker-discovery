"""TCGA-BRCA expression and IHC labels from UCSC Xena."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tnbc_biomarkers.geo import download, sha256_file

EXPRESSION_URL = "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz"
CLINICAL_URL = "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
EXPRESSION_NOTE = (
    "UCSC Xena TCGA.BRCA.sampleMap/HiSeqV2, log2(normalized RSEM count + 1), "
    "one matrix for the BRCA cohort"
)

ER_COLUMN = "breast_carcinoma_estrogen_receptor_status"
PR_COLUMN = "breast_carcinoma_progesterone_receptor_status"
HER2_COLUMN = "lab_proc_her2_neu_immunohistochemistry_receptor_status"


def _status(value: object) -> str | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip().lower()
    missing = {
        "",
        "nan",
        "none",
        "indeterminate",
        "[not available]",
        "[not evaluated]",
        "[unknown]",
    }
    if text in missing:
        return None
    if text in {"positive", "negative"}:
        return text
    return None


def ihc_group(er: object, pr: object, her2: object) -> str | None:
    """TNBC if ER, PR, and HER2 IHC are all Negative.

    Equivocal, indeterminate, and missing HER2 calls are excluded. A tumor is
    non-TNBC only when every receptor is Positive or Negative and at least one
    is Positive.
    """

    er_status = _status(er)
    pr_status = _status(pr)
    her2_status = _status(her2)
    if None in {er_status, pr_status, her2_status}:
        return None
    if er_status == "negative" and pr_status == "negative" and her2_status == "negative":
        return "TNBC"
    return "non-TNBC"


def label_primary_tumors(clinical: pd.DataFrame) -> pd.DataFrame:
    """Primary tumors with a complete ER/PR/HER2 IHC call."""

    required = ["sampleID", "sample_type", ER_COLUMN, PR_COLUMN, HER2_COLUMN]
    missing = [column for column in required if column not in clinical.columns]
    if missing:
        raise ValueError(f"Clinical matrix is missing {missing}")
    primary = clinical.loc[clinical["sample_type"].astype(str) == "Primary Tumor"].copy()
    primary["group"] = [
        ihc_group(er, pr, her2)
        for er, pr, her2 in zip(primary[ER_COLUMN], primary[PR_COLUMN], primary[HER2_COLUMN])
    ]
    columns = ["sampleID", "group", ER_COLUMN, PR_COLUMN, HER2_COLUMN]
    labelled = primary.loc[primary["group"].notna(), columns]
    labelled = labelled.drop_duplicates("sampleID", keep="first").set_index("sampleID")
    labelled.index.name = "sample_id"
    return labelled


def fetch_tcga_brca(cache_dir: Path) -> dict[str, object]:
    """Download the Xena HiSeqV2 matrix and IHC labels, then keep labelled tumors."""

    cache_dir.mkdir(parents=True, exist_ok=True)
    clinical_path = download(CLINICAL_URL, cache_dir / "TCGA_BRCA_clinicalMatrix.tsv")
    expression_path = download(EXPRESSION_URL, cache_dir / "TCGA_BRCA_HiSeqV2.gz")
    clinical = pd.read_csv(clinical_path, sep="\t", low_memory=False)
    labels = label_primary_tumors(clinical)
    matrix = pd.read_csv(expression_path, sep="\t", index_col=0, compression="gzip")
    matrix.index = matrix.index.astype(str).str.strip()
    matrix = matrix.groupby(level=0).mean()
    shared = [sample for sample in labels.index if sample in matrix.columns]
    expression = matrix.loc[:, shared].T
    expression.index.name = "sample_id"
    expression = expression.apply(pd.to_numeric, errors="coerce")
    finite_genes = expression.columns[expression.notna().all(axis=0)]
    expression = expression.loc[:, finite_genes]
    groups = labels.loc[shared, "group"]
    return {
        "expression": expression,
        "groups": groups,
        "clinical": labels.loc[shared],
        "meta": {
            "cohort": "TCGA-BRCA",
            "expression_url": EXPRESSION_URL,
            "clinical_url": CLINICAL_URL,
            "expression_processing": EXPRESSION_NOTE,
            "er_column": ER_COLUMN,
            "pr_column": PR_COLUMN,
            "her2_column": HER2_COLUMN,
            "n_clinical_rows": int(clinical.shape[0]),
            "n_labelled_primary_tumors": int(labels.shape[0]),
            "n_with_expression": int(len(shared)),
            "checksums": {
                "clinical": sha256_file(clinical_path),
                "expression": sha256_file(expression_path),
            },
        },
    }
