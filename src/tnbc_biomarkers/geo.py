"""Download and parse GEO series matrices for a SuperSeries accession."""

from __future__ import annotations

import gzip
import hashlib
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

USER_AGENT = "tnbc-biomarkers/1.0 (Ishitta Sarkar; public GEO reanalysis)"
SOFT_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text&view=brief"


def geo_folder(accession: str) -> str:
    """NCBI FTP grouping directory for a GEO accession.

    GSE76124 is stored under GSE76nnn and GPL570 under GPLnnn.
    """

    prefix = "".join(character for character in accession if character.isalpha())
    digits = "".join(character for character in accession if character.isdigit())
    if len(digits) <= 3:
        grouped = "nnn"
    else:
        grouped = f"{digits[:-3]}nnn"
    return f"{prefix}{grouped}"


def series_matrix_url(accession: str) -> str:
    folder = geo_folder(accession)
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{folder}/{accession}/"
        f"matrix/{accession}_series_matrix.txt.gz"
    )


def platform_annotation_url(platform: str) -> str:
    folder = geo_folder(platform)
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{folder}/{platform}/"
        f"annot/{platform}.annot.gz"
    )


def parse_superseries_members(soft_text: str) -> list[str]:
    """Return SubSeries accessions declared by a SuperSeries SOFT record."""

    members: list[str] = []
    for line in soft_text.splitlines():
        if "SuperSeries of:" not in line:
            continue
        accession = line.split("SuperSeries of:", 1)[1].strip().strip('"')
        if accession.startswith("GSE") and accession not in members:
            members.append(accession)
    return members


def _request(url: str, stream: bool = False, timeout: float = 180) -> requests.Response:
    headers = {"User-Agent": USER_AGENT}
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            response = requests.get(url, headers=headers, stream=stream, timeout=timeout)
            response.raise_for_status()
            return response
        except requests.RequestException as error:
            last_error = error
            if attempt == 2:
                break
            time.sleep(2**attempt)
    raise RuntimeError(f"Failed to fetch {url}") from last_error


def fetch_text(url: str) -> str:
    response = _request(url, stream=False, timeout=120)
    response.encoding = "utf-8"
    return response.text


def download(url: str, dest: Path) -> Path:
    """Download ``url`` to ``dest`` unless a readable file is already cached."""

    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000 and _gzip_readable(dest):
        return dest
    if dest.exists():
        dest.unlink()
    temporary = dest.with_suffix(dest.suffix + ".part")
    response = _request(url, stream=True)
    try:
        with temporary.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1 << 20):
                if chunk:
                    handle.write(chunk)
    finally:
        response.close()
    temporary.replace(dest)
    if not _gzip_readable(dest):
        dest.unlink(missing_ok=True)
        raise RuntimeError(f"Downloaded file is not a readable gzip archive: {dest}")
    return dest


def _gzip_readable(path: Path) -> bool:
    if not str(path).endswith(".gz"):
        return path.exists() and path.stat().st_size > 0
    try:
        with gzip.open(path, "rb") as handle:
            handle.read(64)
    except OSError:
        return False
    return True


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _geo_fields(line: str) -> tuple[str, list[str]]:
    parts = line.rstrip("\n").split("\t")
    key = parts[0][1:] if parts[0].startswith("!") else parts[0]
    values = [part.strip().strip('"') for part in parts[1:]]
    return key, values


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    """Parse a GEO series matrix into sample metadata, a probe matrix, and series fields.

    Sample characteristics in a series matrix are aligned by row index, not by
    key: samples can list age, receptor status, and histology in different
    orders. Each cell is parsed as ``key: value`` on its own.
    """

    series_fields: dict[str, str] = {}
    sample_rows: dict[str, list[str]] = {}
    characteristic_rows: list[list[str]] = []
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                expression = pd.read_csv(handle, sep="\t", low_memory=False)
                break
            if not line.startswith("!"):
                continue
            key, values = _geo_fields(line)
            if key.startswith("Sample_"):
                if key == "Sample_characteristics_ch1":
                    characteristic_rows.append(values)
                elif key not in sample_rows:
                    sample_rows[key] = values
            elif key not in series_fields and values:
                series_fields[key] = values[0]
        else:
            raise ValueError(f"No expression table in {path}")

    expression = expression.rename(columns=lambda name: str(name).strip().strip('"'))
    expression = expression.loc[
        ~expression.iloc[:, 0].astype(str).str.startswith("!"), :
    ].copy()
    probe_ids = expression.iloc[:, 0].astype(str).str.strip().str.strip('"')
    expression = expression.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    expression.index = probe_ids.to_numpy()
    expression.index.name = "probe_id"
    expression = expression.loc[:, ~expression.columns.astype(str).str.startswith("Unnamed")]
    if expression.index.has_duplicates:
        expression = expression[~expression.index.duplicated(keep="first")]

    accessions = sample_rows.get("Sample_geo_accession", list(expression.columns))
    if len(accessions) != expression.shape[1]:
        raise ValueError(
            f"Sample count in metadata ({len(accessions)}) does not match the "
            f"expression table ({expression.shape[1]}) in {path}"
        )
    expression.columns = accessions
    records = []
    for index, gsm in enumerate(accessions):
        characteristics: dict[str, str] = {}
        for row in characteristic_rows:
            if index >= len(row):
                continue
            raw = row[index].strip()
            if not raw or ":" not in raw:
                continue
            key, value = raw.split(":", 1)
            characteristics[key.strip().lower()] = value.strip()
        processing_values = sample_rows.get("Sample_data_processing", [""] * len(accessions))
        records.append(
            {
                "sample_id": gsm,
                "title": sample_rows.get("Sample_title", [""] * len(accessions))[index],
                "series": series_fields.get("Series_geo_accession", ""),
                "data_processing": (
                    processing_values[index] if index < len(processing_values) else ""
                ),
                "characteristics": characteristics,
            }
        )
    samples = pd.DataFrame.from_records(records).set_index("sample_id")
    return samples, expression, series_fields


def annotation_date(path: Path) -> str | None:
    """Read the GEO platform annotation date stamped in the file header."""

    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("!Annotation_date"):
                return line.split("=", 1)[1].strip()
            if line.startswith("!platform_table_begin"):
                break
    return None


def load_gene_symbols(path: Path) -> pd.Series:
    """Load a probe-to-symbol map from a GEO platform annotation file."""

    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                table = pd.read_csv(
                    handle,
                    sep="\t",
                    usecols=["ID", "Gene symbol"],
                    dtype=str,
                    low_memory=False,
                )
                break
        else:
            raise ValueError(f"No platform table in {path}")
    table["ID"] = table["ID"].str.strip()
    table = table.dropna(subset=["ID"]).drop_duplicates("ID", keep="first")
    symbols = table.set_index("ID")["Gene symbol"]
    symbols.index.name = "probe_id"
    return symbols


def normalise_tn_status(value: str | None) -> str | None:
    """Map the GEO triple-negative field onto TNBC or non-TNBC."""

    if value is None:
        return None
    text = value.strip().lower().replace("_", " ").replace("-", " ")
    text = " ".join(text.split())
    if text in {"tn", "tnbc", "triple negative", "yes"}:
        return "TNBC"
    if text in {"not tn", "non tn", "non tnbc", "no"}:
        return "non-TNBC"
    return None


def assign_groups(samples: pd.DataFrame) -> pd.DataFrame:
    """Add receptor columns and the analysis label from sample characteristics."""

    frame = samples.copy()
    groups = []
    er_values = []
    pr_values = []
    her2_values = []
    raw_status = []
    for characteristics in frame["characteristics"]:
        status_raw = None
        for key, value in characteristics.items():
            if "triple" in key and "status" in key:
                status_raw = value
                break
        raw_status.append(status_raw)
        groups.append(normalise_tn_status(status_raw))
        er_values.append(characteristics.get("er"))
        pr_values.append(characteristics.get("pr"))
        her2_values.append(characteristics.get("her2"))
    frame["tn_status_raw"] = raw_status
    frame["group"] = groups
    frame["er"] = er_values
    frame["pr"] = pr_values
    frame["her2"] = her2_values
    frame["receptor_group"] = [
        _receptor_group(er, pr, her2) for er, pr, her2 in zip(er_values, pr_values, her2_values)
    ]
    return frame


def _receptor_group(er: str | None, pr: str | None, her2: str | None) -> str | None:
    values = [er, pr, her2]
    if any(value is None or value.strip() == "" for value in values):
        return None
    if all(str(value).strip().lower().startswith("neg") for value in values):
        return "TNBC"
    return "non-TNBC"


def collapse_to_genes(
    expression: pd.DataFrame, symbols: pd.Series
) -> tuple[pd.DataFrame, pd.Series, dict[str, int]]:
    """Collapse a probe matrix to one row per gene symbol.

    Control probes (AFFX*), probes without a single gene symbol, and ambiguous
    multi-symbol probes are removed. When several probes share a symbol, the
    probe with the highest mean expression is kept. Ties break toward the
    smaller probe identifier.
    """

    probes = expression.copy()
    probes = probes.loc[~probes.index.astype(str).str.upper().str.startswith("AFFX")]
    n_control = int(expression.shape[0] - probes.shape[0])
    mapped = symbols.reindex(probes.index)
    keep_symbol = pd.Series(
        [_single_symbol(value) for value in mapped.tolist()],
        index=mapped.index,
    )
    n_unmapped = int(keep_symbol.isna().sum())
    n_ambiguous = int(
        mapped.fillna("").astype(str).str.contains("///", regex=False).sum()
    )
    probes = probes.loc[keep_symbol.notna()].copy()
    gene_names = keep_symbol.loc[probes.index]
    missing = probes.isna().any(axis=1)
    n_missing = int(missing.sum())
    probes = probes.loc[~missing]
    gene_names = gene_names.loc[probes.index]
    means = probes.mean(axis=1)
    ranking = pd.DataFrame(
        {
            "probe_id": probes.index.astype(str),
            "gene": gene_names.to_numpy(),
            "mean_expression": means.to_numpy(),
        }
    )
    ranking = ranking.sort_values(
        ["gene", "mean_expression", "probe_id"],
        ascending=[True, False, True],
        kind="mergesort",
    )
    chosen = ranking.groupby("gene", as_index=False, sort=False).head(1)
    selected = probes.loc[chosen["probe_id"].to_numpy()]
    selected.index = chosen["gene"].to_numpy()
    selected.index.name = "gene"
    probe_for_gene = pd.Series(chosen["probe_id"].to_numpy(), index=chosen["gene"].to_numpy())
    probe_for_gene.index.name = "gene"
    counts = {
        "n_probes": int(expression.shape[0]),
        "n_control_probes_removed": n_control,
        "n_unmapped_or_ambiguous_probes": n_unmapped,
        "n_ambiguous_symbol_probes": n_ambiguous,
        "n_probes_with_missing_values": n_missing,
        "n_genes": int(selected.shape[0]),
    }
    return selected, probe_for_gene, counts


def _single_symbol(raw: object) -> str | None:
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return None
    text = str(raw).strip()
    if text in {"", "---", "NA", "N/A", "nan", "None"}:
        return None
    if "///" in text or any(character.isspace() for character in text):
        return None
    return text


def fetch_cohort(accession: str, cache_dir: Path) -> dict[str, object]:
    """Resolve ``accession`` and load every SubSeries matrix it points to.

    A series that is not a SuperSeries is loaded directly. Expression from the
    subseries is concatenated on probe identifiers. Files are cached under
    ``cache_dir`` and are not part of the repository.
    """

    cache_dir.mkdir(parents=True, exist_ok=True)
    soft_url = SOFT_URL.format(accession=accession)
    soft_text = fetch_text(soft_url)
    soft_path = cache_dir / f"{accession}_family_soft.txt"
    soft_path.write_text(soft_text, encoding="utf-8")
    members = parse_superseries_members(soft_text)
    if not members:
        members = [accession]
    sample_frames = []
    expression_frames = []
    platforms = set()
    processing = set()
    checksums = {}
    for member in members:
        url = series_matrix_url(member)
        path = download(url, cache_dir / f"{member}_series_matrix.txt.gz")
        checksums[member] = sha256_file(path)
        samples, expression, fields = parse_series_matrix(path)
        platform = fields.get("Series_platform_id", "")
        platforms.add(platform)
        processing.update(value for value in samples["data_processing"].dropna().unique() if value)
        sample_frames.append(samples)
        expression_frames.append(expression)
    if len(platforms) != 1:
        raise RuntimeError(f"Expected one platform, found {sorted(platforms)}")
    platform = next(iter(platforms))
    annotation_path = download(
        platform_annotation_url(platform),
        cache_dir / f"{platform}.annot.gz",
    )
    checksums[platform] = sha256_file(annotation_path)
    symbols = load_gene_symbols(annotation_path)
    dated = annotation_date(annotation_path)
    expression = pd.concat(expression_frames, axis=1, join="inner")
    if expression.columns.has_duplicates:
        raise RuntimeError("Duplicate GEO sample accessions across subseries")
    samples = pd.concat(sample_frames, axis=0)
    samples = assign_groups(samples)
    samples = samples.loc[expression.columns]
    return {
        "accession": accession,
        "subseries": members,
        "platform": platform,
        "soft_url": soft_url,
        "matrix_urls": {member: series_matrix_url(member) for member in members},
        "annotation_url": platform_annotation_url(platform),
        "checksums": checksums,
        "annotation_date": dated,
        "data_processing": sorted(processing),
        "samples": samples,
        "expression": expression,
        "symbols": symbols,
    }
