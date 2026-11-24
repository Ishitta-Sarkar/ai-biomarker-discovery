"""GEO accession layout and series-matrix parsing."""

from pathlib import Path

import pandas as pd

from tnbc_biomarkers.geo import (
    assign_groups,
    collapse_to_genes,
    geo_folder,
    load_gene_symbols,
    normalise_tn_status,
    parse_series_matrix,
    parse_superseries_members,
    platform_annotation_url,
    series_matrix_url,
)

FIXTURES = Path(__file__).parent / "fixtures"


def test_ftp_layout():
    assert geo_folder("GSE76124") == "GSE76nnn"
    assert geo_folder("GSE76275") == "GSE76nnn"
    assert geo_folder("GPL570") == "GPLnnn"
    assert series_matrix_url("GSE76124") == (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76124/"
        "matrix/GSE76124_series_matrix.txt.gz"
    )
    assert platform_annotation_url("GPL570") == (
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
    )


def test_superseries_members_keep_order():
    text = "\n".join(
        [
            "!Series_relation = SuperSeries of: GSE76124",
            "!Series_relation = SuperSeries of: GSE76274",
            "!Series_relation = SuperSeries of: GSE76124",
        ]
    )
    assert parse_superseries_members(text) == ["GSE76124", "GSE76274"]


def test_characteristics_are_parsed_per_sample_not_by_row():
    samples, expression, fields = parse_series_matrix(FIXTURES / "mini_series_matrix.txt")
    labelled = assign_groups(samples)
    assert fields["Series_platform_id"] == "GPL570"
    assert list(expression.columns) == ["GSM1", "GSM2", "GSM3"]
    assert expression.loc["1007_s_at", "GSM2"] == 2.5
    assert "AFFX-BioB" in expression.index
    assert labelled.loc["GSM1", "group"] == "TNBC"
    assert labelled.loc["GSM2", "group"] == "non-TNBC"
    assert labelled.loc["GSM2", "er"] == "Positive"
    assert labelled.loc["GSM3", "group"] == "TNBC"
    assert labelled.loc["GSM3", "her2"] == "Negative"
    assert normalise_tn_status("not TN") == "non-TNBC"


def test_gene_collapse_drops_controls_and_ambiguous_probes():
    symbols = load_gene_symbols(FIXTURES / "mini_platform.annot.txt")
    expression = pd.DataFrame(
        {
            "s1": [1.0, 5.0, 9.0, 3.0, 8.0],
            "s2": [2.0, 5.0, 9.0, 4.0, 1.0],
        },
        index=["1007_s_at", "1053_at", "AFFX-BioB", "2000_at", "2001_at"],
    )
    collapsed, probes, counts = collapse_to_genes(expression, symbols)
    assert list(collapsed.index) == ["DDR1", "RFC2"]
    assert probes["DDR1"] == "2001_at"
    assert counts["n_genes"] == 2
    assert counts["n_control_probes_removed"] == 1
    assert counts["n_ambiguous_symbol_probes"] == 1
