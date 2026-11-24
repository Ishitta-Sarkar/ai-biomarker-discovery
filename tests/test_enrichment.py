"""g:Profiler response parsing and request body."""

import pandas as pd

from tnbc_biomarkers.enrichment import gprofiler_enrich, parse_gprofiler_response


class _Response:
    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self.payload


class _Session:
    def __init__(self, payload):
        self.payload = payload
        self.bodies = []

    def post(self, url, json, headers, timeout):
        self.bodies.append(json)
        return _Response(self.payload)


def test_parser_keeps_significant_terms_and_version():
    payload = {
        "meta": {
            "version": "test-version",
            "timestamp": "2026-11-24T00:00:00+00:00",
            "query_metadata": {},
            "result_metadata": {},
        },
        "result": [
            {
                "source": "GO:BP",
                "native": "GO:0000002",
                "name": "later",
                "p_value": 0.02,
                "term_size": 10,
                "query_size": 3,
                "intersection_size": 2,
                "precision": 0.5,
                "recall": 0.2,
                "significant": True,
            },
            {
                "source": "KEGG",
                "native": "KEGG:00010",
                "name": "earlier",
                "p_value": 0.001,
                "term_size": 20,
                "query_size": 3,
                "intersection_size": 2,
                "precision": 0.5,
                "recall": 0.1,
                "significant": True,
            },
            {
                "source": "REAC",
                "native": "REAC:R1",
                "name": "hidden",
                "p_value": 0.9,
                "term_size": 5,
                "query_size": 3,
                "intersection_size": 1,
                "precision": 0.2,
                "recall": 0.1,
                "significant": False,
            },
        ],
    }
    table, meta = parse_gprofiler_response(payload)
    assert meta["version"] == "test-version"
    assert list(table["native"]) == ["KEGG:00010", "GO:0000002"]
    assert table.iloc[0]["p_value"] == 0.001


def test_request_uses_custom_background_and_unique_genes():
    session = _Session({"result": [], "meta": {"version": "v", "timestamp": "t"}})
    table, meta = gprofiler_enrich(
        ["ESR1", "ESR1", "FOXA1"],
        ["ESR1", "KRT5"],
        ("GO:BP", "KEGG", "REAC"),
        session=session,
    )
    assert isinstance(table, pd.DataFrame)
    assert table.empty
    body = session.bodies[0]
    assert body["domain_scope"] == "custom"
    assert body["query"] == ["ESR1", "FOXA1"]
    assert body["background"] == ["ESR1", "FOXA1", "KRT5"]
    assert meta["background_fallback"] is False
    assert meta["n_query_genes"] == 2


def test_empty_query_does_not_call_the_service():
    session = _Session({"result": [], "meta": {}})
    table, meta = gprofiler_enrich([], ["ESR1"], ("GO:BP",), session=session)
    assert table.empty
    assert session.bodies == []
    assert meta["n_query_genes"] == 0
