"""Pathway enrichment through the g:Profiler public API."""

from __future__ import annotations

import time

import pandas as pd
import requests

from tnbc_biomarkers.geo import USER_AGENT

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
TERM_COLUMNS = [
    "source",
    "native",
    "name",
    "p_value",
    "term_size",
    "query_size",
    "intersection_size",
    "precision",
    "recall",
]


def parse_gprofiler_response(payload: dict) -> tuple[pd.DataFrame, dict]:
    """Turn a g:Profiler JSON body into a term table and a short metadata record."""

    meta = payload.get("meta", {})
    rows = []
    for item in payload.get("result") or []:
        if not item.get("significant", False):
            continue
        rows.append({column: item.get(column) for column in TERM_COLUMNS})
    table = pd.DataFrame(rows, columns=TERM_COLUMNS)
    if not table.empty:
        table = table.sort_values(["p_value", "native"], kind="mergesort").reset_index(drop=True)
    record = {
        "version": meta.get("version"),
        "timestamp": meta.get("timestamp"),
        "query_metadata": meta.get("query_metadata", {}),
        "result_metadata": meta.get("result_metadata", {}),
    }
    return table, record


def gprofiler_enrich(
    genes: list[str],
    background: list[str],
    sources: tuple[str, ...] | list[str],
    session: requests.Session | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Query g:Profiler with a custom background of genes tested in this analysis.

    If the service rejects the custom background, the same query is sent against
    the annotated genome and the metadata records that fallback.
    """

    if not genes:
        empty, record = parse_gprofiler_response({"result": [], "meta": {}})
        record["domain_scope"] = "custom"
        record["n_query_genes"] = 0
        record["n_background"] = len(background)
        record["background_fallback"] = False
        return empty, record
    query = list(dict.fromkeys(genes))
    background_genes = list(dict.fromkeys([*query, *background]))
    scopes = ("custom", "annotated")
    last_error: Exception | None = None
    for scope in scopes:
        body = {
            "organism": "hsapiens",
            "query": query,
            "sources": list(sources),
            "user_threshold": 0.05,
            "significance_threshold_method": "g_SCS",
            "no_evidences": True,
            "domain_scope": scope,
        }
        if scope == "custom":
            body["background"] = background_genes
        try:
            payload = _post_json(body, session=session)
        except requests.RequestException as error:
            last_error = error
            continue
        table, record = parse_gprofiler_response(payload)
        record["domain_scope"] = scope
        record["n_query_genes"] = len(query)
        record["n_background"] = len(background_genes) if scope == "custom" else None
        record["background_fallback"] = scope != "custom"
        return table, record
    raise RuntimeError("g:Profiler enrichment request failed") from last_error


def _post_json(body: dict, session: requests.Session | None) -> dict:
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    sender = session.post if session is not None else requests.post
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            response = sender(GPROFILER_URL, json=body, headers=headers, timeout=180)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as error:
            last_error = error
            if attempt == 2:
                break
            time.sleep(2**attempt)
    raise requests.RequestException("g:Profiler request failed") from last_error
