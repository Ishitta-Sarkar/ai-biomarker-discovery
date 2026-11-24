"""Nested cross-validation for an elastic-net logistic classifier."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from tnbc_biomarkers.config import AnalysisConfig
from tnbc_biomarkers.statistics import welch_test

TNBC = "TNBC"


def select_feature_indices(x_train: np.ndarray, y_train: np.ndarray, k: int) -> np.ndarray:
    """Rank training columns by |Welch t| and return the top ``k`` indexes.

    The statistic is computed only on ``x_train``. Non-finite statistics sort last.
    """

    positive = x_train[y_train == 1]
    negative = x_train[y_train == 0]
    statistic, _, _ = welch_test(positive, negative)
    statistic = np.abs(np.asarray(statistic, dtype=float))
    statistic[~np.isfinite(statistic)] = -1.0
    k = min(k, statistic.size)
    order = np.argsort(-statistic, kind="mergesort")
    return order[:k]


def classification_metrics(
    y_true: np.ndarray,
    scores: np.ndarray,
    predictions: np.ndarray | None = None,
) -> dict[str, float]:
    """Accuracy, balanced accuracy, F1, and ROC AUC for one labelled set."""

    y_true = np.asarray(y_true)
    scores = np.asarray(scores, dtype=float)
    if predictions is None:
        predictions = (scores >= 0.5).astype(int)
    return _fold_metrics(y_true, scores, np.asarray(predictions))


def _fold_metrics(
    y_true: np.ndarray, scores: np.ndarray, predictions: np.ndarray
) -> dict[str, float]:
    metrics = {
        "accuracy": float(accuracy_score(y_true, predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, predictions)),
        "f1": float(f1_score(y_true, predictions, zero_division=0)),
    }
    if len(np.unique(y_true)) < 2 or len(np.unique(scores)) < 2:
        metrics["roc_auc"] = 0.5
    else:
        metrics["roc_auc"] = float(roc_auc_score(y_true, scores))
    return metrics


def _summarise(
    fold_rows: list[dict[str, float]], pooled_y: np.ndarray, pooled_scores: np.ndarray
) -> dict:
    frame = pd.DataFrame(fold_rows)
    metric_names = ["accuracy", "balanced_accuracy", "roc_auc", "f1"]
    mean = {name: float(frame[name].mean()) for name in metric_names}
    std = {name: float(frame[name].std(ddof=1)) if len(frame) > 1 else 0.0 for name in metric_names}
    if len(np.unique(pooled_y)) < 2 or len(np.unique(pooled_scores)) < 2:
        pooled_auc = 0.5
    else:
        pooled_auc = float(roc_auc_score(pooled_y, pooled_scores))
    return {
        "folds": fold_rows,
        "mean": mean,
        "std": std,
        "pooled_roc_auc": pooled_auc,
    }


def _logistic(
    config: AnalysisConfig, elastic_net: bool, inner: StratifiedKFold
) -> LogisticRegressionCV:
    common = dict(
        Cs=np.asarray(config.cs, dtype=float),
        cv=inner,
        scoring="roc_auc",
        class_weight="balanced",
        max_iter=4000,
        tol=1e-3,
        random_state=config.random_state,
        refit=True,
        n_jobs=1,
        use_legacy_attributes=False,
    )
    if elastic_net:
        return LogisticRegressionCV(
            l1_ratios=config.l1_ratios,
            solver="saga",
            **common,
        )
    return LogisticRegressionCV(
        l1_ratios=(0.0,),
        solver="lbfgs",
        **common,
    )


def _scalar(value: object) -> float:
    array = np.asarray(value, dtype=float).reshape(-1)
    return float(array[0])


def nested_classification(
    expression: pd.DataFrame,
    groups: pd.Series,
    config: AnalysisConfig,
) -> dict[str, object]:
    """Compare elastic-net and L2 logistic models with a majority-class baseline.

    Feature selection is repeated inside each outer training fold: the top
    ``config.n_features`` genes by absolute Welch t on that fold only. Inner
    cross-validation selects the inverse regularisation strength ``C`` and, for
    the elastic net, the L1 ratio. The positive class is TNBC.
    """

    groups = groups.loc[expression.index]
    y = (groups == TNBC).astype(int).to_numpy()
    x = expression.to_numpy(dtype=float)
    genes = expression.columns.astype(str).to_numpy()
    if y.sum() < config.outer_splits or (len(y) - y.sum()) < config.outer_splits:
        raise ValueError("Not enough samples per class for the outer stratified folds")
    outer = StratifiedKFold(
        n_splits=config.outer_splits,
        shuffle=True,
        random_state=config.random_state,
    )
    model_names = ("elastic_net", "l2_logistic", "majority_baseline")
    fold_rows: dict[str, list[dict]] = {name: [] for name in model_names}
    pooled_y: list[np.ndarray] = []
    pooled_scores: dict[str, list[np.ndarray]] = {name: [] for name in model_names}
    coefficient_rows: list[pd.DataFrame] = []
    convergence_warnings = 0
    for fold, (train, test) in enumerate(outer.split(x, y), start=1):
        selected = select_feature_indices(x[train], y[train], config.n_features)
        scaler = StandardScaler()
        x_train = scaler.fit_transform(x[train][:, selected])
        x_test = scaler.transform(x[test][:, selected])
        y_train = y[train]
        y_test = y[test]
        inner = StratifiedKFold(
            n_splits=config.inner_splits,
            shuffle=True,
            random_state=config.random_state,
        )
        fitted = {}
        for name, elastic in (("elastic_net", True), ("l2_logistic", False)):
            model = _logistic(config, elastic, inner)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", ConvergenceWarning)
                model.fit(x_train, y_train)
                convergence_warnings += sum(
                    issubclass(item.category, ConvergenceWarning) for item in caught
                )
            fitted[name] = model
            scores = model.predict_proba(x_test)[:, 1]
            predictions = (scores >= 0.5).astype(int)
            metrics = _fold_metrics(y_test, scores, predictions)
            metrics.update(
                {
                    "fold": fold,
                    "C": _scalar(model.C_),
                    "l1_ratio": _scalar(model.l1_ratio_),
                    "n_nonzero": int(np.sum(np.abs(model.coef_) > config.nonzero_coef)),
                }
            )
            fold_rows[name].append(metrics)
            pooled_scores[name].append(scores)
        majority = int(np.bincount(y_train).argmax())
        prevalence = float(y_train.mean())
        majority_scores = np.full(y_test.shape, prevalence)
        majority_pred = np.full(y_test.shape, majority)
        baseline_metrics = _fold_metrics(y_test, majority_scores, majority_pred)
        baseline_metrics.update(
            {
                "fold": fold,
                "majority_class": majority,
                "training_prevalence": prevalence,
            }
        )
        fold_rows["majority_baseline"].append(baseline_metrics)
        pooled_scores["majority_baseline"].append(majority_scores)
        pooled_y.append(y_test)
        coefficients = np.asarray(fitted["elastic_net"].coef_, dtype=float).reshape(-1)
        coefficient_rows.append(
            pd.DataFrame(
                {
                    "gene": genes[selected],
                    "coefficient": coefficients,
                    "fold": fold,
                }
            )
        )
    pooled_y_array = np.concatenate(pooled_y)
    summary: dict[str, object] = {
        "positive_class": TNBC,
        "n_samples": int(len(y)),
        "n_tnbc": int(y.sum()),
        "n_non_tnbc": int(len(y) - y.sum()),
        "outer_splits": config.outer_splits,
        "inner_splits": config.inner_splits,
        "n_features": config.n_features,
        "cs": list(config.cs),
        "l1_ratios": list(config.l1_ratios),
        "random_state": config.random_state,
        "convergence_warnings": convergence_warnings,
        "class_weight": "balanced",
        "decision_threshold": 0.5,
    }
    for name in model_names:
        summary[name] = _summarise(
            fold_rows[name], pooled_y_array, np.concatenate(pooled_scores[name])
        )
    coefficients = pd.concat(coefficient_rows, ignore_index=True)
    summary["coefficients"] = _coefficient_stability(coefficients, config)
    return summary


def _coefficient_stability(coefficients: pd.DataFrame, config: AnalysisConfig) -> pd.DataFrame:
    coefficients = coefficients.copy()
    coefficients["nonzero"] = coefficients["coefficient"].abs() > config.nonzero_coef
    grouped = coefficients.groupby("gene", sort=False)
    stability = grouped.agg(
        n_folds_nonzero=("nonzero", "sum"),
        mean_coefficient=("coefficient", "mean"),
        mean_abs_coefficient=("coefficient", lambda values: float(np.mean(np.abs(values)))),
    ).reset_index()
    stability["n_folds_nonzero"] = stability["n_folds_nonzero"].astype(int)
    stability = stability.sort_values(
        ["n_folds_nonzero", "mean_abs_coefficient", "gene"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    return stability
