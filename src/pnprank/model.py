"""Model fitting with inner cross-validation on pre-cutoff genes only.

Two model families:
- logistic: median imputation with missing-value indicators, standardisation, L2
  logistic regression with balanced class weights; C tuned by inner CV.
- hgb: histogram gradient boosting (native NaN handling), balanced class weights;
  learning rate and depth tuned by inner CV.

Label-derived features are passed in as a callable so they can be recomputed from the
training labels of every fold; nothing here ever sees held-out or post-cutoff labels.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

GRIDS = {
    "logistic": {"logisticregression__C": [0.001, 0.01, 0.1, 1.0]},
    "hgb": {"learning_rate": [0.03, 0.1], "max_depth": [2, 3], "l2_regularization": [0.0, 1.0]},
}


def make_model(kind: str, seed: int):
    if kind == "logistic":
        return make_pipeline(
            SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
            StandardScaler(),
            LogisticRegression(class_weight="balanced", max_iter=5000),
        )
    if kind == "hgb":
        return HistGradientBoostingClassifier(
            class_weight="balanced", max_iter=200, random_state=seed
        )
    raise ValueError(f"unknown model kind {kind}")


def fit_tuned(X: pd.DataFrame, y: np.ndarray, kind: str, seed: int, folds: int):
    """Grid-search hyperparameters by inner stratified CV on (X, y) only; refit on all."""
    cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    gs = GridSearchCV(make_model(kind, seed), GRIDS[kind], scoring="roc_auc", cv=cv, n_jobs=1)
    gs.fit(X, y)
    return gs.best_estimator_, gs.best_params_


LabelFeatureFn = Callable[[set[str]], pd.DataFrame]


def with_label_features(
    X_static: pd.DataFrame, label_fn: LabelFeatureFn | None, known: set[str]
) -> pd.DataFrame:
    if label_fn is None:
        return X_static
    return pd.concat([X_static, label_fn(known)], axis=1)


def cross_validated_scores(
    X_static: pd.DataFrame,
    y: pd.Series,
    kind: str,
    seed: int,
    outer_folds: int,
    inner_folds: int,
    label_fn: LabelFeatureFn | None = None,
) -> tuple[pd.Series, list[dict]]:
    """Out-of-fold scores by nested CV.

    For each outer fold, label-derived features are recomputed from that fold's training
    positives only, hyperparameters are tuned by inner CV on the training fold, and the
    held-out genes are scored. Returns (out-of-fold scores, per-fold best params).
    """
    oof = pd.Series(np.nan, index=y.index)
    params = []
    outer = StratifiedKFold(n_splits=outer_folds, shuffle=True, random_state=seed)
    for k, (tr, te) in enumerate(outer.split(X_static, y)):
        train_idx, test_idx = y.index[tr], y.index[te]
        known = set(train_idx[y.iloc[tr].to_numpy().astype(bool)])
        X = with_label_features(X_static, label_fn, known)
        est, best = fit_tuned(X.loc[train_idx], y.loc[train_idx].to_numpy(), kind, seed + k, inner_folds)
        oof.loc[test_idx] = est.predict_proba(X.loc[test_idx])[:, 1]
        params.append({"fold": k, **best})
    return oof, params
