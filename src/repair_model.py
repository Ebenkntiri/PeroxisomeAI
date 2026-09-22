"""Repaired evaluation core. Biological evaluation is NOT RUN.

This module is a software repair, not a validated predictor or trained model.
It requires an independently reviewed grouping/design record before real use.
All learned profiles, tuning, and sigmoid calibration remain inside training
partitions. Each calibration-training partition has its own hyperparameter search.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedGroupKFold, GridSearchCV
from sklearn.metrics import make_scorer, matthews_corrcoef
from sklearn.utils.validation import check_is_fitted

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {c: i for i, c in enumerate(AA)}


def validate_sequences(sequences, window):
    if not isinstance(window, (int, np.integer)) or window <= 0:
        raise ValueError("window must be a positive integer")
    out = []
    for i, seq in enumerate(sequences):
        if not isinstance(seq, str) or len(seq) < window:
            raise ValueError(f"Record {i}: missing sequence or shorter than {window}")
        if any(c not in AA_INDEX for c in seq):
            raise ValueError(f"Record {i}: invalid/non-uppercase residue; do not silently delete residues")
        out.append(seq)
    if not out:
        raise ValueError("Empty sequence set")
    return np.asarray(out, dtype=object)


class BiProfileEncoder(TransformerMixin, BaseEstimator):
    def __init__(self, window=30, pseudocount=1.0):
        self.window = window
        self.pseudocount = pseudocount

    def fit(self, X, y):
        X = validate_sequences(X, self.window)
        y = np.asarray(y)
        if len(y) != len(X) or set(np.unique(y)) != {0, 1}:
            raise ValueError("Training labels must contain both binary classes")
        if self.pseudocount <= 0:
            raise ValueError("pseudocount must be positive")
        self.profiles_ = []
        for label in (1, 0):
            counts = np.full((self.window, 20), self.pseudocount, dtype=float)
            for seq in X[y == label]:
                for pos, aa in enumerate(seq[-self.window:]):
                    counts[pos, AA_INDEX[aa]] += 1
            self.profiles_.append(counts / counts.sum(axis=1, keepdims=True))
        self.n_fit_records_ = len(X)
        return self

    def transform(self, X):
        check_is_fitted(self, "profiles_")
        X = validate_sequences(X, self.window)
        indices = np.array([[AA_INDEX[c] for c in seq[-self.window:]] for seq in X])
        positions = np.arange(self.window)
        return np.concatenate([np.log(p[positions, indices]) for p in self.profiles_], axis=1)


def checked_splits(X, y, groups, n_splits, seed):
    """Preserve group boundaries and fail explicitly on unusable folds."""
    X, y, groups = np.asarray(X), np.asarray(y), np.asarray(groups)
    if not (len(X) == len(y) == len(groups)):
        raise ValueError("Misaligned sequences, labels, and groups")
    if len(np.unique(groups)) < n_splits:
        raise ValueError("Too few groups for requested split; do not silently weaken the design")
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    splits = list(cv.split(X, y, groups))
    for tr, te in splits:
        if set(groups[tr]) & set(groups[te]):
            raise AssertionError("Group leakage")
        if set(y[tr]) != {0, 1} or set(y[te]) != {0, 1}:
            raise ValueError("A fold lacks a class; revise the prespecified feasible design")
    return splits


def fit_sigmoid(scores, y):
    """Platt-style sigmoid with class-count target smoothing.

    Fits q = sigmoid(a * decision_score + b) on calibration data only.
    This is not evidence of probability calibration in a target proteome.
    """
    scores, y = np.asarray(scores, float), np.asarray(y, int)
    n1, n0 = np.sum(y == 1), np.sum(y == 0)
    if min(n1, n0) == 0:
        raise ValueError("Calibration requires both classes")
    target = np.where(y == 1, (n1 + 1) / (n1 + 2), 1 / (n0 + 2))
    def loss(v):
        z = v[0] * scores + v[1]
        return np.sum(np.logaddexp(0, z) - target * z)
    def jac(v):
        delta = expit(v[0] * scores + v[1]) - target
        return np.array([np.dot(delta, scores), delta.sum()])
    result = minimize(loss, np.array([0., np.log((n1 + 1) / (n0 + 1))]),
                      jac=jac, method="L-BFGS-B", options={"maxiter": 1000, "ftol": 1e-12})
    if not result.success or not np.isfinite(result.x).all():
        raise RuntimeError(f"Sigmoid calibration failed: {result.message}")
    return result.x


def nested_cv(sequences, labels, groups, *, window=30, outer_k=5,
              inner_k=5, calibration_k=3, seed=42, grid=None):
    """Software core; call only after label, sequence and relatedness review.

    Groups must reflect the declared generalization target. Unique dummy groups
    are used only in synthetic tests, never as a substitute for a biological audit.
    Outputs do not establish generalization until the external preflight passes.
    """
    X = validate_sequences(sequences, window)
    y, groups = np.asarray(labels, int), np.asarray(groups)
    if set(y) != {0, 1}:
        raise ValueError("Both binary classes required")
    if grid is None:
        grid = {"svm__C": [0.1, 1, 10, 100, 1000],
                "svm__gamma": [0.001, 0.01, 0.1, 1]}
    prob = np.full(len(y), np.nan)
    folds = np.full(len(y), -1)
    trace = []
    for fold, (outer_tr, outer_te) in enumerate(checked_splits(X, y, groups, outer_k, seed)):
        probabilities = []
        cal_splits = checked_splits(X[outer_tr], y[outer_tr], groups[outer_tr],
                                    calibration_k, seed + 100 + fold)
        for cal_fold, (fit_local, cal_local) in enumerate(cal_splits):
            fit_idx, cal_idx = outer_tr[fit_local], outer_tr[cal_local]
            inner = checked_splits(X[fit_idx], y[fit_idx], groups[fit_idx],
                                   inner_k, seed + 1000 + fold * calibration_k + cal_fold)
            model = Pipeline([("profile", BiProfileEncoder(window)),
                              ("svm", SVC(kernel="rbf", class_weight="balanced"))])
            search = GridSearchCV(model, grid, scoring=make_scorer(matthews_corrcoef),
                                  cv=inner, n_jobs=1, error_score="raise", refit=True)
            search.fit(X[fit_idx], y[fit_idx])
            a, b = fit_sigmoid(search.decision_function(X[cal_idx]), y[cal_idx])
            probabilities.append(expit(a * search.decision_function(X[outer_te]) + b))
            trace.append({"outer_fold": fold, "calibration_fold": cal_fold,
                          "fit_indices": fit_idx.tolist(), "calibration_indices": cal_idx.tolist(),
                          "test_indices": outer_te.tolist(), "best_params": search.best_params_,
                          "sigmoid": [float(a), float(b)],
                          "inner_splits": [{"train": fit_idx[it].tolist(), "validation": fit_idx[iv].tolist()}
                                           for it, iv in inner]})
        prob[outer_te] = np.mean(probabilities, axis=0)
        folds[outer_te] = fold
    if np.isnan(prob).any() or np.any(folds < 0):
        raise AssertionError("Incomplete out-of-fold predictions")
    return {"pred": (prob >= .5).astype(int), "prob": prob, "fold": folds, "trace": trace}
