import numpy as np
from scipy import sparse

def check_counts_layer(
    adata,
    layer: str = "counts",
    sample_n: int = 200_000,
    int_tol: float = 1e-6,
    min_cells: int = 1,
    verbose: bool = True,
):
    """
    English comments only.

    Quick QC for AnnData count layer suitability (e.g., scVI/scANVI).
    Works for large sparse matrices without densifying.

    Returns
    -------
    ok: bool
        Whether the counts layer passes basic requirements (exists, same shape, non-negative, integer-like).
    report: dict
        Key metrics and flags for debugging / Methods reporting.
    """
    report = {
        "layer": layer,
        "ok": True,
        "errors": [],
        "warnings": [],
    }

    # ---------- Basic presence & shape ----------
    if layer not in adata.layers:
        report["ok"] = False
        report["errors"].append(f"Missing adata.layers['{layer}'].")
        return report["ok"], report

    mat = adata.layers[layer]
    report["counts_shape"] = tuple(mat.shape)
    report["X_shape"] = tuple(adata.X.shape)

    if mat.shape != adata.X.shape:
        report["ok"] = False
        report["errors"].append(f"Shape mismatch: layers['{layer}'] {mat.shape} != X {adata.X.shape}.")

    # ---------- Sparse/dtype ----------
    report["is_sparse"] = sparse.issparse(mat)
    report["dtype"] = str(mat.dtype)

    # ---------- Uniqueness checks (important for mapping/integration) ----------
    report["obs_names_unique"] = bool(adata.obs_names.is_unique)
    report["var_names_unique"] = bool(adata.var_names.is_unique)
    if not report["obs_names_unique"]:
        report["ok"] = False
        report["errors"].append("obs_names are not unique (cell IDs must be unique).")
    if not report["var_names_unique"]:
        report["ok"] = False
        report["errors"].append("var_names are not unique (gene IDs must be unique).")

    # ---------- Nonzero stats (fast on sparse) ----------
    if sparse.issparse(mat):
        nnz = mat.nnz
        data = mat.data
        report["nnz"] = int(nnz)
        report["density"] = float(nnz / (mat.shape[0] * mat.shape[1])) if mat.shape[0] * mat.shape[1] > 0 else np.nan
    else:
        # Avoid flattening gigantic dense matrices; sample via random indices if dense (rare here)
        data = mat.ravel()
        nnz = int(np.count_nonzero(data))
        report["nnz"] = nnz
        report["density"] = float(nnz / data.size) if data.size > 0 else np.nan

    if mat.shape[0] < min_cells:
        report["ok"] = False
        report["errors"].append(f"Too few cells: {mat.shape[0]} < {min_cells}.")

    # ---------- Value checks (sampled) ----------
    if nnz == 0:
        report["ok"] = False
        report["errors"].append("Counts layer has nnz=0 (all zeros).")
        return report["ok"], report

    # Sample nonzero entries for integer-like test
    if data.size > sample_n:
        rng = np.random.default_rng(0)
        idx = rng.choice(data.size, size=sample_n, replace=False)
        samp = data[idx]
    else:
        samp = data

    report["sample_nonzero_n"] = int(samp.size)
    report["min_nonzero"] = float(np.min(samp))
    report["max_nonzero"] = float(np.max(samp))
    report["n_negative_in_sample"] = int(np.sum(samp < 0))

    # Non-negative requirement
    if report["n_negative_in_sample"] > 0:
        report["ok"] = False
        report["errors"].append("Found negative values in counts (should be non-negative).")

    # Integer-like requirement (for raw counts)
    frac = np.abs(samp - np.round(samp))
    nonint_rate = float(np.mean(frac > int_tol))
    report["non_integer_rate_sample"] = nonint_rate
    if nonint_rate > 0.001:
        report["ok"] = False
        report["errors"].append(
            f"Counts look non-integer (non_integer_rate_sample={nonint_rate:.4f}). "
            "This may be normalized/log data instead of raw counts."
        )

    # ---------- Per-cell library size / detected genes (fast for sparse) ----------
    if sparse.issparse(mat):
        libsize = np.asarray(mat.sum(axis=1)).ravel()
        detected = mat.getnnz(axis=1)
    else:
        libsize = np.sum(mat, axis=1)
        detected = np.sum(mat > 0, axis=1)

    # Basic quantiles
    qs = [0, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999, 1.0]
    report["libsize_quantiles"] = {f"q{int(q*1000):04d}": float(np.quantile(libsize, q)) for q in qs}
    report["detected_genes_quantiles"] = {f"q{int(q*1000):04d}": float(np.quantile(detected, q)) for q in qs}

    # Heuristic warnings (not hard failures)
    if report["libsize_quantiles"]["q0000"] <= 0:
        report["warnings"].append("Some cells have zero total counts (may be filtered later).")
    if report["detected_genes_quantiles"]["q0000"] <= 0:
        report["warnings"].append("Some cells have zero detected genes (may be filtered later).")

    # ---------- Scanpy compatibility hints ----------
    has_neighbors = "neighbors" in adata.uns
    has_pca_uns = "pca" in adata.uns
    report["has_uns_neighbors"] = bool(has_neighbors)
    report["has_uns_pca"] = bool(has_pca_uns)
    if ("X_pca" in adata.obsm) and (not has_pca_uns):
        report["warnings"].append("X_pca exists in obsm but adata.uns['pca'] is missing (some Scanpy plotting/tools may expect it).")
    if ("X_umap" in adata.obsm) and (not has_neighbors):
        report["warnings"].append("X_umap exists in obsm but adata.uns['neighbors'] is missing (some Scanpy tools expect neighbors graph).")

    if verbose:
        print(f"[check_counts_layer] layer='{layer}' ok={report['ok']}")
        if report["errors"]:
            print("Errors:")
            for e in report["errors"]:
                print("  -", e)
        if report["warnings"]:
            print("Warnings:")
            for w in report["warnings"]:
                print("  -", w)
        print("Counts shape:", report["counts_shape"], "| sparse:", report["is_sparse"], "| dtype:", report["dtype"])
        print("Sample nonzero min/max:", report["min_nonzero"], report["max_nonzero"])
        print("Non-integer rate (sample):", report["non_integer_rate_sample"])
        print("Libsize q50/q99:", report["libsize_quantiles"]["q0500"], report["libsize_quantiles"]["q0990"])
        print("Detected genes q50/q99:", report["detected_genes_quantiles"]["q0500"], report["detected_genes_quantiles"]["q0990"])

    return report["ok"], report
