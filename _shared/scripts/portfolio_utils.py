"""
Portfolio-wide helpers — used across all 10 projects.

Import in any notebook with:
    import sys; sys.path.append('../../_shared/scripts')
    import portfolio_utils as pu
"""

from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc


def setup_plotting(style_path: str | Path = "../../_shared/styles/portfolio.mplstyle") -> None:
    """Apply portfolio-wide matplotlib style and Scanpy figure params."""
    style_path = Path(style_path)
    if style_path.exists():
        plt.style.use(str(style_path))
    sc.settings.set_figure_params(
        dpi=120,
        dpi_save=300,
        frameon=False,
        figsize=(5, 4),
        format="png",
    )
    sc.settings.verbosity = 1


def save_figure(fig, name: str, out_dir: str | Path = "../figures",
                formats=("png", "tiff", "pdf"), dpi: int = 300) -> None:
    """Save a figure in multiple formats (PNG for web, TIFF for journals, PDF for vector)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        path = out_dir / f"{name}.{fmt}"
        save_kwargs = {"dpi": dpi, "bbox_inches": "tight"}
        if fmt == "tiff":
            save_kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, **save_kwargs)
        print(f"  saved {path}")


def quick_qc_summary(adata) -> dict:
    """Return a one-line QC summary dictionary for an AnnData object."""
    return {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "median_counts": float(adata.obs["total_counts"].median()) if "total_counts" in adata.obs else None,
        "median_genes_per_cell": float(adata.obs["n_genes_by_counts"].median()) if "n_genes_by_counts" in adata.obs else None,
        "pct_mt_median": float(adata.obs["pct_counts_mt"].median()) if "pct_counts_mt" in adata.obs else None,
    }
