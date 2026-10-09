"""Shared scanpy/AnnData helpers used by the preprocessing pipeline.

Distilled from the working notebook utilities: only the project-agnostic
routines are kept (QC, doublet calling, merge, marker scoring, plotting).
Interactive / spatial-format specific helpers were dropped.
"""
import os
import warnings

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scanpy as sc

warnings.filterwarnings("ignore")


# --------------------------------------------------------------------------- #
# QC and doublet detection
# --------------------------------------------------------------------------- #
def scrFilter(adata, doublet_rate=0.3, doublet_stdev=0.05,
              min_counts=2, min_cells=3, min_gene_variability_pctl=85,
              n_prin_comps=20):
    """Flag doublets with Scrublet.

    Uses the standalone ``scrublet`` package when installed and otherwise falls
    back to scanpy's built-in ``sc.external.pp.scrublet`` (scanpy >= 1.9).

    Adds ``predicted_doublets`` and ``doublet_scores`` to ``adata.obs`` and
    returns the object unchanged, so the caller decides what to drop.
    """
    try:
        import scrublet as scr
    except ImportError:
        scr = None

    if scr is not None:
        scrub = scr.Scrublet(
            adata.X,
            expected_doublet_rate=doublet_rate,
            stdev_doublet_rate=doublet_stdev,
        )
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=min_counts,
            min_cells=min_cells,
            min_gene_variability_pctl=min_gene_variability_pctl,
            n_prin_comps=n_prin_comps,
        )
        adata.obs["predicted_doublets"] = predicted_doublets
        adata.obs["doublet_scores"] = doublet_scores
        return adata

    # scanpy fallback
    sc.external.pp.scrublet(
        adata,
        expected_doublet_rate=doublet_rate,
        use_approx_neighbors=False,
        n_prin_comps=min(n_prin_comps, min(adata.n_obs, adata.n_vars) - 1) or None,
        verbose=False,
    )
    adata.obs["doublet_scores"] = adata.obs["doublet_score"]
    adata.obs["predicted_doublets"] = adata.obs["predicted_doublet"].astype(bool)
    return adata


def qc_pip(adata, logfile, keep=False,
           min_counts=500, min_cells=10,
           gene_low=250, gene_high=5000, pct_mt=20,
           log10_umi_low=0.8, total_counts_max=25000,
           doublet_rate=0.3):
    """Per-sample QC exactly as described in the Methods.

    Steps
    -----
    1. mitochondrial percentage (genes starting with ``MT-``),
    2. quantile-based gene-count filter plus the fixed thresholds
       (250 < n_genes < 5000, pct_counts_mt < 20),
    3. Scrublet doublet removal.

    Parameters
    ----------
    keep : bool
        If True, return the *unfiltered* object with QC columns attached
        (``predicted_doublets``, ``doublet_scores``, ``qc_keep``) instead of a
        subset, so several samples can be filtered independently and merged
        afterwards.

    Returns
    -------
    AnnData
    """
    with open(logfile, "w") as fh:
        print(adata, file=fh)

        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        adata.obs["log10genes_by_umi"] = (
            np.log10(adata.obs["n_genes_by_counts"]) / np.log10(adata.obs["total_counts"])
        )

        raw = adata.copy() if keep else None

        qtl_5 = np.quantile(adata.obs.n_genes_by_counts, q=0.05)
        qtl_95 = np.quantile(adata.obs.n_genes_by_counts, q=0.95)
        sc.pp.filter_cells(adata, min_counts=min_counts)
        sc.pp.filter_cells(adata, min_genes=qtl_5)
        sc.pp.filter_cells(adata, max_genes=qtl_95)
        sc.pp.filter_genes(adata, min_cells=min_cells)

        adata = adata[
            (adata.obs.pct_counts_mt < pct_mt)
            & (adata.obs.log10genes_by_umi > log10_umi_low)
            & (adata.obs.n_genes_by_counts > gene_low)
            & (adata.obs.n_genes_by_counts < gene_high)
            & (adata.obs.total_counts < total_counts_max),
            :
        ].copy()

        print(adata, file=fh)

        if adata.n_obs == 0:
            # Fail loudly and early: an empty object would otherwise surface as
            # a confusing sklearn error from deep inside Scrublet.
            raise ValueError(
                "QC removed every cell. Thresholds "
                f"(gene_low={gene_low}, gene_high={gene_high}, pct_mt={pct_mt}, "
                f"log10_umi_low={log10_umi_low}, total_counts_max={total_counts_max}) "
                "are incompatible with this sample - check that the input holds raw "
                "counts and that gene symbols match the expected species prefix."
            )

        adata = scrFilter(adata, doublet_rate=doublet_rate)

        if keep:
            submask = raw.obs_names.isin(adata.obs_names)
            for col in ("predicted_doublets", "doublet_scores", "qc_keep"):
                raw.obs[col] = np.nan
            raw.obs.loc[submask, "predicted_doublets"] = adata.obs["predicted_doublets"]
            raw.obs.loc[submask, "doublet_scores"] = adata.obs["doublet_scores"]

            adata = adata[adata.obs["predicted_doublets"] == False, :].copy()  # noqa: E712
            submask = raw.obs_names.isin(adata.obs_names)
            raw.obs["qc_keep"] = 0
            raw.obs.loc[submask, "qc_keep"] = 1

            for col in ("n_genes_by_counts", "log10genes_by_umi", "pct_counts_mt", "total_counts"):
                print(np.quantile(adata.obs[col], q=[0, 1]), file=fh)
            print(raw, file=fh)
            return raw

        adata = adata[adata.obs["predicted_doublets"] == False, :].copy()  # noqa: E712
        for col in ("n_genes_by_counts", "log10genes_by_umi", "pct_counts_mt", "total_counts"):
            print(np.quantile(adata.obs[col], q=[0, 1]), file=fh)
        print(adata, file=fh)
        return adata


def merge_h5ad_to_adata(file_paths):
    """Concatenate several .h5ad files, prefixing barcodes with the file stem."""
    adatas = []
    for file_path in file_paths:
        adata = sc.read_h5ad(file_path)
        stem = os.path.splitext(os.path.basename(file_path))[0]
        adata.obs_names = [f"{stem}_{bc}" for bc in adata.obs_names]
        adatas.append(adata)

    merged = adatas[0].concatenate(adatas[1:], index_unique=None)
    return merged


# --------------------------------------------------------------------------- #
# Metadata helpers
# --------------------------------------------------------------------------- #
def add_level(adata, obsCol="sample", col_order=None):
    """Turn ``adata.obs[obsCol]`` into an ordered categorical."""
    if col_order is None:
        warnings.warn("col_order is None; levels left unchanged.")
        return adata

    observed = set(adata.obs[obsCol])
    wanted = set(col_order)
    intersect = list(wanted.intersection(observed))
    extra_obs = list(observed.difference(wanted))
    extra_wanted = list(wanted.difference(observed))

    if extra_wanted:
        warnings.warn("Some col_order entries are absent from the data and were ignored.")
    if extra_obs:
        warnings.warn("Some data entries are absent from col_order and were appended.")

    new_order = intersect + extra_obs
    adata.obs[obsCol] = pd.Categorical(
        adata.obs[obsCol], categories=new_order, ordered=True
    )
    adata.obs.sort_values(obsCol, inplace=True)
    return adata


def add_metadata(adata, meta_df, lefCol="sample", rigCol="SampleID", col_order=None):
    """Left-join a sample sheet onto ``adata.obs``."""
    if rigCol not in meta_df.columns:
        raise ValueError(f"meta_df must contain a '{rigCol}' column.")

    adata.obs = adata.obs.merge(meta_df, left_on=lefCol, right_on=rigCol, how="left")
    if lefCol != rigCol:
        adata.obs.drop(rigCol, axis=1, inplace=True)

    if col_order is not None:
        adata = add_level(adata, obsCol=lefCol, col_order=col_order)
    return adata


# --------------------------------------------------------------------------- #
# Clustering / annotation pipeline
# --------------------------------------------------------------------------- #
def sc_ana(adata, mharmony=0, preprocess=True, leidenrun=True,
           markers=None, subsample=True, n_cells=100000,
           n_top_genes=3000, batch_key="batch"):
    """Normalise, embed, cluster and (optionally) score-marker annotate.

    Mirrors the Methods: log-normalisation, HVG selection, regression of
    unwanted variation, PCA, Harmony batch correction, neighbour graph, UMAP
    and Leiden clustering.
    """
    if subsample and adata.n_obs > n_cells:
        sc.pp.subsample(adata, n_obs=n_cells)

    if preprocess:
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.layers["data"] = adata.X.copy()
        adata.raw = adata

        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            subset=True,
            layer="counts",
            flavor="seurat_v3",
        )
        sc.pp.regress_out(adata, ["n_counts", "pct_counts_mt"])
        sc.pp.scale(adata, max_value=10)
        adata.layers["scale_data"] = adata.X.copy()
        sc.tl.pca(adata, svd_solver="arpack")

    if mharmony > 0:
        # Harmony embeddings, made robust to the harmonypy API change:
        #   1.x returns Z_corr shaped (n_components, n_cells)
        #   2.x returns Z_corr shaped (n_cells, n_components)
        # scanpy's own harmony_integrate still assumes 1.x and therefore breaks
        # under harmonypy 2.0, so orient the matrix by n_obs instead.
        import harmonypy as harmony

        out = harmony.run_harmony(
            data_mat=np.asarray(adata.obsm["X_pca"]),
            meta_data=adata.obs,
            vars_use=[batch_key],
            max_iter_harmony=mharmony,
        )
        Z = np.asarray(out.Z_corr)
        if Z.shape[0] == adata.n_obs:
            corrected = Z
        elif Z.shape[-1] == adata.n_obs:
            corrected = Z.reshape(Z.shape[0], -1).T
        else:
            raise ValueError(
                f"Harmony returned Z_corr with shape {Z.shape}, which matches "
                f"neither n_obs={adata.n_obs} nor n_comps={adata.obsm['X_pca'].shape[1]}. "
                "Check the installed harmonypy version."
            )
        adata.obsm["X_pca_harmony"] = corrected
        sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    else:
        sc.pp.neighbors(adata, use_rep="X_pca")

    sc.tl.umap(adata)

    if leidenrun:
        sc.tl.leiden(adata)

    if markers is not None:
        for cell_type, genes in markers.items():
            sc.tl.score_genes(adata, gene_list=genes, score_name=f"{cell_type}_score")

    return adata.copy()


def assign_by_max_score(adata, score_cols, out_col, strip_suffix=None):
    """Label each cell by its highest-scoring signature.

    Scores are min-max normalised first so that signatures of different sizes
    stay comparable, then the argmax column name becomes the label.

    Parameters
    ----------
    score_cols : list[str]
        Columns in ``adata.obs`` holding the per-signature scores.
    strip_suffix : str, optional
        Suffix to remove from the column names before writing the label
        (e.g. ``"_major_type"``).
    """
    df = adata.obs[score_cols].copy()
    df = (df - df.min()) / (df.max() - df.min())
    adata.obs[out_col] = df.idxmax(axis=1)
    if strip_suffix:
        adata.obs[out_col] = adata.obs[out_col].str.replace(strip_suffix, "", regex=False)
    return adata


def degTopN(adata, cluName="cluster", CToI=None, n=4, padj=0.05, lfc=0.25):
    """Collect per-cluster marker genes from a ``rank_genes_groups`` result."""
    if CToI is None:
        col = adata.obs[cluName]
        cats = list(col.cat.categories) if hasattr(col, "cat") else list(col.unique())
    else:
        cats = list(CToI)

    full = pd.DataFrame()
    for c in cats:
        df = sc.get.rank_genes_groups_df(adata, group=str(c))
        df["comparison"] = np.repeat(f"{c} vs. rest", len(df))
        full = pd.concat([full, df])

    markers = {}
    for c in cats:
        cid = f"{c} vs. rest"
        sub = full.loc[full.comparison == cid]
        sub = sub.sort_values("logfoldchanges", ascending=False)
        sub = sub[sub.logfoldchanges > lfc]
        sub = sub[sub.pvals_adj < padj]
        markers[c] = sub.names.tolist()[:n]

    return full, markers


def sub_score(adata, submask, gene_listV, score_name):
    """Run ``sc.tl.score_genes`` on a subset and write results back as NaN-padded."""
    cell_adata = adata[submask].copy()
    sc.tl.score_genes(cell_adata, gene_list=gene_listV, score_name=score_name)
    scored = cell_adata.obs[[score_name]]

    adata.obs[score_name] = np.nan
    adata.obs.loc[submask, [score_name]] = scored
    return adata


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def color_pick(category_list):
    """Deterministic colour-blind-safe palette mapped onto category names."""
    np.random.seed(0)
    palette = sns.color_palette("colorblind", len(category_list))
    hex_colors = [mcolors.to_hex(c) for c in palette]
    return {cat: hex_colors[i] for i, cat in enumerate(category_list)}


def split_umap(adata, color, split_by, ncol=2, colorplat=None, nrow=None,
               leg_pos="right margin", basis="umap", **kwargs):
    """One UMAP panel per level of ``split_by``."""
    col = adata.obs[split_by]
    categories = list(col.cat.categories) if hasattr(col, "cat") else list(col.unique())
    if nrow is None:
        nrow = int(np.ceil(len(categories) / ncol))

    fig, axs = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4 * nrow))
    axs = np.array(axs).flatten()
    for i, cat in enumerate(categories):
        sub = adata[adata.obs[split_by] == cat].copy()
        sc.pl.embedding(sub, basis=basis, palette=colorplat, color=color,
                        ax=axs[i], show=False, title=str(cat),
                        legend_loc=leg_pos, **kwargs)
    plt.tight_layout()
    return fig


def plot_umap_metrics(adata, metrics, mmax=None, mmin=None, splitby=None,
                      cpalette=None, ncols=3, figsize=(4, 4), save_fig=False,
                      fig_path="umap_metrics.png", use_rawtf=True,
                      ptsize_scale=1, leg_pos=None):
    """Grid of UMAPs for obs metrics and gene expression, optionally split."""
    if splitby and splitby in adata.obs:
        for category in adata.obs[splitby].unique():
            sub = adata[adata.obs[splitby] == category].copy()
            _plot_single_umap(sub, mmax, mmin, metrics, ncols, figsize,
                              f"{splitby} = {category}", save_fig, fig_path,
                              use_raw=use_rawtf, size_scale=ptsize_scale,
                              legpos=leg_pos, ipalette=cpalette)
    else:
        _plot_single_umap(adata, mmax, mmin, metrics, ncols, figsize, None,
                          save_fig, fig_path, use_raw=use_rawtf,
                          size_scale=ptsize_scale, legpos=leg_pos, ipalette=cpalette)


def _plot_single_umap(adata, mmax, mmin, metrics, ncols, figsize, title,
                      save_fig, fig_path, use_raw, size_scale=1,
                      legpos="best", ipalette=None):
    nrows = int(len(metrics) / ncols) + (len(metrics) % ncols > 0)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=(ncols * figsize[0], nrows * figsize[1]))
    axes = np.array(axes).ravel() if nrows * ncols > 1 else [axes]

    psize = 0.15  # tuned for ~120k cells; scale up for smaller objects
    for i, metric in enumerate(metrics):
        ax = axes[i]
        if metric in adata.obs:
            sc.pl.umap(adata, color=metric, size=size_scale * psize, palette=ipalette,
                       vmax=mmax, vmin=mmin, ax=ax, show=False, title=metric,
                       legend_loc=legpos)
        elif adata.raw is not None and metric in adata.raw.var_names:
            sc.pl.umap(adata, color=metric, size=size_scale * psize, cmap=ipalette,
                       vmax=mmax, vmin=mmin, ax=ax, show=False, title=metric,
                       use_raw=use_raw)
        else:
            print(f"{metric} not found in adata.obs or adata.raw.var_names.")

    plt.suptitle(title or "", fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if save_fig:
        plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.show()


def plot_vln(adata, keys=None, prefix="", stripplot=False, ncols_nrows=(1, 3),
             rotation_angle=45, groupby="sample"):
    """Violin QC plots for a list of metrics, saved next to ``prefix``."""
    if keys is None:
        return None

    ncols, nrows = ncols_nrows
    fig, axs = plt.subplots(nrows, ncols, figsize=(16, 5 * nrows),
                            gridspec_kw={"wspace": 0.25, "hspace": 0.4})
    axs = np.array(axs).flatten()

    for i, key in enumerate(keys):
        ax = axs[i]
        sc.pl.violin(adata, keys=key, groupby=groupby, log=False,
                     stripplot=stripplot, rotation=rotation_angle,
                     ax=ax, show=False)
        for label in ax.get_xticklabels():
            label.set_ha("right")
            label.set_va("top")
            label.set_rotation(rotation_angle)
        ax.set_ylabel(key if i % ncols == 0 else "")
        ax.set_xlabel("")

    plt.tight_layout()
    sns.despine()
    fig.savefig(prefix + "_vlnplot.png", dpi=300)
    return fig
