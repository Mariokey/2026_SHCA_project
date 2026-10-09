"""Per-sample QC -> merge -> annotate pipeline (scanpy).

Distilled from the original working notebook. The notebook itself operated on
an internal Stereo-seq dataset; what is kept here is the generic, reusable
flow that matches the Methods section "Preprocessing of single-cell
RNA-sequencing data":

    raw counts (Cell Ranger output, per sample)
      -> QC filter   : 250 < n_genes < 5000, pct_counts_mt < 20
      -> doublets    : Scrublet (sample-specific)
      -> merge samples
      -> normalise / HVG / PCA / Harmony / UMAP / Leiden
      -> marker-based annotation by max module score

Interactive helpers from the original notebook (dash point-selection, Stereo-seq
`.gef` reading, plotly previews) were deliberately not carried over.

Usage
-----
    python preprocess_pipeline.py \
        --h5ad-dir ./data/per_sample_h5ad \
        --meta ./data/sample_sheet.csv \
        --markers ./ref/marker_signatures.tsv \
        --out ./data/merged_annotated.h5ad

Expected inputs
---------------
* ``--h5ad-dir``   one .h5ad per sample, raw counts, ``obs['sample']`` optional
* ``--meta``       sample sheet with a ``SampleID`` column (extra columns are
                   joined into ``adata.obs``)
* ``--markers``    TSV: first column ``celltype``, second column ``gene``
                   (long format; the per-cell-type gene lists feed
                   ``sc.tl.score_genes``)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "utils"))
import sc_utils as uti  # noqa: E402


def collect_inputs(h5ad_dir):
    return sorted(
        os.path.join(h5ad_dir, f)
        for f in os.listdir(h5ad_dir)
        if f.endswith(".h5ad")
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad-dir", required=True, help="directory of per-sample raw .h5ad files")
    ap.add_argument("--meta", default=None, help="sample sheet CSV with a SampleID column")
    ap.add_argument("--markers", default=None,
                    help="long-format TSV (celltype, gene) of marker signatures")
    ap.add_argument("--out", default="./merged_annotated.h5ad", help="output .h5ad path")
    ap.add_argument("--qc-dir", default="./qc_logs", help="directory for per-sample QC logs")
    ap.add_argument("--harmony", type=int, default=10,
                    help="max Harmony iterations (0 disables batch correction)")
    ap.add_argument("--batch-key", default="sample", help="obs column used as batch")
    ap.add_argument("--n-top-genes", type=int, default=3000, help="HVG count")
    # QC thresholds default to the values stated in the Methods.
    ap.add_argument("--min-genes", type=int, default=250, help="QC: lower gene-count bound")
    ap.add_argument("--max-genes", type=int, default=5000, help="QC: upper gene-count bound")
    ap.add_argument("--max-pct-mt", type=float, default=20.0, help="QC: max mitochondrial percent")
    ap.add_argument("--min-log10-umi", type=float, default=0.8,
                    help="QC: lower bound on log10(genes)/log10(UMI)")
    ap.add_argument("--max-total-counts", type=int, default=25000, help="QC: max total counts")
    ap.add_argument("--doublet-rate", type=float, default=0.3,
                    help="Scrublet expected doublet rate")
    args = ap.parse_args()

    os.makedirs(args.qc_dir, exist_ok=True)
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150)

    # ------------------------------------------------------------------ QC
    files = collect_inputs(args.h5ad_dir)
    if not files:
        raise SystemExit(f"no .h5ad files under {args.h5ad_dir}")
    print(f"[qc] {len(files)} sample files")

    qc_list = []
    for f in files:
        stem = os.path.splitext(os.path.basename(f))[0]
        adata = sc.read_h5ad(f)
        if "sample" not in adata.obs:
            adata.obs["sample"] = stem
        log = os.path.join(args.qc_dir, f"{stem}.qc.log")
        n_before = adata.n_obs
        kept = uti.qc_pip(
            adata,
            logfile=log,
            keep=False,
            gene_low=args.min_genes,
            gene_high=args.max_genes,
            pct_mt=args.max_pct_mt,
            log10_umi_low=args.min_log10_umi,
            total_counts_max=args.max_total_counts,
            doublet_rate=args.doublet_rate,
        )
        print(f"[qc] {stem}: {n_before} -> {kept.n_obs} cells")
        qc_list.append(kept)

    # anndata.concat is the current API; AnnData.concatenate was removed in
    # recent anndata releases.
    import anndata as ad

    merged = ad.concat(
        qc_list,
        join="outer",
        label="sample",
        keys=[os.path.splitext(os.path.basename(f))[0] for f in files],
        index_unique="-",
        merge="same",
    )
    merged.obs_names_make_unique()
    print(f"[merge] {merged.n_obs} cells, {merged.n_vars} genes")

    # ---------------------------------------------------------- sample sheet
    if args.meta:
        meta = pd.read_csv(args.meta)
        merged = uti.add_metadata(merged, meta, lefCol="sample", rigCol="SampleID")
        print(f"[meta] joined columns: {[c for c in meta.columns if c != 'SampleID']}")

    # ------------------------------------------------------- process/cluster
    processed = uti.sc_ana(merged, mharmony=args.harmony, preprocess=True,
                           leidenrun=True, batch_key=args.batch_key,
                           n_top_genes=args.n_top_genes, subsample=False)

    # ----------------------------------------------------------- annotation
    if args.markers:
        sig = pd.read_csv(args.markers, sep="\t")
        sig.columns = [c.strip().lower() for c in sig.columns]
        markers = (sig.groupby("celltype")["gene"]
                      .apply(lambda g: sorted(set(g))).to_dict())
        suffix = "_major_type"
        for cell_type, genes in markers.items():
            sc.tl.score_genes(processed, gene_list=genes,
                              score_name=f"{cell_type}{suffix}")
        score_cols = [f"{ct}{suffix}" for ct in markers]
        processed = uti.assign_by_max_score(processed, score_cols,
                                            out_col="major_type",
                                            strip_suffix=suffix)
        print(f"[anno] assigned {processed.obs['major_type'].nunique()} labels "
              f"by max signature score")

    processed.write(args.out)
    print(f"[out] wrote {args.out}")


if __name__ == "__main__":
    main()
