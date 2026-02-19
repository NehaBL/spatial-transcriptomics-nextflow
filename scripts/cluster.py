import scanpy as sc
import argparse
import logging
import os
import sys


def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )


def validate_input(path):
    if not os.path.exists(path):
        logging.error(f"Input file does not exist: {path}")
        sys.exit(1)


def run_clustering(input_file, output_file, resolution=0.5):
    logging.info("Loading normalized AnnData object...")
    adata = sc.read(input_file)

    logging.info("Running PCA...")
    sc.tl.pca(adata, svd_solver="arpack")

    logging.info("Computing neighborhood graph...")
    sc.pp.neighbors(adata)

    logging.info("Computing UMAP embedding...")
    sc.tl.umap(adata)

    logging.info(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution)

    logging.info("Performing differential expression (Wilcoxon test)...")
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon"
    )

    logging.info(f"Identified {adata.obs['leiden'].nunique()} clusters.")

    logging.info("Saving clustered AnnData object...")
    adata.write(output_file)

    logging.info("Clustering module completed successfully.")


if __name__ == "__main__":
    configure_logging()

    parser = argparse.ArgumentParser(description="Clustering and differential expression module.")
    parser.add_argument("--input", required=True, help="Path to normalized .h5ad file")
    parser.add_argument("--output", required=True, help="Output clustered .h5ad file")
    parser.add_argument("--resolution", type=float, default=0.5)

    args = parser.parse_args()

    validate_input(args.input)
    run_clustering(
        args.input,
        args.output,
        args.resolution
    )
