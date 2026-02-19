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


def run_normalization(input_file, output_file, target_sum=1e4, n_top_genes=3000):
    logging.info("Loading QC-filtered AnnData object...")
    adata = sc.read(input_file)

    logging.info(f"Normalizing counts to target sum = {target_sum}...")
    sc.pp.normalize_total(adata, target_sum=target_sum)

    logging.info("Applying log1p transformation...")
    sc.pp.log1p(adata)

    logging.info("Storing raw normalized data...")
    adata.raw = adata

    logging.info(f"Selecting top {n_top_genes} highly variable genes (Seurat flavor)...")
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        n_top_genes=n_top_genes
    )

    adata = adata[:, adata.var.highly_variable]

    logging.info(f"Remaining genes after HVG selection: {adata.n_vars}")

    logging.info("Saving normalized AnnData object...")
    adata.write(output_file)

    logging.info("Normalization module completed successfully.")


if __name__ == "__main__":
    configure_logging()

    parser = argparse.ArgumentParser(description="Normalization and HVG selection module.")
    parser.add_argument("--input", required=True, help="Path to QC-filtered .h5ad file")
    parser.add_argument("--output", required=True, help="Output normalized .h5ad file")
    parser.add_argument("--target_sum", type=float, default=1e4)
    parser.add_argument("--n_top_genes", type=int, default=3000)

    args = parser.parse_args()

    validate_input(args.input)
    run_normalization(
        args.input,
        args.output,
        args.target_sum,
        args.n_top_genes
    )
