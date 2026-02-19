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
        logging.error(f"Input path does not exist: {path}")
        sys.exit(1)


def run_qc(input_dir, output_path, min_genes=500):
    logging.info("Loading Visium dataset...")
    adata = sc.read_visium(path=input_dir)

    logging.info("Ensuring unique gene names...")
    adata.var_names_make_unique()

    logging.info("Calculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    logging.info(f"Filtering spots with < {min_genes} detected genes...")
    adata = adata[adata.obs["n_genes_by_counts"] > min_genes, :]

    logging.info(f"Remaining spots after filtering: {adata.n_obs}")

    logging.info("Saving filtered AnnData object...")
    adata.write(output_path)

    logging.info("QC completed successfully.")


if __name__ == "__main__":
    configure_logging()

    parser = argparse.ArgumentParser(description="QC module for 10x Visium spatial transcriptomics data.")
    parser.add_argument("--input", required=True, help="Path to Visium directory")
    parser.add_argument("--output", required=True, help="Output filtered .h5ad file")
    parser.add_argument("--min_genes", type=int, default=500)

    args = parser.parse_args()

    validate_input(args.input)
    run_qc(args.input, args.output, args.min_genes)
