# ML-Ready Spatial Multi-Omics Integration Framework (10x Visium)

Reproducible spatial transcriptomics analysis pipeline built using:

- 10x Genomics Visium Human Breast Cancer dataset
- Scanpy for biological analysis
- Nextflow DSL2 for modular workflow orchestration

---

## Overview

This project implements a modular and reproducible spatial transcriptomics workflow using Nextflow DSL2.

The pipeline includes:

1. Quality Control (QC)
2. Library size normalization
3. Log-transformation
4. Highly Variable Gene (HVG) selection
5. Principal Component Analysis (PCA)
6. UMAP embedding
7. Leiden clustering
8. Differential expression analysis

The workflow is modularized into separate Nextflow DSL2 processes to ensure scalability, reproducibility, and clear separation of computational stages.

---

## Dataset

- **Platform:** 10x Genomics Visium
- **Sample:** Human Breast Cancer
- **Data Type:** Spatially resolved transcriptomics count matrix (HDF5 format)
- **Input:** Filtered feature-barcode matrix (`.h5`)

The dataset preserves spatial context, enabling analysis of gene expression patterns within tissue architecture.

---

## Pipeline Architecture

main.nf
├── modules/qc.nf
├── modules/normalize.nf
└── modules/cluster.nf
├── scripts/qc.py
├── scripts/normalize.py
└── scripts/cluster.py


- Each module wraps a Python Scanpy script.
- Nextflow DSL2 handles workflow orchestration.
- Modules are independent and reusable.
- Biological logic is separated from workflow logic.

This structure enables reproducible execution and easy extension.

---

## Methodological Rationale

The workflow follows established single-cell and spatial transcriptomics preprocessing principles:

- Quality control filtering to remove low-complexity spots
- Library size normalization to account for sequencing depth variability
- Log-transformation for variance stabilization
- Highly Variable Gene (HVG) selection to focus on biologically informative features
- PCA for dimensionality reduction
- UMAP for nonlinear manifold embedding
- Leiden clustering for community detection

These steps reflect standard best practices in spatial transcriptomics analysis and ensure biological interpretability.

---

## Biological Validation

Clustering results were biologically validated using established breast cancer marker genes and spatial expression consistency.

Validation steps included:

- Verification of known epithelial and tumor-associated marker genes
- Assessment of spatial clustering coherence across tissue sections
- Differential expression analysis between Leiden clusters

The observed gene expression patterns were consistent with expected tumor microenvironment structure, supporting the biological plausibility of the clustering results.

---

## Reproducibility and Engineering Design

- Modular DSL2-based Nextflow implementation
- Separation of workflow logic and biological scripts
- Reproducible computational stages
- Structured for scalability and future nf-core adaptation
- Dependency tracking via `requirements.txt`

This project demonstrates integration of biological reasoning with computational workflow engineering.

---

## Machine Learning Integration
This pipeline is specifically engineered to transform raw spatial transcriptomics data into standardized, ML-ready tensors. It facilitates the construction of spatial neighbor graphs, providing the necessary structural input for Graph Attention Networks and Transformer-based architectures to analyze cell-cell interaction patterns

---

## How to Run

```bash
nextflow run main.nf
Author
Neha BL
Computational Biology / Data Engineering
