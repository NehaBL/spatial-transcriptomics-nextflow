nextflow.enable.dsl=2

params.input   = null
params.outdir  = "results"
params.min_genes = 500
params.target_sum = 1e4
params.n_top_genes = 3000
params.resolution = 0.5

if (!params.input) {
    error "Please provide --input pointing to Visium directory"
}

workflow {

    Channel
        .fromPath(params.input)
        .set { visium_input }

    qc_output = QC(visium_input)

    norm_output = NORMALIZE(qc_output)

    CLUSTER(norm_output)
}

process QC {

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path visium_dir

    output:
    path "qc_filtered.h5ad"

    script:
    """
    python scripts/qc.py \
        --input ${visium_dir} \
        --output qc_filtered.h5ad \
        --min_genes ${params.min_genes}
    """
}

process NORMALIZE {

    publishDir "${params.outdir}/normalize", mode: 'copy'

    input:
    path qc_file

    output:
    path "normalized.h5ad"

    script:
    """
    python scripts/normalize.py \
        --input ${qc_file} \
        --output normalized.h5ad \
        --target_sum ${params.target_sum} \
        --n_top_genes ${params.n_top_genes}
    """
}

process CLUSTER {

    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path norm_file

    output:
    path "clustered.h5ad"

    script:
    """
    python scripts/cluster.py \
        --input ${norm_file} \
        --output clustered.h5ad \
        --resolution ${params.resolution}
    """
}
