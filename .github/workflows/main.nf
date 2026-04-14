process MAYLIS_PREDICT {

    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/${params.container_owner}/modules_containers/maylis:latest"

    input:
    tuple val(meta), path(quant_file)

    output:
    tuple val(meta), path("${meta.id}_predictions.png"), emit: prediction_figure
    tuple val(meta), path("${meta.id}_predictions.tsv"),    emit: prediction_tsv

    script:
    def outdir = "${task.workDir}/${meta.id}"
    """

    Rscript /opt/scripts/Maylis_predict.R \\
        ${quant_file} \\
        ${outdir}
    """
}