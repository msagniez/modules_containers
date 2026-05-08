process ATTENTIONAML_PREDICT {

    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/${params.container_owner}/modules_containers/attentionaml:latest"

    input:
    tuple val(meta), path(quant_file)

    output:
    tuple val(meta), path("${meta.id}_predictions.csv"),    emit: prediction_csv

    script:
    def outdir = "${task.workDir}/${meta.id}"
    """

    python3 /opt/AttentionAML/AttentionAML_predict.py \\
        ${quant_file} \\
        ${outdir}
    """
}