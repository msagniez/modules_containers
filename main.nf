process MNM_PREDICT {

    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/${params.container_owner}/modules_containers/mnm:latest"

    input:
    tuple val(meta), path(quant_file)

    output:
    tuple val(meta), path("${meta.id}/MnM_subtype-predictions.csv"), emit: subtype_predictions
    tuple val(meta), path("${meta.id}/MnM_lineage-predictions.csv"),    emit: type_predictions

    script:
    def outdir = "${task.workDir}/${meta.id}"
    """

    Rscript ${projectDir}/scripts/MnM_predict.R \\
        ${quant_file} \\
        ${outdir}
    """
}