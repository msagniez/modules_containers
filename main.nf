process ALLCATCHR {
    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/${params.container_owner}/modules_containers/allcatchr:latest"

    input:
    tuple val(meta), path(quant_dir)

    output:
    tuple val(meta), path("${meta.id}/*_lineage-prediction.tsv"),  emit: lineage
    tuple val(meta), path("${meta.id}/*_ball-prediction.tsv"),     emit: ball
    tuple val(meta), path("${meta.id}/*_ball-prediction_plot*"),   emit: ball_plot
    tuple val(meta), path("${meta.id}/*_bcrabl1-prediction.tsv"),  emit: bcrabl1
    tuple val(meta), path("${meta.id}/*_tall-prediction.tsv"),     emit: tall
    tuple val(meta), path("${meta.id}/*_tall-prediction_plot*"),   emit: tall_plot

    script:
    def prefix = meta.id
    """
    mkdir -p ${prefix}

    Rscript /opt/scripts/ALLcatchR_predict.R \\
        ${quant_dir} \\
        ${prefix}
    """
}