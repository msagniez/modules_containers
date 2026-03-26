process MDALL_PREDICT {

    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/${params.container_owner}/modules_containers/mdall:latest"

    input:
    tuple val(meta), path(quant_dir)
    val feature_ns

    output:
    tuple val(meta), path("${meta.id}_umap.pdf"),                emit: umap
    tuple val(meta), path("${meta.id}_expression-boxplot.pdf"),  emit: boxplot
    tuple val(meta), path("${meta.id}_prediction-heatmap.pdf"),  emit: heatmap
    path "MD-ALL_predictions.csv",                               emit: predictions
    path "failed_samples.txt", optional: true,                   emit: failed

    script:
    def prefix = meta.id
    def feature_arg = feature_ns ? "${feature_ns}" : ""
    """
    Rscript /opt/scripts/MDALL_predict.R \\
        ${quant_dir} \\
        ${prefix} \\
        ${feature_arg}
    """
}