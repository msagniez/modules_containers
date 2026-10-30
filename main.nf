process MNM_PREDICT {

    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/msagniez/modules-containers/mnm:latest"

    input:
    tuple val(meta), path(quant_file), path(outdir)

    output:
    tuple val(meta), path("${meta.id}/MnM_subtype-predictions.csv"),    emit: subtype_predictions
    tuple val(meta), path("${meta.id}/MnM_lineage-predictions.csv"),    emit: type_predictions
    path "versions.yml",                                                emit: versions

    script:
    def outdir = meta.id
    """
    mkdir -p ${outdir}

    Rscript /opt/scripts/MnM_predict.R \\
        ${quant_file} \\
        ${outdir}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mnm: \$(Rscript -e 'cat(as.character(packageVersion("MnM")))')
    END_VERSIONS
    """
}
