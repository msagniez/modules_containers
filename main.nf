process FUSILLI_PREDICT {

    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/${params.container_owner}/modules_containers/fusilli:latest"

    input:
    tuple val(meta), path(quant_file)

    output:
    tuple val(meta), path("${meta.id}_fusilli.txt"),    emit: fusilli_fusion_list

    script:
    def outdir = "${task.workDir}/${meta.id}"
    """

    fusilli \\
      -p ${quant_file} \\
      --nfm \\
      -o ${outdir}
    """
}