process GENE_SUMMARIZE {
    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/${params.container_owner}/modules_containers/tximport:latest"

    input:
    tuple val(meta), path(quant_files)
    path(reference)                      // either a GTF file or a tx2gene file
    val(use_gene_names)                  // true/false: use --gene-names flag (default: false)
    val(ignore_tx_version)               // true/false: use --ignore-tx-version flag (default: false)
    val(type)                            // type of quantifier used: salmon, oarfish, kallisto, rsem, etc. (default: oarfish)

    output:
    tuple val(meta), path("gene_counts.tsv"), emit: counts
    tuple val(meta), path("gene_tpm.tsv"),    emit: tpm
    path("versions.yml"),                     emit: versions

    script:
    def gene_names     = use_gene_names    ? "--gene-names"        : ""
    def ignore_version = ignore_tx_version ? "--ignore-tx-version" : ""
    def tx2gene_flag   = (reference.name.endsWith(".txt") || reference.name.endsWith(".tsv")) ? "--tx2gene" : ""
    def type_flag      = type ? "--type ${type}" : ""
    """
    Rscript /opt/scripts/gene_summarize.R \\
        ${gene_names} \\
        ${ignore_version} \\
        ${type_flag} \\
        ${tx2gene_flag} \\
        ${reference} \\
        ${quant_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/^R version //; s/ .*\$//')
        tximport: \$(Rscript -e 'cat(as.character(packageVersion("tximport")))')
    END_VERSIONS
    """

    stub:
    """
    touch gene_counts.tsv
    touch gene_tpm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        tximport: stub
    END_VERSIONS
    """
}