process OARFISH_QUANTIFY_READ {
    // TODO : SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/msagniez/oarfish:latest'
    // TODO : SET LEVEL OF RESSOURCES
    tag "$meta.id"
    label 'process_medium'
    label 'process_medium_high_cpu'
    label 'process_low_memory'
    label 'process_low_time'

    input:
    tuple val(meta),
        val(biotype), // ont-cdna or ont-drna or pac-bio or pac-bio-hifi
        path(reads),  // input reads as fastq or fastq.gz
        path(ref)     // reference transcriptomeas fasta/fa or fasta.gz/fa.gz (transcripts/isoforms sequences, gencode, ensembl, etc.)

    output:
    tuple val(meta),
        path("*.quant"),
        emit: quant
    tuple val(meta),
        path("*.tsv"),
        emit: ambig
    tuple val(meta),
        path("*.json"),
        emit: meta_info
    path "versions.yml",
        emit: versions

    script:
    def args = task.ext?.args ?: ''
    def prefix = task.ext?.prefix ?: "${meta.id}"
    def threads = task.cpus
    """
    oarfish \\
        -j ${threads} \\
        ${args} \\
        --reads ${reads} \\
        --annotated ${ref} \\
        --seq-tech ${biotype} \\
        -o ${prefix} \\
        --filter-group no-filters \\
        --model-coverage

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(echo \$(oarfish --version 2>&1) | sed 's/^.*oarfish //; s/Using.*\$//')
    END_VERSIONS
    """
}

process OARFISH_QUANTIFY_ALIGNMENT {
    // TODO : SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/msagniez/oarfish:latest'
    // TODO : SET LEVEL OF RESSOURCES
    tag "$meta.id"
    label 'process_medium'
    label 'process_medium_high_cpu'
    label 'process_low_memory'
    label 'process_low_time'

    input:
    tuple val(meta),
        path(bam)     // input mapping to reference transcriptome as bam or bam.gz

    output:
    tuple val(meta),
        path("*.quant"),
        emit: quant
    tuple val(meta),
        path("*.tsv"),
        emit: ambig
    tuple val(meta),
        path("*.json"),
        emit: meta_info
    path "versions.yml",
        emit: versions

    script:
    def args = task.ext?.args ?: ''
    def prefix = task.ext?.prefix ?: "${meta.id}"
    def threads = task.cpus
    """
    oarfish \\
        -j ${threads} \\
        $args \\
        -a ${bam} \\
        -o ${prefix} \\
        --filter-group no-filters \\
        --model-coverage

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(echo \$(oarfish --version 2>&1) | sed 's/^.*oarfish //; s/Using.*\$//')
    END_VERSIONS
    """
}

process OARFISH_QUANTIFY_GENOMEREADPROJECTION {
    // TODO : SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/msagniez/oarfish:latest'
    // TODO : SET LEVEL OF RESSOURCES
    tag "$meta.id"
    label 'process_medium'
    label 'process_medium_high_cpu'
    label 'process_low_memory'
    label 'process_low_time'

    input:
    tuple val(meta),
        val(biotype), // ont-cdna or ont-drna or pac-bio or pac-bio-hifi
        path(reads),  // input reads as fastq or fastq.gz
        path(ref),    // genome reference as fasta or fasta.gz (hg19, hg38, mm10, etc.)
        path(gtf)     // transcriptome annotation as gtf or gtf.gz (gencode, ensembl, etc.)

    output:
    tuple val(meta),
        path("*.quant"),
        emit: quant
    tuple val(meta),
        path("*.tsv"),
        emit: ambig
    tuple val(meta),
        path("*.json"),
        emit: meta_info
    path "versions.yml",
        emit: versions

    script:
    def args = task.ext?.args ?: ''
    def prefix = task.ext?.prefix ?: "${meta.id}"
    def threads = task.cpus
    """
    oarfish \\
        -j ${threads} \\
        $args \\
        --reads ${reads} \\
        --genome ${ref} \\
        --annotation ${gtf} \\
        --seq-tech ${biotype} \\
        -o ${prefix} \\
        --filter-group no-filters

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(echo \$(oarfish --version 2>&1) | sed 's/^.*oarfish //; s/Using.*\$//')
    END_VERSIONS
    """
}

process OARFISH_QUANTIFY_GENOMEALIGNMENTPROJECTION {
    // TODO : SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/msagniez/oarfish:latest'
    // TODO : SET LEVEL OF RESSOURCES
    tag "$meta.id"
    label 'process_medium'
    label 'process_medium_high_cpu'
    label 'process_low_memory'
    label 'process_low_time'

    input:
    tuple val(meta),
        val(biotype),   // ont-cdna or ont-drna or pac-bio or pac-bio-hifi
        path(bam),      // spliced-aligned reads as bam or bam.gz
        path(ref),      // reference genome used for alignment as fasta or fasta.gz (hg19, hg38, mm10, etc.)
        path(gtf)       // transcriptome annotation as gtf or gtf.gz (gencode, ensembl, etc.)

    output:
    tuple val(meta),
        path("*.quant"),
        emit: quant
    tuple val(meta),
        path("*.tsv"),
        emit: ambig
    tuple val(meta),
        path("*.json"),
        emit: meta_info
    path "versions.yml",
        emit: versions

    script:
    def args = task.ext?.args ?: ''
    def prefix = task.ext?.prefix ?: "${meta.id}"
    def threads = task.cpus
    """
    oarfish \\
        -j ${threads} \\
        $args \\
        --genome-alignments ${bam} \\
        --annotation ${gtf} \\
        --genome-fasta ${ref} \\
        -o ${prefix} \\
        --filter-group no-filters

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(echo \$(oarfish --version 2>&1) | sed 's/^.*oarfish //; s/Using.*\$//')
    END_VERSIONS
    """
}