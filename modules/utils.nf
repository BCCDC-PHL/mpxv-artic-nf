process articDownloadScheme{
    tag params.schemeRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scheme", mode: "copy"

    output:
    path "scheme/**/${params.schemeVersion}/*.reference.fasta" , emit: reffasta
    path "scheme/**/${params.schemeVersion}/${params.scheme}.bed" , emit: bed
    path "scheme" , emit: scheme

    script:
    """
    git clone ${params.schemeRepoURL} scheme
    """
}

process get_bed_ref {
    label 'process_single'

    container 'jitesoft/alpine:3.20.2'

    input:
        path scheme_dir
        val scheme_name
        val scheme_version
    output:
        path "scheme.bed", emit: bed
        path "reference.fasta", emit: ref
        path "reference.gff3", emit: gff

    """
    cp ${scheme_name}/${scheme_version}/${scheme_name}.scheme.bed scheme.bed
    cp ${scheme_name}/${scheme_version}/${scheme_name}.reference.fasta reference.fasta
    cp ${scheme_name}/${scheme_version}/${scheme_name}.reference.gff3 reference.gff3
    """
}

process performHostFilter {
    cpus 1

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_hostfiltered_R*.fastq.gz", mode: 'copy'

    input:
        tuple val(sampleName), path(forward), path(reverse)
    output:
        tuple val(sampleName), path("${sampleName}_hostfiltered_R1.fastq.gz"), path("${sampleName}_hostfiltered_R2.fastq.gz"), emit: fastqPairs

    script:
        """
        bwa mem -t ${task.cpus} ${params.composite_ref} ${forward} ${reverse} | \
            filter_non_human_reads.py -c ${params.viral_contig_name} > ${sampleName}.viral_and_nonmapping_reads.bam
        samtools sort -n ${sampleName}.viral_and_nonmapping_reads.bam | \
             samtools fastq -1 ${sampleName}_hostfiltered_R1.fastq.gz -2 ${sampleName}_hostfiltered_R2.fastq.gz -s ${sampleName}_singletons.fastq.gz -
        """
}

process addCodonPositionToVariants {
    cpus 1

    executor 'local'

    tag { sampleId }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.with_codon_pos.tsv", mode: 'copy'

    input:
        tuple val(sampleName), path(variants), path(gff)
    output:
        tuple val(sampleName), path("${sampleName}.variants.with_codon_pos.tsv")

    script:
    def gff_arg = gff.name == 'NO_FILE' ? "" : "-g ${gff}"
        """
        ivar_variants_add_codon_position.py ${gff_arg} ${variants} > ${sampleName}.variants.with_codon_pos.tsv
        """
}