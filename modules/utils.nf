process performHostFilter {
    cpus 4

    tag { sampleId }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleId}_hostfiltered_R*.fastq.gz", mode: 'copy'

    input:
        tuple val(sampleId), path(forward), path(reverse)
    output:
        tuple val(sampleId), path("${sampleId}_hostfiltered_R1.fastq.gz"), path("${sampleId}_hostfiltered_R2.fastq.gz"), emit: fastqPairs

    script:
        """
        bwa mem -t ${task.cpus} ${params.composite_ref} ${forward} ${reverse} | \
            filter_non_human_reads.py -c ${params.viral_contig_name} > ${sampleId}.viral_and_nonmapping_reads.bam
        samtools sort -n ${sampleId}.viral_and_nonmapping_reads.bam | \
             samtools fastq -1 ${sampleId}_hostfiltered_R1.fastq.gz -2 ${sampleId}_hostfiltered_R2.fastq.gz -s ${sampleId}_singletons.fastq.gz -
        """
}

process downsampleAmplicons {
    cpus 2

    tag { sampleId }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleId}.mapped.primertrimmed.downsampled.sorted{.bam,.bam.bai}", mode: 'copy'

    input:
        tuple val(sampleId), path(trimmed_bam), path(trimmed_bam_index), path(bedfile)
    output:
        tuple val(sampleId), path("${sampleId}.mapped.primertrimmed.downsampled.sorted.bam"), path("${sampleId}.mapped.primertrimmed.downsampled.sorted.bam.bai")

    script:
        """
        downsample_amplicons.py --bed ${bedfile} --depth ${params.downsampleDepth} --mapping-quality ${params.downsampleMappingQuality} --amplicon-subdivisions ${params.downsampleAmpliconSubdivisions}  ${trimmed_bam} | \
            samtools sort - -o ${sampleId}.mapped.primertrimmed.downsampled.sorted.bam
        samtools index ${sampleId}.mapped.primertrimmed.downsampled.sorted.bam
        """
}
