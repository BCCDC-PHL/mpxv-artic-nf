process performHostFilter {
    cpus 4

    tag { sampleId }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleId}_hostfiltered_R*.fastq.gz", mode: 'copy'

    input:
        tuple(val(sampleId), path(forward), path(reverse))
    output:
        tuple sampleId, path("${sampleId}_hostfiltered_R1.fastq.gz"), path("${sampleId}_hostfiltered_R2.fastq.gz"), emit: fastqPairs

    script:
        """
        bwa mem -t ${task.cpus} ${params.composite_ref} ${forward} ${reverse} | \
            filter_non_human_reads.py -c ${params.viral_contig_name} > ${sampleId}.viral_and_nonmapping_reads.bam
        samtools sort -n ${sampleId}.viral_and_nonmapping_reads.bam | \
             samtools fastq -1 ${sampleId}_hostfiltered_R1.fastq.gz -2 ${sampleId}_hostfiltered_R2.fastq.gz -s ${sampleId}_singletons.fastq.gz -
        """
}

process downsampleFastq {
    cpus 2

    tag { sampleId }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleId}_hostfiltered_downsampled_R*.fastq.gz", mode: 'copy'

    input:
        tuple val(sampleId), path(forward), path(reverse)
    output:
        tuple sampleId, path("${sampleId}_hostfiltered_downsampled_R1.fastq.gz"), path("${sampleId}_hostfiltered_downsampled_R2.fastq.gz"), emit: fastqPairs

    script:
        """
        zcat ${forward} | seqkit sample --rand-seed 11 --number ${params.downsample} -o ${sampleId}_hostfiltered_downsampled_R1.fastq.gz
	zcat ${reverse} | seqkit sample --rand-seed 11 --number ${params.downsample} -o ${sampleId}_hostfiltered_downsampled_R2.fastq.gz
        """
}
