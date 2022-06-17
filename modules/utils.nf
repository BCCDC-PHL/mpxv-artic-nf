process performHostFilter {
    cpus 1

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
    cpus 1

    tag { sampleId }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleId}.mapped.primertrimmed.downsampled.sorted{.bam,.bam.bai}", mode: 'copy'

    input:
        tuple val(sampleId), path(trimmed_bam), path(trimmed_bam_index), path(bedfile), path(ref)
    output:
        tuple val(sampleId), path("${sampleId}.mapped.primertrimmed.downsampled.sorted.bam"), path("${sampleId}.mapped.primertrimmed.downsampled.sorted.bam.bai"), emit: alignment
	path("downsampling_summary.csv"), emit: summary

    script:
        """
        samtools faidx ${ref}
        downsample_amplicons.py \
          --bed ${bedfile} \
          --min-depth ${params.downsampleMinDepth} \
          --mapping-quality ${params.downsampleMappingQuality} \
          --amplicon-subdivisions ${params.downsampleAmpliconSubdivisions} \
          --genome-size \$(cut -d \$'\\t' -f 2 ${ref}.fai) \
          ${trimmed_bam} 2> downsampling_summary_no_sample_id.csv | \
            samtools sort - -o ${sampleId}.mapped.primertrimmed.downsampled.sorted.bam
        samtools index ${sampleId}.mapped.primertrimmed.downsampled.sorted.bam
        paste -d ',' <(echo "sample_id" && echo "${sampleId}") downsampling_summary_no_sample_id.csv > downsampling_summary.csv
        """
}

process downsampledBamToFastq {
    cpus 1

    tag { sampleId }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleId}_hostfiltered_downsampled_R{1,2}.fastq.gz", mode: 'copy'

    input:
        tuple val(sampleId), path(downsampled_bam), path(downsampled_bam_index) 
    output:
        tuple val(sampleId), path("${sampleId}_hostfiltered_downsampled_R1.fastq.gz"), path("${sampleId}_hostfiltered_downsampled_R2.fastq.gz"), emit: fastqPairs

    script:
        """
        samtools collate -u ${downsampled_bam} -o ${sampleId}.collate.bam ${sampleId}
        samtools fastq -n -0 /dev/null -1 ${sampleId}_hostfiltered_downsampled_R1.fastq.gz -2 ${sampleId}_hostfiltered_downsampled_R2.fastq.gz -s ${sampleId}_hostfiltered_downsampled_S.fastq.gz ${sampleId}.collate.bam
        """
}

process addCodonPositionToVariants {
    cpus 1

    executor 'local'

    tag { sampleId }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleId}.variants.with_codon_pos.tsv", mode: 'copy'

    input:
        tuple val(sampleId), path(variants), path(gff)
    output:
        tuple val(sampleId), path("${sampleId}.variants.with_codon_pos.tsv")

    script:
    def gff_arg = gff.name == 'NO_FILE' ? "" : "-g ${gff}"
        """
        ivar_variants_add_codon_position.py ${gff_arg} ${variants} > ${sampleId}.variants.with_codon_pos.tsv
        """
}