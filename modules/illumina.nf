process performHostFilter {

    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_hostfiltered_R*.fastq.gz", mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("${sampleName}_hostfiltered_R1.fastq.gz"), path("${sampleName}_hostfiltered_R2.fastq.gz"), emit: fastqPairs

    script:
    """
    bwa mem -t ${task.cpus} ${params.composite_ref} ${forward} ${reverse} | \
      filter_non_human_reads.py -c ${params.viral_contig_name} > ${sampleName}.viral_and_nonmapping_reads.bam
    samtools sort -@ ${task.cpus} -n ${sampleName}.viral_and_nonmapping_reads.bam | \
      samtools fastq -1 ${sampleName}_hostfiltered_R1.fastq.gz -2 ${sampleName}_hostfiltered_R2.fastq.gz -s ${sampleName}_singletons.fastq.gz -
    """
}

process normalizeDepth {

    tag { sampleName }

    memory '16 GB'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_norm_R{1,2}.fq.gz', mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("${sampleName}_norm_R1.fq.gz"), path("${sampleName}_norm_R2.fq.gz")

    script:
    """
    bbnorm.sh \
      target=${params.normalizationTargetDepth} \
      mindepth=${params.normalizationMinDepth} \
      in=${forward} \
      in2=${reverse} \
      out=${sampleName}_norm_R1.fq.gz \
      out2=${sampleName}_norm_R2.fq.gz
    """
}

process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    cpus 1

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), optional: true

    script:
    """
    if [[ \$(gunzip -c ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      cp ${forward} ${sampleName}_hostfiltered_val_1.fq.gz
      cp ${reverse} ${sampleName}_hostfiltered_val_2.fq.gz
    else
      trim_galore --paired $forward $reverse
    fi
    """
}

process filterResidualAdapters {
    /**
    * Discard reads that contain residual adapter sequences that indicate trimming may have failed
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output untrim_filter_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*{1,2}_posttrim_filter.fq.gz', mode: 'copy'

    cpus 1

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("*1_posttrim_filter.fq.gz"), path("*2_posttrim_filter.fq.gz")

    script:
    """
    filter_residual_adapters.py --input_R1 $forward --input_R2 $reverse
    """
}

process indexReference {
    /**
    * Indexes reference fasta file in the scheme repo using bwa.
    */

    tag { ref }

    input:
    path(ref)

    output:
    tuple path('ref.fa'), path('ref.fa.*')

    script:
    """
    ln -s ${ref} ref.fa
    bwa index ref.fa
    """
}

process readMapping {
    /**
    * Maps trimmed paired fastq using BWA (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input 
    * @output 
    */

    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted{.bam,.bam.bai}", mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse), path(ref), path("*")

    output:
    tuple val(sampleName), path("${sampleName}.sorted.bam"), path("${sampleName}.sorted.bam.bai")

    script:
    """
    bwa mem -t ${task.cpus} ${ref} ${forward} ${reverse} | \
    samtools sort -o ${sampleName}.sorted.bam
    samtools index ${sampleName}.sorted.bam
    """
}

process trimPrimerSequences {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped{.bam,.bam.bai}", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted{.bam,.bam.bai}", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(bam_index), path(bedfile)

    output:
    tuple val(sampleName), path("${sampleName}.mapped.bam"), path("${sampleName}.mapped.bam.bai"), emit: mapped
    tuple val(sampleName), path("${sampleName}.mapped.primertrimmed.sorted.bam"), path("${sampleName}.mapped.primertrimmed.sorted.bam.bai" ), emit: ptrim

    script:
    """
    samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
    samtools index ${sampleName}.mapped.bam
    ivar trim -e -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.keepLen} -q ${params.qualThreshold} -f ${params.primer_pairs_tsv} -p ivar.out
    samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
    samtools index ${sampleName}.mapped.primertrimmed.sorted.bam
    """
}

process callConsensusFreebayes {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.consensus.fa", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.norm.vcf", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(bam_index), path(ref)

    output:
    tuple val(sampleName), path("${sampleName}.consensus.fa"), emit: consensus
    tuple val(sampleName), path("${sampleName}.variants.norm.vcf"), emit: variants

    script:
    """
    # the sed is to fix the header until a release is made with this fix
    # https://github.com/freebayes/freebayes/pull/549
    freebayes -p 1 \
              -f ${ref} \
              -F 0.2 \
              -C 1 \
              --pooled-continuous \
              --min-coverage ${params.varMinDepth} \
              --gvcf --gvcf-dont-use-chunk true ${bam} | sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > ${sampleName}.gvcf

    # make depth mask, split variants into ambiguous/consensus
    # NB: this has to happen before bcftools norm or else the depth mask misses any bases exposed during normalization
    process_gvcf.py -d ${params.varMinDepth} \
                    -l ${params.varMinFreqThreshold} \
                    -u ${params.varFreqThreshold} \
                    -m ${sampleName}.mask.txt \
                    -v ${sampleName}.variants.vcf \
                    -c ${sampleName}.consensus.vcf ${sampleName}.gvcf

    # normalize variant records into canonical VCF representation
    for v in "variants" "consensus"; do
        bcftools norm -f ${ref} ${sampleName}.\$v.vcf > ${sampleName}.\$v.norm.vcf
    done

    # split the consensus sites file into a set that should be IUPAC codes and all other bases, using the ConsensusTag in the VCF
    for vt in "ambiguous" "fixed"; do
        cat ${sampleName}.consensus.norm.vcf | awk -v vartag=ConsensusTag=\$vt '\$0 ~ /^#/ || \$0 ~ vartag' > ${sampleName}.\$vt.norm.vcf
        bgzip -f ${sampleName}.\$vt.norm.vcf
        tabix -f -p vcf ${sampleName}.\$vt.norm.vcf.gz
    done

    # apply ambiguous variants first using IUPAC codes. this variant set cannot contain indels or the subsequent step will break
    bcftools consensus -f ${ref} -I ${sampleName}.ambiguous.norm.vcf.gz > ${sampleName}.ambiguous.fa

    # apply remaninng variants, including indels
    bcftools consensus -f ${sampleName}.ambiguous.fa -m ${sampleName}.mask.txt ${sampleName}.fixed.norm.vcf.gz | sed s/${params.viral_contig_name}/${sampleName}/ > ${sampleName}.consensus.fa
    """
}

process alignConsensusToReference {
    /**
    * Aligns consensus sequence against reference using mafft. Uses the --keeplength
    * flag to guarantee that all alignments remain the same length as the reference.
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.consensus.aln.fa", mode: 'copy'

    input:
    tuple val(sampleName), path(consensus), path(reference)

    output:
    tuple val(sampleName), path("${sampleName}.consensus.aln.fa")

    script:
    // Convert multi-line fasta to single line
    awk_string = '/^>/ {printf("\\n%s\\n", $0); next; } { printf("%s", $0); }  END { printf("\\n"); }'
    """
    mafft \
      --preservecase \
      --keeplength \
      --add \
      ${consensus} \
      ${reference} \
      > ${sampleName}.with_ref.multi_line.alignment.fa
    awk '${awk_string}' ${sampleName}.with_ref.multi_line.alignment.fa > ${sampleName}.with_ref.single_line.alignment.fa
    tail -n 2 ${sampleName}.with_ref.single_line.alignment.fa > ${sampleName}.consensus.aln.fa
    """
}

process annotateVariantsVCF {
    /**
    */

    tag { sampleName }

    cpus 1

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.norm.consequence.{vcf,tsv}", mode: 'copy'

    input:
        tuple val(sampleName), path(vcf), path(ref), path(gff)

    output:
        tuple val(sampleName), path("${sampleName}.variants.norm.consequence.vcf"), emit: vcf
        tuple val(sampleName), path("${sampleName}.variants.norm.consequence.tsv"), emit: tsv

    script:
        """
        bcftools csq -f ${ref} -g ${gff} ${vcf} -Ov -o ${sampleName}.variants.norm.consequence.vcf
        bcftools csq -f ${ref} -g ${gff} ${vcf} -Ot -o ${sampleName}.variants.norm.consequence.tsv
        """
}
