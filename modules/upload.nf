process collateSamples {

    tag { sampleName }

    publishDir "${params.outdir}/nml_upload/", pattern: "${params.prefix}/${sampleName}*", mode: 'copy'

    input:
    tuple val(sampleName), path(fasta), path(fastq_r1), path(fastq_r2)

    output:
    path("${params.prefix}/${sampleName}*")

    script:
    """
    mkdir ${params.prefix} && cp ${sampleName}* ${params.prefix}
    """
}
