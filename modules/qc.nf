process makeQCCSV {

    tag { sampleName }

    label 'process_single'

    container 'community.wave.seqera.io/library/biopython_matplotlib_pandas:2b7509b8798f6d71'

    conda "conda-forge::matplotlib=3.9.2", "conda-forge::pandas=1.2.2", "conda-forge::biopython=1.84"

    publishDir "${params.outdir}/qc_plots", pattern: "${sampleName}.depth.png", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(bam_index), path(fasta), path(ref), path(primer_bed), path(primer_pairs)

    output:
    path "${params.prefix}.${sampleName}.qc.csv", emit: csv
    path "${sampleName}.depth.png"

    script:
    """
    qc.py --outfile ${params.prefix}.${sampleName}.qc.csv --sample ${sampleName} --ref ${ref} --bam ${bam} --fasta ${fasta} --primer-bed ${primer_bed} --min-depth ${params.varMinDepth}
    """
}


process writeQCSummaryCSV {
    tag { params.prefix }

    label 'process_single'

    container 'jitesoft/alpine:3.20.2'

    input:
    val lines

    exec:
    file("${params.outdir}/${params.prefix}.qc.csv").withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }   
    }
}

