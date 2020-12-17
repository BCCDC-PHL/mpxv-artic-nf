#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

include {performHostFilter} from './modules/utils.nf'


workflow {
  Channel.fromFilePairs( "${params.fastqSearchPath}/*_R{1,2}*.fastq*", type: 'file', maxDepth: 1).map{ it -> [it[0] - ~/_\w+/, it[1][0], it[1][1]] }.set{ ch_filePairs }

  
  
  main:
    ch_filePairs.view()
    performHostFilter(ch_filePairs)
}
