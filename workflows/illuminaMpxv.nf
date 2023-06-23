#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { articDownloadScheme }       from '../modules/utils.nf'
include { normalizeDepth }            from '../modules/illumina.nf'
include { performHostFilter }         from '../modules/illumina.nf'
include { readTrimming }              from '../modules/illumina.nf'
include { filterResidualAdapters }    from '../modules/illumina.nf'
include { indexReference }            from '../modules/illumina.nf'
include { readMapping }               from '../modules/illumina.nf'
include { trimPrimerSequences }       from '../modules/illumina.nf'
include { callConsensusFreebayes }    from '../modules/illumina.nf'
include { annotateVariantsVCF }       from '../modules/illumina.nf'
include { alignConsensusToReference } from '../modules/illumina.nf'

include { makeQCCSV }         from '../modules/qc.nf'
include { writeQCSummaryCSV } from '../modules/qc.nf'

include { collateSamples }    from '../modules/upload.nf'

workflow prepareReferenceFiles {
    // Get reference fasta
    if (params.ref) {
      Channel.fromPath(params.ref)
              .set{ ch_refFasta }
    } else {
      articDownloadScheme()
      articDownloadScheme.out.reffasta
                          .set{ ch_refFasta }
    }


    /* Either get BWA aux files from reference 
       location or make them fresh */
    
    if (params.ref) {
      // Check if all BWA aux files exist, if not, make them
      bwaAuxFiles = []
      refPath = new File(params.ref).getAbsolutePath()
      new File(refPath).getParentFile().eachFileMatch( ~/.*.bwt|.*.pac|.*.ann|.*.amb|.*.sa/) { bwaAuxFiles << it }
     
      if ( bwaAuxFiles.size() == 5 ) {
        Channel.fromPath( bwaAuxFiles )
               .set{ ch_bwaAuxFiles }

        ch_refFasta.combine(ch_bwaAuxFiles.collect().toList())
                   .set{ ch_preparedRef }
      } else {
        indexReference(ch_refFasta)
        indexReference.out
                      .set{ ch_preparedRef }
      }
    } else {
      indexReference(ch_refFasta)
      indexReference.out
                    .set{ ch_preparedRef }
    }
  
    /* If bedfile is supplied, use that,
       if not, get it from ARTIC github repo */ 
 
    if (params.bed) {
      ch_bedFile = Channel.fromPath(params.bed)
    } else {
      ch_bedFile = articDownloadScheme.out.bed
    }

    ch_gff = Channel.fromPath(params.gff)

    ch_primerPairs = Channel.fromPath(params.primer_pairs_tsv)

    emit:
      bwaindex = ch_preparedRef
      bedfile = ch_bedFile
      primer_pairs = ch_primerPairs
      gff = ch_gff
}


workflow sequenceAnalysis {
    take:
      ch_filePairs
      ch_preparedRef
      ch_bedFile
      ch_primerPairs
      ch_gff

    main:

      if (!params.skip_normalize_depth) {
        ch_reads_to_hostfilter = normalizeDepth(ch_filePairs)
      } else {
        ch_reads_to_hostfilter = ch_filePairs
      }

      performHostFilter(ch_reads_to_hostfilter)

      readTrimming(performHostFilter.out.fastqPairs)

      filterResidualAdapters(readTrimming.out)

      readMapping(filterResidualAdapters.out.combine(ch_preparedRef))

      trimPrimerSequences(readMapping.out.combine(ch_bedFile))

      callConsensusFreebayes(trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] }))

      if (params.gff != 'NO_FILE') {
        annotateVariantsVCF(callConsensusFreebayes.out.variants.combine(ch_preparedRef.map{ it[0] }).combine(ch_gff))
      }

      alignConsensusToReference(callConsensusFreebayes.out.consensus.combine(ch_preparedRef.map{ it[0] }))

      makeQCCSV(trimPrimerSequences.out.ptrim.join(callConsensusFreebayes.out.consensus, by: 0)
                                   .combine(ch_preparedRef.map{ it[0] })
				   .combine(ch_bedFile)
				   .combine(ch_primerPairs))

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .set { qc }

      writeQCSummaryCSV(qc.toList())

      collateSamples(callConsensusFreebayes.out.consensus.join(performHostFilter.out.fastqPairs))

    emit:
      qc_pass = collateSamples.out
}

workflow mpxvIllumina {
    take:
      ch_filePairs

    main:
      prepareReferenceFiles()
      sequenceAnalysis(ch_filePairs, prepareReferenceFiles.out.bwaindex, prepareReferenceFiles.out.bedfile, prepareReferenceFiles.out.primer_pairs, prepareReferenceFiles.out.gff)
}


