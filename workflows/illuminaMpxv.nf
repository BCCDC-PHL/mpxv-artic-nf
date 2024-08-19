#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { get_bed_ref }       from '../modules/utils.nf'
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
    // if (params.ref) {
    //   Channel.fromPath(params.ref)
    //           .set{ ch_refFasta }
    // } else {
    //   articDownloadScheme()
    //   articDownloadScheme.out.reffasta
    //                       .set{ ch_refFasta }
    // }

    scheme_dir_name = "primer-schemes"
    schemes = """./data/${scheme_dir_name}/${params.scheme_name}"""
    scheme_dir = file(projectDir.resolve(schemes), type:'file', checkIfExists:true)

    get_bed_ref(params.schemeDir, params.scheme, params.schemeVersion)
    ch_bedFile = get_bed_ref.out.bed
    ch_refFasta = get_bed_ref.out.ref
    ch_gff = get_bed_ref.out.gff

    /* Either get BWA aux files from reference 
      location or make them fresh */
  
    // Check if all BWA aux files exist, if not, make them
    indexReference(ch_refFasta)
    indexReference.out
                  .set{ ch_preparedRef }

    emit:
      bwaindex = ch_preparedRef
      bedfile = ch_bedFile
      gff = ch_gff
}


workflow sequenceAnalysis {
    take:
      ch_filePairs
      ch_preparedRef
      ch_bedFile
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
      sequenceAnalysis(ch_filePairs, prepareReferenceFiles.out.bwaindex, prepareReferenceFiles.out.bedfile, prepareReferenceFiles.out.gff)
}


