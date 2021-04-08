#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme } from '../modules/artic.nf' 
include {readTrimming} from '../modules/illumina.nf' 
include {filterResidualAdapters} from '../modules/illumina.nf' 
include {indexReference} from '../modules/illumina.nf'
include {readMapping} from '../modules/illumina.nf' 
include {trimPrimerSequences} from '../modules/illumina.nf'
include {callConsensusFreebayes} from '../modules/illumina.nf'
include {annotateVariantsVCF} from '../modules/illumina.nf'
include {callVariants} from '../modules/illumina.nf'
include {makeConsensus} from '../modules/illumina.nf' 
include {cramToFastq} from '../modules/illumina.nf'
include {alignConsensusToReference} from '../modules/illumina.nf'
include {trimUTRFromAlignment} from '../modules/illumina.nf'
include {performHostFilter} from '../modules/utils.nf'
include {downsampleAmplicons} from '../modules/utils.nf'
include {downsampledBamToFastq} from '../modules/utils.nf'
include {addCodonPositionToVariants} from '../modules/utils.nf'

include {makeQCCSV} from '../modules/qc.nf'
include {writeQCSummaryCSV} from '../modules/qc.nf'

include {bamToCram} from '../modules/out.nf'

include {collateSamples} from '../modules/upload.nf'

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
 
    if (params.bed ) {
      Channel.fromPath(params.bed)
             .set{ ch_bedFile }

    } else {
      articDownloadScheme.out.bed
                         .set{ ch_bedFile }
    }

    Channel.fromPath(params.gff)
           .set{ ch_gff }

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

      performHostFilter(ch_filePairs)

      readTrimming(performHostFilter.out)

      filterResidualAdapters(readTrimming.out)

      readMapping(filterResidualAdapters.out.combine(ch_preparedRef))

      trimPrimerSequences(readMapping.out.combine(ch_bedFile))

      downsampleAmplicons(trimPrimerSequences.out.ptrim.combine(ch_bedFile))

      downsampledBamToFastq(downsampleAmplicons.out.alignment)

      callVariants(downsampleAmplicons.out.alignment.combine(ch_preparedRef.map{ it[0] }).combine(ch_gff))

      addCodonPositionToVariants(callVariants.out.combine(ch_gff))

      callConsensusFreebayes(downsampleAmplicons.out.alignment.combine(ch_preparedRef.map{ it[0] }))

      if (params.gff != 'NO_FILE') {
        annotateVariantsVCF(callConsensusFreebayes.out.variants.combine(ch_preparedRef.map{ it[0] }).combine(ch_gff))
      }
      
      makeConsensus(downsampleAmplicons.out.alignment)

      alignConsensusToReference(makeConsensus.out.combine(ch_preparedRef.map{ it[0] }))

      trimUTRFromAlignment(alignConsensusToReference.out)

      makeQCCSV(downsampleAmplicons.out.alignment.join(callConsensusFreebayes.out.consensus, by: 0)
                                   .combine(ch_preparedRef.map{ it[0] }))

      downsampleAmplicons.out.summary.collectFile(keepHeader: true, sort: { it.text }, name: "${params.prefix}.downsampling.csv", storeDir: "${params.outdir}")

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
    		       }
                       .set { qc }

      writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

      collateSamples(callConsensusFreebayes.out.consensus.join(downsampledBamToFastq.out.fastqPairs))

      if (params.outCram) {
        bamToCram(trimPrimerSequences.out.mapped.map{ it[0] } 
                        .join (trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] })) )

      }

    emit:
      qc_pass = collateSamples.out
}

workflow ncovIllumina {
    take:
      ch_filePairs

    main:
      // Build or download fasta, index and bedfile as required
      prepareReferenceFiles()

      // Actually do analysis
      sequenceAnalysis(ch_filePairs, prepareReferenceFiles.out.bwaindex, prepareReferenceFiles.out.bedfile, prepareReferenceFiles.out.gff)
}

workflow ncovIlluminaCram {
    take:
      ch_cramFiles
    main:
      // Convert CRAM to fastq
      cramToFastq(ch_cramFiles)

      // Run standard pipeline
      ncovIllumina(cramToFastq.out)
}

