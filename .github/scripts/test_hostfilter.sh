#!/bin/bash
set -eo pipefail
export PATH=/opt/conda/bin:$PATH

git clone https://github.com/BCCDC-PHL/artic-ncov2019.git

export $REF_FILE=$GITHUB_WORKSPACE/.github/data/fastas/composite_nCoV-2019_chr1_regions.fa
echo ref file: $REF_FILE

echo Nextflow run tests... >> artifacts/test_artifact.log
NXF_VER=20.10.0 nextflow run ./test.nf \
       -profile conda \
       --cache $REPO/conda_cache_dir \
       --composite_ref $REF_FILE \
       --fastqSearchPath $GITHUB_WORKSPACE/.github/data/fastqs
       --illumina \
       --viral_contig_name 'MN908947.3:1-500' \
       --outdir test_hostfilter_out
cp .nextflow.log artifacts/nextflow_tests.nextflow.log

# get first cram file that macthes in given folder 
cram_exists=($(find results/*_bamToCram -name '*.cram'))

if [ -f "$cram_exists" ]; then
      echo Cram file generated successfully >> artifacts/test_artifact.log
      # clean tests
      rm -rf results && rm -rf work && rm -rf .nextflow*
      exit 0
else
      echo Cram file generation failed >> artifacts/test_artifact.log
      # clean-up for following tests
      echo cleaning tests
      rm -rf results && rm -rf work && rm -rf .nextflow*
      exit 1
fi
