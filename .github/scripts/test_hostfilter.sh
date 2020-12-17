#!/bin/bash
set -eo pipefail
export PATH=/opt/conda/bin:$PATH

git clone https://github.com/BCCDC-PHL/artic-ncov2019.git

export $REF_FILE=artic-ncov2019/primer_schemes/nCoV-2019/V1200/nCoV-2019.reference.fasta
export $bed_FILE=artic-ncov2019/primer_schemes/nCoV-2019/V1200/nCoV-2019.bed
echo ref file: $REF_FILE
echo bed file: $BED_FILE

# run current pull request code
# --sanger profile: there are only 2 available cpus in the github runner execution
sed -i s'/cpus = 4/cpus = 2/'g conf/coguk/sanger.config
singularity --version
echo Nextflow run --illumina mode with --ref, --bed --cram and outCram... >> artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run ./main.nf \
       -profile conda \
       --cache 
       --ref $REF_FILE \
       --bed $BED_FILE \
       --cram \
       --outCram \
       --directory $PWD/.github/data/crams/ \
       --illumina \
       --prefix test
cp .nextflow.log artifacts/Outcram_ref_bed.nextflow.log

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
