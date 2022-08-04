# mpxv-artic-nf
A Nextflow pipeline for running the ARTIC network's fieldbioinformatics tools (https://github.com/artic-network/fieldbioinformatics), with a focus on monkeypox virus (mpxv).

![push master](https://github.com/BCCDC-PHL/mpxv-artic-nf/actions/workflows/push_master.yml/badge.svg)

#### Introduction

This pipeline is based on the [BCCDC-PHL/ncov2019-artic-nf](https://github.com/BCCDC-PHL/ncov2019-artic-nf) pipeline, which is a fork of the [connor-lab/ncov2019-artic-nf](https://github.com/connor-lab/ncov2019-artic-nf) pipeline. It has been modified to support analysis of monkeypox virus.

```mermaid
flowchart TD
  ref[ref.fa]
  composite_ref[composite_ref.fa]
  primers[primer.bed]
  primer_pairs[primer_pairs.tsv]
  fastq[fastq_dir]
  fastq --> performHostFilter(performHostFilter)
  composite_ref --> performHostFilter
  performHostFilter(performHostFilter) --> normalizeDepth(normalizeDepth)
  readTrimming(readTrimming) --> filterResidualAdapters(filterResidualAdapters) 
  normalizeDepth(normalizeDepth) --> readTrimming(readTrimming)
  filterResidualAdapters --> readMapping(readMapping)
  ref --> readMapping(readMapping)
  readMapping(readMapping) --> trimPrimerSequences(trimPrimerSequences)
  primers --> trimPrimerSequences(trimPrimerSequences)
  trimPrimerSequences(trimPrimerSequences) --> callConsensusFreebayes(callConsensusFreebayes)
  callConsensusFreebayes(callConsensusFreebayes) --> alignConsensusToReference(alignConsensusToReference)
  ref --> alignConsensusToReference
  trimPrimerSequences --> makeQCCSV(makeQCCSV)
  callConsensusFreebayes --> makeQCCSV
  callConsensusFreebayes --> consensus[consensus.fa]
  callConsensusFreebayes --> variants[variants.vcf]
  ref --> makeQCCSV
  primers --> makeQCCSV
  primer_pairs --> makeQCCSV
  makeQCCSV --> qcCSV(qc.csv)
  makeQCCSV --> depthPNG(depth.png)
```

#### Quick-start

```
nextflow run BCCDC-PHL/mpxv-artic-nf -profile conda \
  --prefix "output_file_prefix" \
  --bed /path/to/primers.bed \
  --ref /path/to/ref.fa \
  --primer_pairs_tsv /path/to/primer_pairs_tsv \
  --composite_ref /path/to/human_and_mpxv_composite_ref \
  --directory /path/to/reads \
  --outdir /path/to/outputs
```


#### Installation
An up-to-date version of Nextflow is required because the pipeline is written in DSL2. Following the instructions at https://www.nextflow.io/ to download and install Nextflow should get you a recent-enough version. 


#### Conda
The repo contains a environment.yml files which automatically build the correct conda env if `-profile conda` is specifed in the command. Although you'll need `conda` installed, this is probably the easiest way to run this pipeline.

--cache /some/dir can be specified to have a fixed, shared location to store the conda build for use by multiple runs of the workflow.

#### Config

Important config options are:

| Option                           | Default  | Description                                                                                                         |
|:---------------------------------|---------:|--------------------------------------------------------------------------------------------------------------------:|
| `normalizationTargetDepth`       | `200`    | Target depth of coverage to normalize to prior to alignment                                                         |
| `normalizationMinDepth`          | `5`      | Minimum depth of coverage to normalize to prior to alignment                                                        |
| `keepLen`                        | `50`     | Length of reads to keep after primer trimming                                                                       |
| `qualThreshold`                  | `20`     | Sliding window quality threshold for keeping reads after primer trimming                                            |
| `varMinFreqThreshold`            | `0.25`   | Allele frequency threshold for ambiguous variant                                                                    |
| `varFreqThreshold`               | `0.75`   | Allele frequency threshold for unambiguous variant                                                                  |
| `varMinDepth`                    | `10`     | Minimum coverage depth to call variant                                                                              |


#### QC
A script to do some basic QC is provided in `bin/qc.py`. It measures the % of reference bases are covered by `varMinDepth`, and the longest stretch of consensus sequence with no `N` bases. This script does not make a QC pass/fail call.

#### Output
A subdirectory for each process in the workflow is created in `--outdir`. A `nml_upload` subdirectory containing dehosted fastq files and consensus sequences is included. 
