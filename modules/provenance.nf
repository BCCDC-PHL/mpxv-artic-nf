process collect_provenance {

  tag { sampleName }

  executor 'local'

  publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_*_provenance.yml", mode: 'copy'

  input:
  tuple val(sampleName), path(provenance_files)

  output:
  tuple val(sampleName), file("${sampleName}_*_provenance.yml")

  script:
  """
  cat ${provenance_files} > ${sampleName}_\$(date +%Y%m%d%H%M%S)_provenance.yml
  """
}


process pipeline_provenance {

  tag { pipeline_name + " / " + pipeline_version }

  executor 'local'

  input:
  tuple val(pipeline_name), val(pipeline_version), val(analysis_start)

  output:
  file("pipeline_provenance.yml")

  script:
  """
  printf -- "- pipeline_name: ${pipeline_name}\\n"             >> pipeline_provenance.yml
  printf -- "  pipeline_version: ${pipeline_version}\\n"       >> pipeline_provenance.yml
  printf -- "  timestamp_analysis_start: ${analysis_start}\\n" >> pipeline_provenance.yml
  """
}