process hash_files {

    tag { sampleName + " / " + file_type }

    input:
    tuple  val(sampleName), path(files_to_hash), val(file_type)

    output:
    tuple  val(sampleName), path("${sampleName}_${file_type}.sha256.csv"), emit: csv
    tuple  val(sampleName), path("${sampleName}_${file_type}_provenance.yml"), emit: provenance

    script:
    """
    shasum -a 256 ${files_to_hash} | tr -s ' ' ',' > ${sampleName}_${file_type}.sha256.csv
    while IFS=',' read -r hash filename; do
      printf -- "- input_filename: \$filename\\n"  >> ${sampleName}_${file_type}_provenance.yml;
      printf -- "  file_type: ${file_type}\\n"     >> ${sampleName}_${file_type}_provenance.yml;
      printf -- "  sha256: \$hash\\n"              >> ${sampleName}_${file_type}_provenance.yml;
    done < ${sampleName}_${file_type}.sha256.csv
    """

}