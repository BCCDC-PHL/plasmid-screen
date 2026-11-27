process fastplong {

    tag { sample_id }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastplong:0.6.0--h43da1c4_0':
        'community.wave.seqera.io/library/fastplong:0.6.0--0b50463ff3c633a7' }"

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastplong.*", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_RL.filtered.fastq.gz", mode: 'copy', enabled: params.publish_filtered_reads && !params.dehost

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_fastplong.csv"),            emit: csv
    tuple val(sample_id), path("${sample_id}_fastplong.json"),           emit: report_json
    tuple val(sample_id), path("${sample_id}_fastplong.html"),           emit: report_html
    tuple val(sample_id), path("${sample_id}_RL.filtered.fastq.gz"),     emit: filtered_reads
    tuple val(sample_id), path("${sample_id}_fastplong_provenance.yml"), emit: provenance

    script:
    // Reads arrive with adapters, barcodes and primers already removed, so
    // adapter trimming stays off: filtering drops whole reads, nothing is trimmed.
    length_limit_arg = params.max_length > 0 ? "--length_limit ${params.max_length}" : ""
    """
    printf -- "- process_name: fastplong\\n"                                          >> ${sample_id}_fastplong_provenance.yml
    printf -- "  tools:\\n"                                                           >> ${sample_id}_fastplong_provenance.yml
    printf -- "    - tool_name: fastplong\\n"                                         >> ${sample_id}_fastplong_provenance.yml
    printf -- "      tool_version: \$(fastplong --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastplong_provenance.yml
    printf -- "      parameters:\\n"                                                  >> ${sample_id}_fastplong_provenance.yml
    printf -- "        - parameter: --disable_adapter_trimming\\n"                    >> ${sample_id}_fastplong_provenance.yml
    printf -- "          value: null\\n"                                              >> ${sample_id}_fastplong_provenance.yml
    printf -- "        - parameter: --length_required\\n"                             >> ${sample_id}_fastplong_provenance.yml
    printf -- "          value: ${params.min_length}\\n"                              >> ${sample_id}_fastplong_provenance.yml
    printf -- "        - parameter: --mean_qual\\n"                                   >> ${sample_id}_fastplong_provenance.yml
    printf -- "          value: ${params.min_mean_quality}\\n"                        >> ${sample_id}_fastplong_provenance.yml
    printf -- "        - parameter: --length_limit\\n"                                >> ${sample_id}_fastplong_provenance.yml
    printf -- "          value: ${params.max_length}\\n"                              >> ${sample_id}_fastplong_provenance.yml

    fastplong \
      --thread ${task.cpus} \
      --in ${reads} \
      --out ${sample_id}_RL.filtered.fastq.gz \
      --disable_adapter_trimming \
      --length_required ${params.min_length} \
      --mean_qual ${params.min_mean_quality} \
      ${length_limit_arg} \
      --report_title "fastplong report: ${sample_id}" \
      --json ${sample_id}_fastplong.json \
      --html ${sample_id}_fastplong.html

    parse_fastplong_report.py \
      --sample-id ${sample_id} \
      --json ${sample_id}_fastplong.json \
      --html ${sample_id}_fastplong.html \
      > ${sample_id}_fastplong.csv
    """
}


process fastplong_stats {

    tag { sample_id }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastplong:0.6.0--h43da1c4_0':
        'community.wave.seqera.io/library/fastplong:0.6.0--0b50463ff3c633a7' }"

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_post-dehosting_fastplong.*", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_post-dehosting_fastplong.csv"),            emit: csv
    tuple val(sample_id), path("${sample_id}_post-dehosting_fastplong.json"),           emit: report_json
    tuple val(sample_id), path("${sample_id}_post-dehosting_fastplong.html"),           emit: report_html
    tuple val(sample_id), path("${sample_id}_post-dehosting_fastplong_provenance.yml"), emit: provenance

    script:
    // Measures an already-filtered read set, so every filter is off and the report's
    // two sides are equal. --out is still given because fastplong only writes the
    // after-filtering sections when it produces output; the copy is discarded with
    // the work directory.
    """
    printf -- "- process_name: fastplong_post_dehosting\\n"                       >> ${sample_id}_post-dehosting_fastplong_provenance.yml
    printf -- "  tools:\\n"                                                       >> ${sample_id}_post-dehosting_fastplong_provenance.yml
    printf -- "    - tool_name: fastplong\\n"                                     >> ${sample_id}_post-dehosting_fastplong_provenance.yml
    printf -- "      tool_version: \$(fastplong --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_post-dehosting_fastplong_provenance.yml
    printf -- "      parameters:\\n"                                              >> ${sample_id}_post-dehosting_fastplong_provenance.yml
    printf -- "        - parameter: --disable_adapter_trimming\\n"                >> ${sample_id}_post-dehosting_fastplong_provenance.yml
    printf -- "          value: null\\n"                                          >> ${sample_id}_post-dehosting_fastplong_provenance.yml
    printf -- "        - parameter: --disable_quality_filtering\\n"               >> ${sample_id}_post-dehosting_fastplong_provenance.yml
    printf -- "          value: null\\n"                                          >> ${sample_id}_post-dehosting_fastplong_provenance.yml
    printf -- "        - parameter: --disable_length_filtering\\n"                >> ${sample_id}_post-dehosting_fastplong_provenance.yml
    printf -- "          value: null\\n"                                          >> ${sample_id}_post-dehosting_fastplong_provenance.yml

    fastplong \
      --thread ${task.cpus} \
      --in ${reads} \
      --out ${sample_id}_passthrough.fastq.gz \
      --disable_adapter_trimming \
      --disable_quality_filtering \
      --disable_length_filtering \
      --report_title "fastplong report: ${sample_id} (post-dehosting)" \
      --json ${sample_id}_post-dehosting_fastplong.json \
      --html ${sample_id}_post-dehosting_fastplong.html

    parse_fastplong_report.py \
      --sample-id ${sample_id} \
      --json ${sample_id}_post-dehosting_fastplong.json \
      --html ${sample_id}_post-dehosting_fastplong.html \
      > ${sample_id}_post-dehosting_fastplong.csv
    """
}