process abricate {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_abricate.tsv", mode: 'copy'

    cpus 1

    input:
      tuple val(sample_id),  path(assemblies)

    output:
      tuple val(sample_id), path("${sample_id}_abricate.tsv"), emit: report
      tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
      """
      printf -- "- process_name: abricate\\n" > ${sample_id}_abricate_provenance.yml
      printf -- "  tool_name: abricate\\n  tool_version: \$(abricate --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_abricate_provenance.yml
      printf -- "  parameters:\\n" >> ${sample_id}_abricate_provenance.yml
      printf -- "  - parameter: db\\n    value: ncbi\\n" >> ${sample_id}_abricate_provenance.yml
      abricate --db ncbi ${assemblies} > ${sample_id}_abricate.tsv
      """
}
