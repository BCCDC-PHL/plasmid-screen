process abricate {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_abricate.tsv", mode: 'copy'

    cpus 1

    input:
      tuple val(sample_id),  path(assembly)

    output:
      tuple val(sample_id), path("${sample_id}_abricate.tsv")

    script:
      """
      abricate --db ncbi ${assembly} > ${sample_id}_abricate.tsv
      """
}
