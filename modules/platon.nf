process platon {

    tag { sample_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}*.{chromosome.fasta,plasmid.fasta,tsv,json}", mode: 'copy'

    cpus 4

    input:
      tuple val(sample_id),  path(assembly), path(platon_db)

    output:
      tuple val(sample_id), path("${sample_id}*.chromosome.fasta"), path("${sample_id}*.plasmid.fasta"), path("${sample_id}*.tsv"), path("${sample_id}*.json")

    script:
      """
      platon --threads ${task.cpus} --db ${platon_db} ${assembly}
      """
}
