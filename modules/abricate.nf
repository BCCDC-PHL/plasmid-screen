process abricate {

    tag { assembly.baseName + " / " + database }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${assembly.baseName}_${database}_abricate.tsv", mode: 'copy'

    cpus 1

    input:
      tuple val(sample_id),  path(assembly), val(database)

    output:
      tuple val(sample_id), path("${assembly.baseName}_${database}_abricate.tsv")

    script:
      """
      abricate --db ${database} ${assembly} > ${assembly.baseName}_${database}_abricate.tsv
      """
}
