process join_mob_typer_and_abricate_reports {

    tag { mob_typer_report.baseName + " / " + abricate_report.baseName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${abricate_report.baseName}_resistance_plasmids.tsv", mode: 'copy'

    executor 'local'

    input:
      tuple val(sample_id),  path(mob_typer_report), path(abricate_report)

    output:
      tuple val(sample_id), path("${abricate_report.baseName}_resistance_plasmids.tsv")

    script:
      """
      join_mob_typer_and_abricate_reports.py --sample-id ${sample_id} --mob-typer-report ${mob_typer_report} --abricate-report ${abricate_report} > ${abricate_report.baseName}_resistance_plasmids.tsv
      """
}
