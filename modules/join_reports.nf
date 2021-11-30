process join_mob_typer_and_abricate_reports {

    tag { mob_typer_report.baseName + " / " + abricate_report.baseName }

    publishDir "${params.outdir}/${sample_id}", pattern: "${abricate_report.baseName}_resistance_plasmids.tsv", mode: 'copy'

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


process join_resistance_plasmid_and_snp_reports {

    tag { resistance_plasmid_report.baseName + " / " + snp_report.baseName }

    publishDir "${params.outdir}/${sample_id}", pattern: "${resistance_plasmid_report.baseName}_with_snps.tsv", mode: 'copy'

    executor 'local'

    input:
      tuple val(sample_id),  path(resistance_plasmid_report), val(plasmid_id), path(snp_report)

    output:
      tuple val(sample_id), path("${resistance_plasmid_report.baseName}_with_snps.tsv")

    script:
      """
      join_resistance_plasmid_and_snp_reports.py --sample-id ${sample_id} --resistance-plasmid-report ${resistance_plasmid_report} --snp-report ${snp_report} > ${resistance_plasmid_report.baseName}_with_snps.tsv
      """
}
