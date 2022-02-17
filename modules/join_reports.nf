process join_mob_typer_and_abricate_reports {

    tag { sample_id }

    executor 'local'

    input:
      tuple val(sample_id),  path(mob_typer_report), path(abricate_report)

    output:
      tuple val(sample_id), path("${sample_id}_resistance_plasmids.tsv")

    script:
      """
      join_mob_typer_and_abricate_reports.py --sample-id ${sample_id} --mob-typer-report ${mob_typer_report} --abricate-report ${abricate_report} > ${sample_id}_resistance_plasmids.tsv
      """
}

process select_resistance_chromosomes {

    tag { sample_id }

    executor 'local'

    input:
      tuple val(sample_id),  path(resistance_gene_report)

    output:
      tuple val(sample_id), path("${sample_id}_resistance_chromosomes.tsv")

    script:
      """
      head -n 1 ${resistance_gene_report} > ${sample_id}_resistance_chromosomes.tsv
      awk -F \$'\\t' '\$2 ~ "chromosome"' ${resistance_gene_report} >> ${sample_id}_resistance_chromosomes.tsv
      """
}

process join_resistance_plasmid_and_snp_reports {

    tag { sample_id + " / " + plasmid_id + " / " + resistance_gene }

    executor 'local'

    input:
      tuple val(sample_id), path(resistance_plasmid_report), val(plasmid_id), val(resistance_gene), path(snp_report), path(coverage_report)

    output:
      tuple val(sample_id), path("${sample_id}_${plasmid_id}_${resistance_gene}_resistance_plasmids.tsv")

    script:
      """
      join_resistance_plasmid_and_snp_reports.py --sample-id ${sample_id} --resistance-plasmid-report ${resistance_plasmid_report} --snp-report ${snp_report} --coverage-report ${coverage_report} > ${sample_id}_${plasmid_id}_${resistance_gene}_resistance_plasmids.tsv
      """
}

process concatenate_resistance_reports {

    tag { sample_id }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_resistance_gene_report.tsv", mode: 'copy'

    executor 'local'

    input:
      tuple val(sample_id),  path(resistance_reports)

    output:
      tuple val(sample_id), path("${sample_id}_resistance_gene_report.tsv")

    script:
      """
      head -n 1 ${resistance_reports[0]} > ${sample_id}_resistance_gene_report.tsv
      tail -qn+2 ${resistance_reports} >> ${sample_id}_resistance_gene_report.tsv
      """
}