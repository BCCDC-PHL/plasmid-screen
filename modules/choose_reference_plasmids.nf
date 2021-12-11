process choose_reference_plasmids {

    tag { sample_id }

    executor 'local'

    input:
      tuple val(sample_id), path(combined_abricate_mobtyper_report)

    output:
      path("${sample_id}_reference_plasmids.csv")

    script:
      """
      choose_reference_plasmids.py -r ${combined_abricate_mobtyper_report} > ${sample_id}_reference_plasmids.csv
      """
}
