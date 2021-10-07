process choose_reference_plasmid_for_sample {

    executor 'local'

    input:
      tuple val(sample_id), path(combined_abricate_mobtyper_report)

    output:
      path("plasmid_cluster_for_sample.tsv")

    script:
      """
      get_plasmid_cluster_id_for_sample.py -r ${combined_abricate_mobtyper_report} > plasmid_cluster_for_sample.tsv
      """
}
