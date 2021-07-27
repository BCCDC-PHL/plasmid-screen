process get_reference_plasmid {

    tag { sample_id + " / " + plasmid_cluster_id + " / " + reference_plasmid_id }

    executor 'local'

    input:
      tuple val(sample_id), val(plasmid_cluster_id), val(reference_plasmid_id), path(mob_db)

    output:
      tuple val(sample_id), val(reference_plasmid_id), path("${reference_plasmid_id}.fa")

    script:
      
      """
      seqkit grep -r -p "${reference_plasmid_id}" ${mob_db}/ncbi_plasmid_full_seqs.fas > ${reference_plasmid_id}.fa 
      """
}
