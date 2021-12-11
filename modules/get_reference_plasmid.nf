process get_reference_plasmid {

    tag { sample_id + " / " + reference_plasmid_id + " / " + resistance_gene }

    executor 'local'

    input:
      tuple val(sample_id), val(reference_plasmid_id), val(resistance_gene), path(mob_db)

    output:
      tuple val(sample_id), val(reference_plasmid_id), val(resistance_gene), path("${reference_plasmid_id}.fa")

    script:
      
      """
      seqkit grep -r -p "${reference_plasmid_id}" ${mob_db}/ncbi_plasmid_full_seqs.fas > ${reference_plasmid_id}.fa 
      """
}
