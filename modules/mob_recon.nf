process mob_recon {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}*.{fasta,tsv}"
 
    input:
      tuple val(sample_id), path(assembly)

    output:
      tuple val(sample_id), path("${sample_id}*.fasta"), emit: sequences
      tuple val(sample_id), path("${sample_id}_mobtyper_plasmid_report.tsv"), emit: mobtyper_reports, optional: true
      tuple val(sample_id), path("${sample_id}_mobtyper_contig_report.tsv"), emit: contig_reports
      tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance
      
    script:
      """
      
      printf -- "- process_name: mob_recon\\n" > ${sample_id}_mob_recon_provenance.yml
      printf -- "  tool_name: mob_recon\\n  tool_version: \$(mob_recon --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_mob_recon_provenance.yml
      printf -- "  - parameter: database_directory\\n    value: \$(realpath ${params.mob_db})\\n" >> ${sample_id}_mob_recon_provenance.yml
      printf -- "  - parameter: filter_db\\n    value: \$(realpath ${params.mob_filter_db})\\n" >> ${sample_id}_mob_recon_provenance.yml
      printf -- "  - parameter: min_con_cov\\n    value: 95\\n" >> ${sample_id}_mob_recon_provenance.yml

      mob_recon \
        -s ${sample_id} \
        --num_threads ${task.cpus} \
        --unicycler_contigs \
        --database_directory ${params.mob_db} \
        --filter_db ${params.mob_filter_db} \
        --min_con_cov 95 \
        --run_overhang \
        --run_typer \
        --infile ${assembly} \
        --outdir ${sample_id} \
        --force

      for assembly in ${sample_id}/*.fasta;
      do
        sed -i -r 's/^>${sample_id}_([[:digit:]]+)_(.+)/>${sample_id}_\\1 \\2/' \${assembly};
      done
      rename plasmid ${sample_id}_plasmid ${sample_id}/plasmid*.fasta || true
      cp ${sample_id}/${sample_id}_plasmid*.fasta . || true
      cp ${sample_id}/chromosome.fasta ${sample_id}_chromosome.fasta || true
      cp ${sample_id}/mobtyper_results.txt ${sample_id}_mobtyper_plasmid_report.tsv || true
      cp ${sample_id}/contig_report.txt ${sample_id}_mobtyper_contig_report.tsv
      """
}
