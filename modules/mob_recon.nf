process mob_recon {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}*.{fasta,txt}"
 
    input:
      tuple val(sample_id), path(assembly), path(mob_plasmid_db)

    output:
      tuple val(sample_id), path("${sample_id}*.fasta"), emit: sequences
      tuple val(sample_id), path("${sample_id}_mobtyper*.txt"), emit: mobtyper_reports
      tuple val(sample_id), path("${sample_id}_contig*.txt"), emit: contig_reports
      
    script:
      """
      mob_recon \
        -s ${sample_id} \
        --num_threads ${task.cpus} \
        --database ${mob_plasmid_db} \
        -g ${params.mob_filter_db} \
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
      rename chromosome ${sample_id}_chromosome ${sample_id}/chromosome.fasta || true
      rename mobtyper ${sample_id}_mobtyper ${sample_id}/mobtyper_results.txt || true
      rename contig ${sample_id}_contig ${sample_id}/contig_report.txt
      mv ${sample_id}/* .
      """
}
