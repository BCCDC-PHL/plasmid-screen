process mob_recon {

    tag { sample_id }

    cpus 4

    publishDir "${params.outdir}/${sample_id}/mob_recon", mode: 'copy', pattern: "${sample_id}/${sample_id}*"
 
    input:
      tuple val(sample_id), path(assembly)

    output:
      tuple val(sample_id), path("${sample_id}/${sample_id}*")

    script:
      """
      mob_recon -s ${sample_id} --num_threads ${task.cpus} --unicycler_contigs --run_overhang --run_typer --infile ${assembly} --outdir ${sample_id}
      sed -i 's/^>/>${sample_id}_/' ${sample_id}/*.fasta
      rename plasmid ${sample_id}_plasmid ${sample_id}/plasmid*.fasta
      rename chromosome ${sample_id}_chromosome ${sample_id}/chromosome.fasta
      rename mobtyper ${sample_id}_mobtyper ${sample_id}/mobtyper_results.txt
      """
}
