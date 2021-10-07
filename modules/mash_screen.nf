process mash_screen {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.mash_screen.tsv", mode: 'copy'

    input:
      tuple val(sample_id), path(reads_r1), path(reads_r2), path(mob_db)

    output:
      tuple val(sample_id), path("${sample_id}.mash_screen.tsv")

    script:
      """
      mash screen -p ${task.cpus} -i ${params.mashthreshold} ${mob_db}/ncbi_plasmid_full_seqs.fas.msh ${reads_r1} ${reads_r2} | \
        sort -nrk1,1 > ${sample_id}.mash_screen.tsv
      """
}
