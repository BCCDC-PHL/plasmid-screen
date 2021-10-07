process unicycler {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.fa", mode: 'copy'

    input:
      tuple val(sample_id), path(reads_1), path(reads_2)

    output:
      tuple val(sample_id), path("${sample_id}.fa")

    script:
      """
      unicycler --threads ${task.cpus} -1 ${reads_1} -2 ${reads_2} -o ${sample_id}_assembly
      sed 's/^>/>${sample_id}_contig/' ${sample_id}_assembly/assembly.fasta > ${sample_id}.fa
      """
}
