process assemble {
    /**
    * Assembles reads using unicycler (https://github.com/rrwick/Unicycler)
    * @input tuple val(sample_id), path(forward), path(reverse)
    * @output 
    */

    tag { sample_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_assembly.fasta", mode: 'copy'

    cpus 16

    input:
      tuple val(sample_id), path(forward), path(reverse)

    output:
      tuple val(sample_id), path("${sample_id}_assembly.fasta")

    script:
      """
      unicycler --threads ${task.cpus} -1 ${forward} -2 ${reverse} -o .
      cp assembly.fasta ${sample_id}_assembly.fasta
      """
}
