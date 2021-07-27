process assemble {
    /**
    * @input tuple val(sample_id), path(forward), path(reverse)
    * @output 
    */

    tag { sample_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_assembly.fasta", mode: 'copy'

    cpus 16

    input:
      tuple val(sample_id), path(forward), path(reverse)

    output:
      tuple val(sample_id), path("${sample_id}_*.fasta")

    script:
      """
      spades.py --threads ${task.cpus} -1 ${forward} -2 ${reverse} -o .
      sed -r 's/^>NODE_([[:digit:]]+)_(.+)/>${sample_id}_\\1 \\2/' contigs.fasta > ${sample_id}_assembly.fasta
      """
}
