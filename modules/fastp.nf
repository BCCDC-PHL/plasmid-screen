process trim_reads {
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/qc", pattern: '${sample_id}_fastp.json', mode: 'copy'

    input:
      tuple val(sample_id), path(reads_1), path(reads_2)

    output:
      tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads
      tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json

    script:
      """
      fastp --thread ${task.cpus} --in1 ${reads_1} --in2 ${reads_2} --cut_tail --out1 ${sample_id}_trimmed_R1.fastq.gz --out2 ${sample_id}_trimmed_R2.fastq.gz --unpaired1 ${sample_id}_unpaired.fastq.gz --unpaired2 ${sample_id}_unpaired.fastq.gz
      mv fastp.json ${sample_id}_fastp.json
      """
}