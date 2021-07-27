process align_reads_to_reference_plasmid {

    tag { sample_id + " / " + plasmid_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_${plasmid_id}.sorted{.bam,.bam.bai}", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${reference_plasmid}", mode: 'copy'
    
    cpus 8

    input:
      tuple val(sample_id), path(reads_r1), path(reads_r2), val(plasmid_id), path(reference_plasmid)

    output:
      tuple val(sample_id), val(plasmid_id), path(reference_plasmid), path("${sample_id}_${plasmid_id}.sorted.bam"), path("${sample_id}_${plasmid_id}.sorted.bam.bai")

    script:
      // SAM flag 1540 = 'read unmapped' + 'read fails quality checks' + 'read is optical duplicate'
      """
      samtools faidx ${reference_plasmid}
      bwa index ${reference_plasmid}
      bwa mem -a -Y -M -t ${task.cpus} ${reference_plasmid} ${reads_r1} ${reads_r2} | \
        samtools view -@ ${task.cpus} -h -F 1540 | \
        samtools sort -@ ${task.cpus} -n  | \
        samtools fixmate -@ ${task.cpus} -m - - | \
        samtools sort -@ ${task.cpus} | \
        samtools markdup - - > ${sample_id}_${plasmid_id}.sorted.bam
      samtools index ${sample_id}_${plasmid_id}.sorted.bam
      """
}
