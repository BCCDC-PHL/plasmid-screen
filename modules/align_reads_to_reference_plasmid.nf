process align_reads_to_reference_plasmid {

    tag { sample_id + " / " + plasmid_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_${plasmid_id}.sorted{.bam,.bam.bai}", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${reference_plasmid}", mode: 'copy'

    input:
      tuple val(sample_id), path(reads_r1), path(reads_r2), val(plasmid_id), path(reference_plasmid)

    output:
      tuple val(sample_id), val(plasmid_id), path(reference_plasmid), path("${sample_id}_${plasmid_id}.sorted.bam"), path("${sample_id}_${plasmid_id}.sorted.bam.bai"), emit: alignment
      tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}_coverage.csv"), emit: coverage
      tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
      // SAM flag 1540 = 'read unmapped' + 'read fails quality checks' + 'read is optical duplicate'
      """
      printf -- "- process_name: align_reads_to_reference_plasmid\\n" > ${sample_id}_bwa_samtools_provenance.yml
      printf -- "  tool_name: bwa\\n  tool_version: \$(bwa 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n" >> ${sample_id}_bwa_samtools_provenance.yml
      printf -- "  parameters:\\n" >> ${sample_id}_bwa_samtools_provenance.yml
      printf -- "  - parameter: alignment_algorithm\\n    value: mem\\n" >> ${sample_id}_bwa_samtools_provenance.yml
      printf -- "- process_name: align_reads_to_reference_plasmid\\n" >> ${sample_id}_bwa_samtools_provenance.yml
      printf -- "  tool_name: samtools\\n  tool_version: \$(samtools --version | grep 'samtools' | cut -d ' ' -f 2)\\n" >> ${sample_id}_bwa_samtools_provenance.yml
      printf -- "  parameters:\\n" >> ${sample_id}_bwa_samtools_provenance.yml
      printf -- "  - parameter: exclude_flags\\n    value: 1540\\n" >> ${sample_id}_bwa_samtools_provenance.yml

      samtools faidx ${reference_plasmid}
      bwa index ${reference_plasmid}
      bwa mem -a -Y -M -t ${task.cpus} ${reference_plasmid} ${reads_r1} ${reads_r2} | \
        samtools view -@ ${task.cpus} -h --exclude-flags 1540 | \
        samtools sort -@ ${task.cpus} -n  | \
        samtools fixmate -@ ${task.cpus} -m - - | \
        samtools sort -@ ${task.cpus} | \
        samtools markdup - - > ${sample_id}_${plasmid_id}.sorted.bam
      samtools index ${sample_id}_${plasmid_id}.sorted.bam

      samtools depth -aa ${sample_id}_${plasmid_id}.sorted.bam | calculate_depth.py --sample-id ${sample_id} --plasmid-id ${plasmid_id} > ${sample_id}_${plasmid_id}_coverage.csv
      """
}
