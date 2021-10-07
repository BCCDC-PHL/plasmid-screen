process calculate_coverage {

    tag { sample_id + " / " + plasmid_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_${plasmid_id}_coverage.csv", mode: 'copy'

    input:
      tuple val(sample_id), val(plasmid_id), path(reference_plasmid), path(alignment), path(alignment_index)

    output:
      tuple val(sample_id), path("${sample_id}_${plasmid_id}_coverage.csv")

    script:
      """
      samtools depth -aa ${alignment} | calculate_depth.py > ${sample_id}_${plasmid_id}_coverage.csv
      """
}