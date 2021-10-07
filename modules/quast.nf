process quast {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/qc", pattern: "${sample_id}_quast_report.tsv", mode: 'copy'

    input:
      tuple val(sample_id), path(assembly)

    output:
      tuple val(sample_id), path("${sample_id}_quast_report.tsv")

    script:
      """
      quast --threads ${task.cpus} ${assembly}
      echo "assembly	num_contigs_gt_0_bp	num_contigs_gt_1000_bp	num_contigs_gt_5000_bp	num_contigs_gt_10000_bp	num_contigs_gt_25000_bp	num_contigs_gt_50000_bp	total_length_gt_0_bp	total_length_gt_1000_bp	total_length_gt_5000_bp	total_length_gt_10000_bp	total_length_gt_25000_bp	total_length_gt_50000_bp	num_contigs	largest_contig	total_length	gc_percent	N50	N75	L50	L75	num_Ns_per_100_kbp" > ${sample_id}_quast_report.tsv
      tail -n+2 quast_results/latest/transposed_report.tsv >> ${sample_id}_quast_report.tsv
      """
}
