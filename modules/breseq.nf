process breseq {

    tag { sample_id + " / " + plasmid_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_${plasmid_id}*{gd,json,txt,csv}", mode: 'copy'

    input:
      tuple val(sample_id), path(reads_1), path(reads_2), val(plasmid_id), path(reference_plasmid)

    output:
      tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}_breseq_cnv.gd"), emit: cnv_gd
      tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}_num_cnvs.csv"), emit: num_cnv
      tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}_breseq_summary.json"), emit: summary_json
      tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}_breseq_log.txt"), emit: log

    script:
      """
      breseq --num-processors ${task.cpus} -r ${reference_plasmid} ${reads_1} ${reads_2} --cnv
      cp *_copy_number_variation/*.cn_evidence.gd ${sample_id}_${plasmid_id}_breseq_cnv.gd
      cp output/summary.json ${sample_id}_${plasmid_id}_breseq_summary.json
      cp output/log.txt ${sample_id}_${plasmid_id}_breseq_log.txt
      echo 'sample_id,alignment_ref,num_snps' > ${sample_id}_${plasmid_id}_num_cnvs.csv
      echo '${sample_id},${reference_plasmid},' \$(grep -v '^#' ${sample_id}_${plasmid_id}_breseq_cnv.gd | wc -l) | tr -d ' ' >> ${sample_id}_${plasmid_id}_num_cnvs.csv
      """
}
