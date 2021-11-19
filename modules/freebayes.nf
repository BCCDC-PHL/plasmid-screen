process call_snps {

  tag { sample_id + ' / ' + plasmid_id }

  publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_${plasmid_id}*.{vcf,csv}", mode: 'copy'

  input:
    tuple val(sample_id), val(plasmid_id), path(reference_plasmid), path(alignment), path(alignment_index)

  output:
    tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}.snps.vcf"), emit: vcf
    tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}_num_snps.csv"), emit: num_snps

  script:
  """
  freebayes -f ${reference_plasmid} --ploidy 1 --min-base-quality 20 --min-mapping-quality 60 --min-coverage 10 --min-alternate-fraction 0.8 --min-repeat-entropy 1.0 ${alignment} | bcftools view -Ov --include 'INFO/TYPE="snp"' > ${sample_id}_${plasmid_id}.snps.vcf
  echo 'sample_id,alignment_ref,num_snps' > ${sample_id}_${plasmid_id}_num_snps.csv
  echo '${sample_id},${reference_plasmid},' \$(grep -v '^#' ${sample_id}_${plasmid_id}.snps.vcf | wc -l) | tr -d ' ' >> ${sample_id}_${plasmid_id}_num_snps.csv
  """
}