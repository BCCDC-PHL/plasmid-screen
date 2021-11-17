process call_snps {

  tag { sample_id + ' / ' + plasmid_id }

  publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_${plasmid_id}.snps.vcf", mode: 'copy'

  input:
    tuple val(sample_id), val(plasmid_id), path(reference_plasmid), path(alignment), path(alignment_index)

  output:
    tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}.snps.vcf")

  script:
  """
  freebayes -f ${reference_plasmid} --ploidy 1 --min-base-quality 20 --min-mapping-quality 60 --min-coverage 10 --min-alternate-fraction 0.8 --min-repeat-entropy 1.0 ${alignment} | bcftools view -Ov --types "snps,mnps" > ${sample_id}_${plasmid_id}.snps.vcf
  """
}