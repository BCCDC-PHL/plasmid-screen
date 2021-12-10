process call_snps {

  tag { sample_id + ' / ' + plasmid_id }

  publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_${plasmid_id}*.vcf", mode: 'copy'

  input:
    tuple val(sample_id), val(plasmid_id), path(reference_plasmid), path(alignment), path(alignment_index)

  output:
    tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}.snps.vcf"), emit: vcf
    tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}_num_snps.csv"), emit: num_snps
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance    

  script:
  """
  printf -- "- process_name: call_snps\\n" > ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  tool_name: freebayes\\n  tool_version: \$(freebayes --version 2>&1 | cut -d ' ' -f 3 | tr -d 'v')\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  process_tags:\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "    ref_plasmid_id: ${plasmid_id}\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  parameters:\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  - parameter: ploidy\\n    value: 1\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  - parameter: min_base_quality\\n    value: 20\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  - parameter: min_mapping_quality\\n    value: 60\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  - parameter: min_coverage\\n    value: 10\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  - parameter: min_alternate_fraction\\n    value: 0.8\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  - parameter: min_repeat_entropy\\n    value: 1.0\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "- process_name: call_snps\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  process_tags:\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "    ref_plasmid_id: ${plasmid_id}\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  tool_name: bcftools\\n  tool_version: \$(bcftools --version | grep 'bcftools' | cut -d ' ' -f 2)\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  subcommand: view\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  parameters:\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml
  printf -- "  - parameter: include\\n    value: INFO/TYPE=snp\\n" >> ${sample_id}_${plasmid_id}_freebayes_provenance.yml

  freebayes -f ${reference_plasmid} --ploidy 1 --min-base-quality 20 --min-mapping-quality 60 --min-coverage 10 --min-alternate-fraction 0.8 --min-repeat-entropy 1.0 ${alignment} | bcftools view -Ov --include 'INFO/TYPE="snp"' > ${sample_id}_${plasmid_id}.snps.vcf
  echo 'sample_id,alignment_ref,num_snps' > ${sample_id}_${plasmid_id}_num_snps.csv
  echo '${sample_id},${reference_plasmid},' \$(grep -v '^#' ${sample_id}_${plasmid_id}.snps.vcf | wc -l) | tr -d ' ' >> ${sample_id}_${plasmid_id}_num_snps.csv
  """
}