process MINIMAP2_ALIGN {

    tag { sample_id + " / " + plasmid_id + " / " + resistance_gene }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data' :
        'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce' }"

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_${plasmid_id}.sorted{.bam,.bam.bai}", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${reference_plasmid}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(plasmid_id), val(resistance_gene), path(reference_plasmid)

    output:
    tuple val(sample_id), val(plasmid_id), val(resistance_gene), path(reference_plasmid), path("${sample_id}_${plasmid_id}.sorted.bam"), path("${sample_id}_${plasmid_id}.sorted.bam.bai"), emit: alignment
    tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}_coverage.csv"), emit: coverage
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: MINIMAP2_ALIGN:align_long_reads_to_reference_plasmid\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "  process_tags:\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "    ref_plasmid_id: ${plasmid_id}\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "    resistance_gene: ${resistance_gene}\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "  tools:\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "    - tool_name: minimap2\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "      tool_version: \$(minimap2 --version)\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "      parameters:\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "        - parameter: preset\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "          value: map-ont\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "    - tool_name: samtools\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml
    printf -- "      tool_version: \$(samtools --version | grep 'samtools' | cut -d ' ' -f 2)\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_minimap2_samtools_provenance.yml

    minimap2 -ax map-ont -t ${task.cpus} ${reference_plasmid} ${reads} | \\
      samtools sort -@ ${task.cpus} -o ${sample_id}_${plasmid_id}.sorted.bam -
    samtools index -@ ${task.cpus} ${sample_id}_${plasmid_id}.sorted.bam
    samtools depth -aa ${sample_id}_${plasmid_id}.sorted.bam | calculate_depth.py --sample-id ${sample_id} --plasmid-id ${plasmid_id} > ${sample_id}_${plasmid_id}_coverage.csv
    """
}
