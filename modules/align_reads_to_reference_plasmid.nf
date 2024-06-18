process align_reads_to_reference_plasmid {

    tag { sample_id + " / " + plasmid_id + " / " + resistance_gene }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_${plasmid_id}.sorted{.bam,.bam.bai}", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${reference_plasmid}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2), val(plasmid_id), val(resistance_gene), path(reference_plasmid)

    output:
    tuple val(sample_id), val(plasmid_id), val(resistance_gene), path(reference_plasmid), path("${sample_id}_${plasmid_id}.sorted.bam"), path("${sample_id}_${plasmid_id}.sorted.bam.bai"), emit: alignment
    tuple val(sample_id), val(plasmid_id), path("${sample_id}_${plasmid_id}_coverage.csv"), emit: coverage
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
    // SAM flag 1540 = 'read unmapped' + 'read fails quality checks' + 'read is optical duplicate'
    """
    printf -- "- process_name: align_reads_to_reference_plasmid\\n"          >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "  process_tags:\\n"                                           >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "    ref_plasmid_id: ${plasmid_id}\\n"                         >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "    resistance_gene: ${resistance_gene}\\n"                   >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "  tools:\\n"                                                  >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "    - tool_name: bwa\\n"                                      >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "      tool_version: \$(bwa 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "      subcommand: mem\\n"                                     >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "      parameters:\\n"                                         >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "        - parameter: output_all_alignments\\n"                >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "          value: true\\n"                                     >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "        - parameter: use_soft_clipping_for_supplementary_alignments\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "          value: true\\n"                                     >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "        - parameter: mark_shorter_split_hits_as_secondary\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "          value: true\\n"                                     >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "    - tool_name: samtools\\n"                                 >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "      tool_version: \$(samtools --version | grep 'samtools' | cut -d ' ' -f 2)\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "      subcommand: view\\n"                                    >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "      parameters:\\n" >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "        - parameter: exclude_flags\\n"                        >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml
    printf -- "          value: 1540\\n"                                     >> ${sample_id}_${plasmid_id}_${resistance_gene}_bwa_samtools_provenance.yml

    samtools faidx ${reference_plasmid}

    bwa index ${reference_plasmid}

    bwa mem -a -Y -M -t ${task.cpus} ${reference_plasmid} ${reads_r1} ${reads_r2} | \
        samtools view -h --exclude-flags 1540 | \
        samtools sort -n  | \
        samtools fixmate -m - - | \
        samtools sort | \
        samtools markdup - - > ${sample_id}_${plasmid_id}.sorted.bam
      samtools index ${sample_id}_${plasmid_id}.sorted.bam

      samtools depth -aa ${sample_id}_${plasmid_id}.sorted.bam | calculate_depth.py --sample-id ${sample_id} --plasmid-id ${plasmid_id} > ${sample_id}_${plasmid_id}_coverage.csv
      """
}
