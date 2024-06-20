process mash_screen {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_mash_screen.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2), path(mob_db)

    output:
    tuple val(sample_id), path("${sample_id}_mash_screen.tsv"), emit: mash_screen
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: mash_screen\\n"                  >> ${sample_id}_mash_provenance.yml
    printf -- "  tools:\\n"                                      >> ${sample_id}_mash_provenance.yml
    printf -- "    - tool_name: mash\\n"                        >> ${sample_id}_mash_provenance.yml
    printf -- "      tool_version: \$(mash --version)\\n"       >> ${sample_id}_mash_provenance.yml
    printf -- "      parameters:\\n"                             >> ${sample_id}_mash_provenance.yml
    printf -- "        - name: threshold\\n"                    >> ${sample_id}_mash_provenance.yml
    printf -- "          value: ${params.mashthreshold}\\n"     >> ${sample_id}_mash_provenance.yml

    mash screen -p ${task.cpus} -i ${params.mashthreshold} ${mob_db}/ncbi_plasmid_full_seqs.fas.msh ${reads_r1} ${reads_r2} | \
        sort -nrk1,1 > ${sample_id}_mash_screen.tsv
    """
}
