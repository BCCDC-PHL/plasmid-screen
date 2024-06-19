process quast {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_quast.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_quast.csv"), emit: csv
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: quast\\n" >> ${sample_id}_quast_provenance.yml
    printf -- "  tools:\\n"              >> ${sample_id}_quast_provenance.yml
    printf -- "    - tool_name: quast\\n"    >> ${sample_id}_quast_provenance.yml
    printf -- "      tool_version: \$(quast --version | cut -d ' ' -f 2 | tr -d 'v')\\n" >> ${sample_id}_quast_provenance.yml
    
    quast --threads ${task.cpus} ${assembly}

    cp quast_results/latest/transposed_report.tsv ${sample_id}_quast.tsv

    parse_quast_report.py ${sample_id}_quast.tsv > ${sample_id}_quast.csv
    """
}
