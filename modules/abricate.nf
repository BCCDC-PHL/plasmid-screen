process abricate {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_abricate.tsv", mode: 'copy'

    input:
    tuple val(sample_id),  path(assemblies)

    output:
    tuple val(sample_id), path("${sample_id}_abricate.tsv"), emit: report
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
    db = "ncbi"
    """
    printf -- "- process_name: abricate\\n"  >> ${sample_id}_abricate_provenance.yml
    printf -- "  tools:\\n"                  >> ${sample_id}_abricate_provenance.yml
    printf -- "    - tool_name: abricate\\n" >> ${sample_id}_abricate_provenance.yml
    printf -- "      tool_version: \$(abricate --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_abricate_provenance.yml
    printf -- "      parameters:\\n"         >> ${sample_id}_abricate_provenance.yml
    printf -- "        - parameter: db\\n"   >> ${sample_id}_abricate_provenance.yml
    printf -- "          value: ${db}\n"     >> ${sample_id}_abricate_provenance.yml
    
    abricate \
	--threads ${task.cpus} \
	--db ${db} \
	--nopath \
	${assemblies} > ${sample_id}_abricate.tsv
    """
}
