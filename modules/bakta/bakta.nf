process BAKTA {
    tag { "${meta.id}_${meta.plasmid_id}" }
    
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/50/50b75335f6394ae83fd05f364db27ee2eb75f4170e3525bb2aea47ad717a9e64/data'
        : 'community.wave.seqera.io/library/bakta_diamond:7830b94718da4f96'}"
    
    errorStrategy 'ignore'

    publishDir "${params.outdir}/${meta.id}", pattern: "${meta.id}*.{gbk,gff}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_bakta.gbk"),            emit: gbk
    tuple val(meta), path("*_bakta.gff"),            emit: gff
    tuple val(meta), path("*_bakta.json"),           emit: json
    tuple val(meta), path("*_bakta.log"),            emit: log
    tuple val(meta), path("*_bakta_provenance.yml"), emit: provenance

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.plasmid_id}"
    def locustag = params.use_sample_id_as_annotation_locustag ? "--locust-tag \"${meta.id}\"" : ""
    """
    printf -- "- process_name: bakta\\n"                                     >> ${prefix}_bakta_provenance.yml
    printf -- "  tools:\\n"                                                  >> ${prefix}_bakta_provenance.yml
    printf -- "    - tool_name: bakta\\n"                                    >> ${prefix}_bakta_provenance.yml
    printf -- "      tool_version: \$(bakta --version | cut -d ' ' -f 2)\\n" >> ${prefix}_bakta_provenance.yml
    printf -- "      parameters:\\n"                                         >> ${prefix}_bakta_provenance.yml
    printf -- "        - parameter: --db\\n"                                 >> ${prefix}_bakta_provenance.yml
    printf -- "          value: ${params.bakta_db}\\n"                       >> ${prefix}_bakta_provenance.yml
    printf -- "        - parameter: --keep-contig-headers\\n"                >> ${prefix}_bakta_provenance.yml
    printf -- "          value: null\\n"                                     >> ${prefix}_bakta_provenance.yml

    mkdir tmp

    bakta \
        --threads ${task.cpus} \
        --tmp-dir ./tmp \
        --debug \
        --db ${params.bakta_db} \
        --keep-contig-headers \
        ${locustag} \
        --prefix "${prefix}" \
        ${assembly}

    cp ${prefix}.gff3 ${prefix}_bakta.gff
    cp ${prefix}.gbff ${prefix}_bakta.gbk
    cp ${prefix}.json ${prefix}_bakta.json
    cp ${prefix}.log  ${prefix}_bakta.log
    """
}