process unicycler {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_unicycler.fa", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_unicycler.fa"), emit: assembly
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: unicycler\\n"  >> ${sample_id}_unicycler_provenance.yml
    printf -- "  tools:\\n"                   >> ${sample_id}_unicycler_provenance.yml
    printf -- "    - tool_name: unicycler\\n" >> ${sample_id}_unicycler_provenance.yml
    printf -- "      tool_version: \$(unicycler --version | cut -d ' ' -f 2 | tr -d 'v')\\n"
    
    unicycler \
	--threads ${task.cpus} \
	-1 ${reads_1} \
	-2 ${reads_2} \
	-o ${sample_id}_assembly

    sed 's/^>/>${sample_id}_contig/' ${sample_id}_assembly/assembly.fasta > ${sample_id}_unicycler.fa
    """
}
