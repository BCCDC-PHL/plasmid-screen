process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.{json,csv}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastp.csv"), emit: csv
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: trim_reads\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "  tools:\\n"                   >> ${sample_id}_fastp_provenance.yml
    printf -- "    - tool_name: fastp\\n"     >> ${sample_id}_fastp_provenance.yml
    printf -- "      tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "      parameters:\\n"          >> ${sample_id}_fastp_provenance.yml
    printf -- "        - parameter: cut_tail\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "          value: true\\n"         >> ${sample_id}_fastp_provenance.yml

    fastp \
	--thread ${task.cpus} \
	--in1 ${reads_1} \
	--in2 ${reads_2} \
	--cut_tail \
	--out1 ${sample_id}_trimmed_R1.fastq.gz \
	--out2 ${sample_id}_trimmed_R2.fastq.gz \
	--unpaired1 ${sample_id}_unpaired.fastq.gz \
	--unpaired2 ${sample_id}_unpaired.fastq.gz

    mv fastp.json ${sample_id}_fastp.json

    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
    """
}
