process select_resistance_reconstructions {

    tag { sample_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "selected/${sample_id}*.fasta", mode: 'copy'

    cpus 1

    executor 'local'

    input:
      tuple val(sample_id),  path(assembly), path(abricate_report)

    output:
      tuple val(sample_id), path("selected/${sample_id}*.fasta") optional true

    script:
      """
      mkdir selected
      grep --no-filename 'CARBAPENEM' ${abricate_report} | cut -f 1 > selected.tsv || true
      while read -r selected;
      do
        mv \${selected} selected
      done < selected.tsv
      """
}
