process trimReads {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sample_id, path(forward), path(reverse))
    * @output trimgalore_out tuple(sample_id, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sample_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    cpus 2

    input:
      tuple val(grouping_key), path(reads)

    output:
      tuple val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz")

    script:
      if (grouping_key =~ '_S[0-9]+_') {
        sample_id = grouping_key.split("_S[0-9]+_")[0]
      } else if (grouping_key =~ '_') {
        sample_id = grouping_key.split("_")[0]
      } else {
        sample_id = grouping_key
      }
      """
      trim_galore --cores 2 --paired ${reads}
      """
}