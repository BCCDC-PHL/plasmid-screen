#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { trimReads } from './modules/trim_galore.nf'
include { assemble } from './modules/unicycler.nf'
include { mob_recon } from './modules/mob_recon.nf'
include { abricate } from './modules/abricate.nf'
include { join_mob_typer_and_abricate_reports } from './modules/join_reports.nf'

workflow {
  ch_fastq = Channel.fromFilePairs( "${params.fastq_input}/*_{R1,R2}_*.fastq.gz" )
  ch_abricate_dbs = Channel.of(["ncbi"])
  
  main:
    trimReads(ch_fastq)
    assemble(trimReads.out)
    mob_recon(assemble.out)
    
    // This takes this:              [sample_id, [sample_id_chromosome.fasta, sample_id_plasmid.fasta, sample_id_mob_typer_result.txt]]
    // ...and transforms it to this: [sample_id, sample_id_chromosome.fasta, ncbi]
    //                               [sample_id, sample_id_plasmid.fasta, ncbi]
    ch_abricate_input = mob_recon.out.map{ sample_id, files -> files.collect{ [sample_id, it] } }.flatMap{ sampleid_files -> sampleid_files.findAll{ it[1].name =~ /fasta/ }}.combine(ch_abricate_dbs)
    abricate(ch_abricate_input)
    ch_mob_typer_reports = mob_recon.out.map{ sample_id, files -> files.collect{ [sample_id, it] } }.flatMap{ sampleid_files -> sampleid_files.findAll{ it[1].name =~ /mobtyper_results/ }}
    ch_join_reports_input = ch_mob_typer_reports.cross(abricate.out).map{ it -> [it[0][0], it[0][1], it[1][1]] }
    join_mob_typer_and_abricate_reports(ch_join_reports_input).map{ it -> it[1] }collectFile(keepHeader: true, sort: { it.text }, name: "combined_abricate_mobtyper_report.tsv", storeDir: "${params.outdir}")
}