#!/usr/bin/env nextflow

import sun.nio.fs.UnixPath

// enable dsl2
nextflow.enable.dsl = 2

include { trim_reads } from './modules/fastp.nf'
include { unicycler } from './modules/unicycler.nf'
include { platon } from './modules/platon.nf'
include { mash_screen } from './modules/mash_screen.nf'
include { quast } from './modules/quast.nf'
include { mob_recon } from './modules/mob_recon.nf'
include { abricate } from './modules/abricate.nf'
include { join_mob_typer_and_abricate_reports } from './modules/join_reports.nf'
include { select_resistance_contigs } from './modules/select_resistance_contigs.nf'
include { select_resistance_reconstructions } from './modules/select_resistance_reconstructions.nf'
include { choose_reference_plasmid_for_sample } from './modules/choose_reference_plasmid_for_sample.nf'
include { get_reference_plasmid } from './modules/get_reference_plasmid.nf'
include { align_reads_to_reference_plasmid } from './modules/align_reads_to_reference_plasmid.nf'
include { calculate_coverage } from './modules/calculate_coverage.nf'
include { call_snps } from './modules/freebayes.nf'

workflow {
  ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  
  ch_mob_db = Channel.fromPath(params.mob_db)
  
  main:
    trim_reads(ch_fastq)

    if (params.pre_assembled) {
      ch_assemblies = Channel.fromPath( params.assembly_search_path ).map{ it -> [it.baseName, it] }.unique{ it -> it[0] }
    } else {
      ch_assemblies = unicycler(trim_reads.out.reads)
    }

    quast(ch_assemblies)

    ch_mash_screen = mash_screen(trim_reads.out.reads.combine(ch_mob_db))

    ch_mob_recon = mob_recon(ch_assemblies.combine(ch_mob_db))

    // pass reconstructed plasmids as [sample_id, [seq1, seq2, seq3...]]
    ch_mob_recon_sequences = mob_recon.out.sequences.map{ it -> [it[0], it[1..-1][0]] }

    ch_abricate = abricate(ch_mob_recon_sequences)

    select_resistance_contigs(mob_recon.out.sequences.join(ch_abricate))

    select_resistance_reconstructions(mob_recon.out.sequences.join(ch_abricate))
    
    ch_join_reports_input = mob_recon.out.mobtyper_reports.cross(ch_abricate).map{ it -> [it[0][0], it[0][1], it[1][1]] }

    ch_combined_abricate_mobtyper_report = join_mob_typer_and_abricate_reports(ch_join_reports_input)

    ch_reference_plasmid_id_for_sample = choose_reference_plasmid_for_sample(ch_combined_abricate_mobtyper_report).splitCsv(header: true, sep: '\t')

    ch_reference_plasmid_for_sample = get_reference_plasmid(ch_reference_plasmid_id_for_sample.map{ it -> [it['sample_id'], it['plasmid_cluster_id'], it['reference_plasmid_id']] }.combine(ch_mob_db))

    align_reads_to_reference_plasmid(trim_reads.out.reads.join(ch_reference_plasmid_for_sample))

    calculate_coverage(align_reads_to_reference_plasmid.out)

    ch_above_coverage_threshold = calculate_coverage.out.filter{ it -> file(it[1]).readLines()[0].split(',')[2].toFloat() > params.min_plasmid_coverage_breadth }.map{ it -> it[0] }

    call_snps(ch_above_coverage_threshold.join(align_reads_to_reference_plasmid.out))
    
}