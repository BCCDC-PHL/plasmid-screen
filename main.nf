#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files as hash_files_fastq } from './modules/hash_files.nf'
include { hash_files as hash_files_assemblies } from './modules/hash_files.nf'
include { trim_reads } from './modules/fastp.nf'
include { unicycler } from './modules/unicycler.nf'
include { mash_screen } from './modules/mash_screen.nf'
include { quast } from './modules/quast.nf'
include { parse_quast_report } from './modules/quast.nf'
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
include { join_resistance_plasmid_and_snp_reports } from './modules/join_reports.nf'
include { collect_provenance } from './modules/provenance.nf'
include { pipeline_provenance } from './modules/provenance.nf'

workflow {
  ch_start_time = Channel.of(LocalDateTime.now())
  ch_pipeline_name = Channel.of(workflow.manifest.name)
  ch_pipeline_version = Channel.of(workflow.manifest.version)

  ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

  ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }

  ch_mob_db = Channel.fromPath(params.mob_db)

  main:
    hash_files_fastq(ch_fastq.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))

    trim_reads(ch_fastq)

    if (params.pre_assembled) {
      ch_assemblies = Channel.fromPath( params.assembly_search_path ).map{ it -> [it.baseName.split('_')[0], it] }.unique{ it -> it[0] }
      hash_files_assemblies(ch_assemblies.map{ it -> [it[0], [it[1]]] }.combine(Channel.of("assembly_input")))
    } else {
      unicycler(trim_reads.out.reads)
      ch_assemblies = unicycler.out.assembly
    }

    quast(ch_assemblies)

    parse_quast_report(quast.out.report)

    mash_screen(trim_reads.out.reads.combine(ch_mob_db))

    ch_mob_recon = mob_recon(ch_assemblies)

    // pass reconstructed plasmids as [sample_id, [seq1, seq2, seq3...]]
    ch_mob_recon_sequences = mob_recon.out.sequences.map{ it -> [it[0], it[1..-1][0]] }

    abricate(ch_mob_recon_sequences)

    select_resistance_contigs(mob_recon.out.sequences.join(abricate.out.report))

    select_resistance_reconstructions(mob_recon.out.sequences.join(abricate.out.report))

    ch_join_reports_input = mob_recon.out.mobtyper_reports.cross(abricate.out.report).map{ it -> [it[0][0], it[0][1], it[1][1]] }

    ch_combined_abricate_mobtyper_report = join_mob_typer_and_abricate_reports(ch_join_reports_input)

    ch_reference_plasmid_id_for_sample = choose_reference_plasmid_for_sample(ch_combined_abricate_mobtyper_report).splitCsv(header: true, sep: '\t')

    ch_reference_plasmid_for_sample = get_reference_plasmid(ch_reference_plasmid_id_for_sample.map{ it -> [it['sample_id'], it['plasmid_cluster_id'], it['reference_plasmid_id']] }.combine(ch_mob_db))

    align_reads_to_reference_plasmid(trim_reads.out.reads.join(ch_reference_plasmid_for_sample))

    calculate_coverage(align_reads_to_reference_plasmid.out.alignment)

    ch_above_coverage_threshold = calculate_coverage.out.filter{ it -> file(it[1]).readLines()[1].split(',')[2].toFloat() > params.min_plasmid_coverage_breadth }.map{ it -> it[0] }

    call_snps(ch_above_coverage_threshold.join(align_reads_to_reference_plasmid.out.alignment))

    join_resistance_plasmid_and_snp_reports(ch_combined_abricate_mobtyper_report.join(call_snps.out.num_snps))

    ch_provenance = mob_recon.out.provenance
    ch_provenance = ch_provenance.join(abricate.out.provenance).map{ it -> [it[0], [it[1]] << it[2]] }
    ch_provenance = ch_provenance.join(hash_files_fastq.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    if (params.pre_assembled) {
      ch_provenance = ch_provenance.join(hash_files_assemblies.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }
    ch_provenance = ch_provenance.join(align_reads_to_reference_plasmid.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(call_snps.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(quast.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(ch_fastq.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }
    collect_provenance(ch_provenance)
}
