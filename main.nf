#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files as hash_files_fastq } from './modules/hash_files.nf'
include { hash_files as hash_files_assemblies } from './modules/hash_files.nf'
include { trim_reads } from './modules/fastp.nf'
include { fastp_json_to_csv } from './modules/fastp.nf'
include { unicycler } from './modules/unicycler.nf'
include { mash_screen } from './modules/mash_screen.nf'
include { quast } from './modules/quast.nf'
include { parse_quast_report } from './modules/quast.nf'
include { mob_recon } from './modules/mob_recon.nf'
include { abricate } from './modules/abricate.nf'
include { join_mob_typer_and_abricate_reports } from './modules/join_reports.nf'
include { select_resistance_chromosomes } from './modules/join_reports.nf'
include { select_resistance_contigs } from './modules/select_resistance_contigs.nf'
include { select_resistance_reconstructions } from './modules/select_resistance_reconstructions.nf'
include { choose_reference_plasmids } from './modules/choose_reference_plasmids.nf'
include { get_reference_plasmid } from './modules/get_reference_plasmid.nf'
include { align_reads_to_reference_plasmid } from './modules/align_reads_to_reference_plasmid.nf'
include { calculate_coverage } from './modules/calculate_coverage.nf'
include { call_snps } from './modules/freebayes.nf'
include { join_resistance_plasmid_and_snp_reports } from './modules/join_reports.nf'
include { concatenate_resistance_reports } from './modules/join_reports.nf'
include { collect_provenance } from './modules/provenance.nf'
include { pipeline_provenance } from './modules/provenance.nf'

workflow {
  ch_start_time = Channel.of(LocalDateTime.now())
  ch_pipeline_name = Channel.of(workflow.manifest.name)
  ch_pipeline_version = Channel.of(workflow.manifest.version)

  ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
  } else {
    ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  }

  ch_mob_db = Channel.fromPath(params.mob_db)

  main:
    hash_files_fastq(ch_fastq.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))

    trim_reads(ch_fastq)

    fastp_json_to_csv(trim_reads.out.json)

    if (params.pre_assembled) {
      if (params.samplesheet_input != 'NO_FILE') {
        ch_assemblies = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['ASSEMBLY']] }
      } else {
        ch_assemblies = Channel.fromPath( params.assembly_search_path ).map{ it -> [it.baseName.split('_')[0], it] }.unique{ it -> it[0] }
      }
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

    // select_resistance_contigs(mob_recon.out.sequences.join(abricate.out.report))

    // select_resistance_reconstructions(mob_recon.out.sequences.join(abricate.out.report))

    ch_join_reports_input = mob_recon.out.mobtyper_reports.cross(abricate.out.report).map{ it -> [it[0][0], it[0][1], it[1][1]] }

    ch_combined_abricate_mobtyper_report = join_mob_typer_and_abricate_reports(ch_join_reports_input)

    select_resistance_chromosomes(ch_combined_abricate_mobtyper_report)
    
    ch_reference_plasmid_id = choose_reference_plasmids(ch_combined_abricate_mobtyper_report).map{ it -> file(it).text }.splitCsv(header: true).map{ it -> [it.sample_id, it.reference_plasmid_id] }
  
    ch_reference_plasmid = get_reference_plasmid(ch_reference_plasmid_id.combine(ch_mob_db))

    align_reads_to_reference_plasmid(trim_reads.out.reads.cross(ch_reference_plasmid).map{ it -> it[0] + it[1].drop(1) })

    // ch_above_coverage_threshold = align_reads_to_reference_plasmid.out.coverage.filter{ it -> file(it[2]).readLines()[1].split(',')[5].toFloat() > params.min_plasmid_coverage_breadth }.map{ it -> [it[0], it[1]] }

    call_snps(align_reads_to_reference_plasmid.out.alignment)

    join_resistance_plasmid_and_snp_reports(ch_combined_abricate_mobtyper_report.cross(call_snps.out.num_snps.join(align_reads_to_reference_plasmid.out.coverage, by: [0, 1])).map{ it -> it[0] + it[1].drop(1) })
    concatenate_resistance_reports(select_resistance_chromosomes.out.join(join_resistance_plasmid_and_snp_reports.out, remainder: true).groupTuple().map{ it -> [it[0], (it[1] - null) + (it[2] - null)] })


    ch_provenance = mob_recon.out.provenance
    ch_provenance = ch_provenance.join(abricate.out.provenance).map{ it -> [it[0], [it[1]] << it[2]] }
    ch_provenance = ch_provenance.join(trim_reads.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(hash_files_fastq.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    if (params.pre_assembled) {
      ch_provenance = ch_provenance.join(hash_files_assemblies.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }
    ch_provenance = ch_provenance.join(align_reads_to_reference_plasmid.out.provenance, remainder: true).map{ it -> it.collect{ x -> x ? x : [] }}.map{ it -> [it[0], it[1] << it[2]] }.groupTuple().map{ it -> [it[0], it[1].flatten()] }
    ch_provenance = ch_provenance.join(call_snps.out.provenance, remainder: true).map{ it -> it.collect{ x -> x ? x : [] }}.map{ it -> [it[0], it[1] << it[2]] }.groupTuple().map{ it -> [it[0], it[1].flatten()] }
    ch_provenance = ch_provenance.join(quast.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(ch_fastq.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }
    collect_provenance(ch_provenance)
}
