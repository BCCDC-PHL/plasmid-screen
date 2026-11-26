#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { PLASMID_ANNOTATION } from '../subworkflows/local/plasmid_annotation.nf'


include { hash_files as hash_files_fastq }          from '../modules/hash_files.nf'
include { hash_files as hash_files_assemblies }     from '../modules/hash_files.nf'
include { fastp }                                   from '../modules/fastp.nf'
include { fastplong }                               from '../modules/fastplong/fastplong.nf'
include { unicycler }                               from '../modules/unicycler.nf'
include { mash_screen }                             from '../modules/mash_screen.nf'
include { quast }                                   from '../modules/quast.nf'
include { mob_recon }                               from '../modules/mob_recon.nf'
include { abricate as abricate_ncbi }               from '../modules/abricate.nf'
include { abricate as abricate_plasmidfinder }      from '../modules/abricate.nf'
include { join_mob_typer_and_abricate_reports }     from '../modules/join_reports.nf'
include { select_resistance_chromosomes }           from '../modules/join_reports.nf'
include { select_resistance_contigs }               from '../modules/select_resistance_contigs.nf'
include { select_resistance_reconstructions }       from '../modules/select_resistance_reconstructions.nf'
include { choose_reference_plasmids }               from '../modules/choose_reference_plasmids.nf'
include { get_reference_plasmid }                   from '../modules/get_reference_plasmid.nf'
include { align_short_reads_to_reference_plasmid }  from '../modules/align_short_reads_to_reference_plasmid.nf'
include { MINIMAP2_ALIGN }                          from '../modules/minimap2_align/minimap2_align.nf'
include { calculate_coverage }                      from '../modules/calculate_coverage.nf'
include { call_snps }                               from '../modules/freebayes.nf'
include { join_resistance_plasmid_and_snp_reports } from '../modules/join_reports.nf'
include { concatenate_resistance_reports }          from '../modules/join_reports.nf'
include { collect_provenance }                      from '../modules/provenance.nf'
include { pipeline_provenance }                     from '../modules/provenance.nf'

workflow PLASMID_SCREEN {
    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_workflow_metadata = Channel.value([
        workflow.sessionId,
        workflow.runName,
        workflow.manifest.name,
        workflow.manifest.version,
        workflow.start,
    ])

    ch_provenance = channel.empty()
    ch_assembly_provenance = channel.empty()

    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)
    if (!(params.mode in ['short_reads', 'long_reads'])) {
        error "Invalid --mode '${params.mode}'. Choose 'short_reads' or 'long_reads'."
    }

    ch_samplesheet = samplesheet.map { row ->
        tuple(
            [id: row[0].id.toString()],
            row[1] ? file(row[1]) : null,
            row[2] ? file(row[2]) : null,
            row[3] ? file(row[3]) : null,
            row[4] ? file(row[4]) : null
        )
    }
    

    ch_sample_assemblies = ch_samplesheet
        .filter { row -> row[3] }
        .map { row -> tuple(row[0].id, row[3]) }

    if (params.mode == 'short_reads') {
        ch_fastq_short = ch_samplesheet.map { row ->
            if (!row[1] || !row[2]) {
                throw new IllegalArgumentException("Short-read mode requires fastq_1 and fastq_2 for sample ${row[0].id}")
            }
            tuple(row[0].id, row[1], row[2])
        }
        ch_fastq_long = Channel.empty()
        ch_short_without_assembly = ch_samplesheet
            .filter { row -> !row[3] }
            .map { row -> tuple(row[0].id, true) }
    } else {
        ch_fastq_short = Channel.empty()
        ch_fastq_long = ch_samplesheet.map { row ->
            if (!row[3] || !row[4]) {
                throw new IllegalArgumentException("Long-read mode requires assembly and fastq for sample ${row[0].id}")
            }
            tuple(row[0].id, row[4])
        }
        ch_short_without_assembly = Channel.empty()
    }

    ch_mob_db = Channel.fromPath(params.mob_db)

    ch_input_samples = ch_fastq_short.map { it -> it[0] }
        .mix(ch_fastq_long.map { it -> it[0] })

    ch_fastq_hash_inputs = ch_fastq_short
        .map { it -> [it[0], it[1..-1]] }
        .mix(ch_fastq_long.map { it -> [it[0], [it[1]]] })
    hash_files_fastq(ch_fastq_hash_inputs.combine(Channel.of("fastq_input")))
    hash_files_assemblies(ch_sample_assemblies.map { it -> [it[0], [it[1]]] }.combine(Channel.of("assembly_input")))


    

    if (params.mode == 'long_reads') {
        fastplong(ch_fastq_long)
        ch_read_provenance = fastplong.out.provenance
        ch_assemblies = ch_sample_assemblies
        ch_reads = fastplong.out.filtered_reads
        
    } else {
        fastp(ch_fastq_short)
        ch_read_provenance = fastp.out.provenance
        ch_fastq_short_for_assembly = fastp.out.reads
            .join(ch_short_without_assembly)
            .map { it -> tuple(it[0], it[1], it[2]) }
        unicycler(ch_fastq_short_for_assembly)
        ch_assembly_provenance = hash_files_assemblies.out.provenance.mix(unicycler.out.provenance)
        ch_assemblies = ch_sample_assemblies.mix(unicycler.out.assembly)
        ch_reads = fastp.out.reads
    }

    quast(ch_assemblies)
    ch_mash_reads = ch_reads.map { row -> tuple(row[0], row[1..-1]) }
    mash_screen(ch_mash_reads.combine(ch_mob_db))
    ch_mob_recon = mob_recon(ch_assemblies)

    // pass reconstructed plasmids as [sample_id, [seq1, seq2, seq3...]]
    ch_mob_recon_sequences = mob_recon.out.plasmids.map { it -> [it[0], it[1..-1][0]] }
    
    abricate_ncbi(ch_mob_recon_sequences.combine(Channel.of("ncbi")))
    abricate_plasmidfinder(ch_mob_recon_sequences.combine(Channel.of("plasmidfinder")))

    ch_join_reports_input = mob_recon.out.mobtyper_reports.cross(abricate_ncbi.out.report)
        .map { it -> [it[0][0], it[0][1], it[0][2], it[1][1]] }
    ch_combined_abricate_mobtyper_report = join_mob_typer_and_abricate_reports(ch_join_reports_input)
    select_resistance_chromosomes(ch_combined_abricate_mobtyper_report)

    ch_reference_plasmid_id = choose_reference_plasmids(ch_combined_abricate_mobtyper_report)
        .map { it -> file(it).text }
        .splitCsv(header: true)
        .map { it -> [it.sample_id, it.reference_plasmid_id, it.resistance_gene_id] }
    ch_reference_plasmid = get_reference_plasmid(ch_reference_plasmid_id.combine(ch_mob_db))
    
    
    if (params.mode == 'long_reads'){
        MINIMAP2_ALIGN(ch_reads.join(ch_reference_plasmid, by: 0))
        ch_alignment = MINIMAP2_ALIGN.out.alignment
        ch_coverage = MINIMAP2_ALIGN.out.coverage
    }
    else {
        align_short_reads_to_reference_plasmid(ch_reads.join(ch_reference_plasmid, by: 0))
        ch_alignment = align_short_reads_to_reference_plasmid.out.alignment
        ch_coverage = align_short_reads_to_reference_plasmid.out.coverage
    }
    call_snps(ch_alignment)
    
    join_resistance_plasmid_and_snp_reports(ch_combined_abricate_mobtyper_report.cross(call_snps.out.num_snps.join(ch_coverage, by: [0, 1])).map{ it -> it[0].toList() + it[1].drop(1) })
    
    concatenate_resistance_reports(select_resistance_chromosomes.out.join(join_resistance_plasmid_and_snp_reports.out, remainder: true).groupTuple().map{ it -> [it[0], (it[1] - null) + (it[2] - null)]}.map{ it -> [it[0]] + [it[1].unique{ path -> path.getFileName() }] })

    plasmids_annotation_ch = mob_recon.out.plasmids.flatMap { sample_id, fasta_list ->
        def fasta_files = fasta_list instanceof List ? fasta_list : fasta_list ? [fasta_list] : []
        fasta_files.collect { fasta ->
            def fasta_name = fasta.getFileName().toString()
            def plasmid_id = fasta_name
                .replaceFirst("^${sample_id}_plasmid_", "")
                .replaceFirst(/\.[^.]+$/, "")
            tuple([id: sample_id, plasmid_id: plasmid_id], fasta)
        }
    }
    PLASMID_ANNOTATION(plasmids_annotation_ch)
    plasmid_annotation_provenance = PLASMID_ANNOTATION.out.plasmid_annotation_provenance
        .map { meta, provenance -> tuple(meta.id, provenance) }
        .groupTuple()

    if (params.collect_outputs) {
        if (params.mode == 'short_reads') {
            fastp.out.csv.map { it -> it[1] }.collectFile(
                name: params.collected_outputs_prefix + "_fastp.csv",
                storeDir: params.outdir,
                keepHeader: true,
                sort: true,
            )
        }
        quast.out.csv.map { it -> it[1] }.collectFile(
            name: params.collected_outputs_prefix + "_quast.csv",
            storeDir: params.outdir,
            keepHeader: true,
            sort: true,
        )
        concatenate_resistance_reports.out.map { it -> it[1] }.collectFile(
            name: params.collected_outputs_prefix + "_resistance_gene_report.tsv",
            storeDir: params.outdir,
            keepHeader: true,
            sort: true,
        )
    }
    
    // Collect per-sample provenance from the input and analysis processes.
    ch_provenance = ch_input_samples.combine(ch_pipeline_provenance).map { it -> [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_files_fastq.out.provenance).map { it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(ch_read_provenance, remainder: true)
        .map { it -> [it[0], it[1] + (it[2] ? [it[2]] : [])] }
    ch_provenance = ch_provenance.join(ch_assembly_provenance, remainder: true).map { it -> [it[0], it[1] + (it[2] ? [it[2]] : [])] }
    
    ch_provenance = ch_provenance.join(quast.out.provenance).map { it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(mash_screen.out.provenance, remainder: true).map { it -> [it[0], it[1] + (it[2] ? [it[2]] : [])] }
    ch_provenance = ch_provenance.join(mob_recon.out.provenance).map { it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(abricate_ncbi.out.provenance).map { it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(abricate_plasmidfinder.out.provenance).map { it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(plasmid_annotation_provenance, remainder: true)
        .map { it -> [it[0], (it[1] ?: []) + (it[2] ?: [])] }
    if (params.mode == 'short_reads') {
        ch_provenance = ch_provenance.join(align_short_reads_to_reference_plasmid.out.provenance, remainder: true)
            .map { it -> it.collect { value -> value ? value : [] } }
            .map { it -> [it[0], it[1] << it[2]] }
            .groupTuple()
            .map { it -> [it[0], it[1].flatten()] }
        ch_provenance = ch_provenance.join(call_snps.out.provenance, remainder: true)
            .map { it -> it.collect { value -> value ? value : [] } }
            .map { it -> [it[0], it[1] << it[2]] }
            .groupTuple()
            .map { it -> [it[0], it[1].flatten()] }
    }
    ch_provenance = ch_provenance.map { [it[0], it[1].unique { path -> path.getFileName() }] }
    collect_provenance(ch_provenance)
}
