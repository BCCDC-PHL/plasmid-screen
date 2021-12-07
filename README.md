# plasmid-screen

## Usage

This pipeline is designed to take either raw reads alone, or assemblies plus raw reads as input. If only reads are provided, they will be assembled with [unicycler](https://github.com/rrwick/Unicycler).

```
nextflow run BCCDC-PHL/plasmid-screen \
  --fastq_input </path/to/fastqs> \
  --mob_db </path/to/mob-suite-db> \
  --outdir </path/to/outdir> 
```

If assemblies are already available, they can be provided by adding the `--pre_assembled` flag, and supplying the assemblies to the `--assembly_input` flag.

```
nextflow run BCCDC-PHL/plasmid-screen \
  --pre_assembled \
  --assembly_input </path/to/assemblies> \
  --fastq_input </path/to/fastqs> \
  --mob_db </path/to/mob-suite-db> \
  --outdir </path/to/outdir> 
```

## Outputs

The main output of the pipeline is the 'Resistance plasmid report', which summarizes where the resistance gene was located (contig and position), the quality of the resitance gene match (% identity and
% coverage) and a characterization of the plasmid reconstruction. The report includes the following fields:

```
sample_id
assembly_file
resistance_gene_contig
num_contigs_in_reconstruction
reconstruction_size
resistance_gene
gene_start
gene_end
percent_resistance_gene_coverage
percent_resistance_gene_identity
replicon_types
mob_suite_primary_cluster_id
mob_suite_secondary_cluster_id
mash_nearest_neighbor
mash_neighbor_distance
alignment_ref_plasmid
num_snps_vs_ref_plasmid
```

### Provenance

Each analysis will create a provenance.yml file for each sample. The filename of the `provenance.yml` file includes a timestamp with format `YYYYMMDDHHMMSS` to ensure
that a unique file will be produced if a sample is re-analyzed and outputs are stored to the same directory.

Example provenance output:

```yml
- process_name: mob_recon
  tool_name: mob_recon
  tool_version: 3.0.3
  - parameter: database_directory
    value: /path/to/mob_db
  - parameter: filter_db
    value: /path/to/mob_filter_db
  - parameter: min_con_cov
    value: 95
- process_name: abricate
  tool_name: abricate
  tool_version: 1.0.1
  parameters:
  - parameter: db
    value: ncbi
- input_filename: sample-01_R1.fastq.gz
  input_path: /path/to/sample-01_R1_001.fastq.gz
  sha256: b0534592d61321243897e842a9ea655d396d4496cbf6d926b6c6fea8e06aa98d
- input_filename: sample-01_R2.fastq.gz
  input_path: /path/to/sample-01_R2_001.fastq.gz
  sha256: cc66309103da91e337143eb649196d84ed3ebe2ff08a45b197cd4151d137a167
- input_filename: sample-01.fa
  input_path: /path/to/sample-01.fa
  sha256: 6fffb542711ee301ef1185a403a74fed36c066872e3fbfb7aa5c81464243bd00
- process_name: align_reads_to_reference_plasmid
  tool_name: bwa
  tool_version: 0.7.17-r1188
  parameters:
  - parameter: alignment_algorithm
    value: mem
- process_name: align_reads_to_reference_plasmid
  tool_name: samtools
  tool_version: 1.13
  parameters:
  - parameter: exclude_flags
    value: 1540
- process_name: call_snps
  tool_name: freebayes
  tool_version: 1.3.5
  parameters:
  - parameter: ploidy
    value: 1
  - parameter: min_base_quality
    value: 20
  - parameter: min_mapping_quality
    value: 60
  - parameter: min_coverage
    value: 10
  - parameter: min_alternate_fraction
    value: 0.8
  - parameter: min_repeat_entropy
    value: 1.0
- process_name: call_snps
  tool_name: bcftools
  tool_version: 1.12
- process_name: quast
  tool_name: quast
  tool_version: 5.0.2
- pipeline_name: BCCDC-PHL/plasmid-screen
  pipeline_version: 0.1.0
- timestamp_analysis_start: 2021-12-06T16:12:31.252055
```