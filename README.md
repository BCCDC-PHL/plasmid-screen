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

Alternatively, a 'samplesheet.csv' file may be provided with fields `ID`, `R1`, `R2`:

```csv
ID,R1,R2
sample-01,/path/to/sample-01_R1.fastq.gz,/path/to/sample-01_R2.fastq.gz
sample-02,/path/to/sample-02_R1.fastq.gz,/path/to/sample-02_R2.fastq.gz
...
```

```
nextflow run BCCDC-PHL/plasmid-screen \
  --samplesheet_input </path/to/samplesheet.csv> \
  --mob_db </path/to/mob-suite-db> \
  --outdir </path/to/outdir> 
```

...or if assemblies are available, the `samplesheet.csv` file may also include the field `ASSEMBLY`:

```csv
ID,R1,R2,ASSEMBLY
sample-01,/path/to/sample-01_R1.fastq.gz,/path/to/sample-01_R2.fastq.gz,/path/to/sample-01.fa
sample-02,/path/to/sample-02_R1.fastq.gz,/path/to/sample-02_R2.fastq.gz,/path/to/sample-01.fa
...
```

```
nextflow run BCCDC-PHL/plasmid-screen \
  --pre_assembled \
  --samplesheet_input </path/to/samplesheet.csv> \
  --mob_db </path/to/mob-suite-db> \
  --outdir </path/to/outdir> 
```

## Outputs

The main output of the pipeline is the 'Resistance gene report', which summarizes where the resistance gene was located (contig and position), the quality of the resitance gene match (% identity and
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
depth_coverage_threshold
percent_ref_plasmid_coverage_above_depth_threshold
num_snps_vs_ref_plasmid
```

### Additional Output Files

For each sample, the following output files are created:

```
sample-01/
├── sample-01_20211207163723_provenance.yml
├── sample-01_abricate.tsv
├── sample-01_chromosome.fasta
├── sample-01_fastp.csv
├── sample-01_mash_screen.tsv
├── sample-01_mobtyper_contig_report.tsv
├── sample-01_mobtyper_plasmid_report.tsv
├── sample-01_resistance_gene_report.tsv
├── sample-01_NC_019152.1.snps.vcf
├── sample-01_NC_019152.1.sorted.bam
├── sample-01_NC_019152.1.sorted.bam.bai
├── sample-01_plasmid_AA023.fasta
├── sample-01_plasmid_AA026.fasta
├── sample-01_quast.csv
└── NC_019152.1.fa
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
- process_name: trim_reads
  tool_name: fastp
  tool_version: 0.22.0
  parameters:
  - parameter: cut_tail
    value: true
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
  process_tags:
    ref_plasmid_id: NC_019152.1
  tool_name: bwa
  subcommand: mem
  tool_version: 0.7.17-r1188
  parameters:
  - parameter: alignment_algorithm
    value: mem
- process_name: align_reads_to_reference_plasmid
  process_tags:
    ref_plasmid_id: NC_019152.1
  tool_name: samtools
  subcommand: view
  tool_version: 1.13
  parameters:
  - parameter: exclude_flags
    value: 1540
- process_name: call_snps
  process_tags:
    ref_plasmid_id: NC_019152.1
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
  process_tags:
    ref_plasmid_id: NC_019152.1
  tool_name: bcftools
  subcommand: view
  tool_version: 1.12
  parameters:
  - parameter: include
    value: INFO/TYPE=snp
- process_name: quast
  tool_name: quast
  tool_version: 5.0.2
- pipeline_name: BCCDC-PHL/plasmid-screen
  pipeline_version: 0.1.0
- timestamp_analysis_start: 2021-12-06T16:12:31.252055
```