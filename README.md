# plasmid-screen

## Usage

```
nextflow run BCCDC-PHL/plasmid-screen \
  --fastq_input </path/to/fastqs> \
  --mob_db </path/to/mob-suite-db> \
  --outdir </path/to/outdir> 
```

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
