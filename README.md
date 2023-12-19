![image](pics/wasteflow_logo.png)

## Introduction

The excretion of pathogens by humans during active infection allows the contemporaneous monitoring of outbreaks in wastewater.
Wastewater includes liquids drained from toilets, showers, and dishwashers. Sewage is treated across 5 wastewater treatment plants in
Metro Vancouver before becoming effluent. Samples are collected prior to treatment from the non-solid fraction, RNA is extracted, and genomic material
is sequenced. This pipeline, WasteFlow, written by JMC in Nextflow aims to provide a standardized approach to further resolve presence or absence of SARS-CoV-2 in wastewater, into lineages infecting the population to inform precautionary measures and therapeutic strategies.

## Table of Contents

- [Introduction](#introduction)
- [Quick-Start Guide](#quick-start%guide)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Input](#input)
- [Output](#output)
- [Workflow](#workflow)
- [References](#references)

## Quick-Start Guide

```
nextflow run main.nf -profile conda --data_dir /path/to/pe/fastq/files/ --ref ./resources/cov2_ref.fasta --primers ./resources/articV5.3.bed --conda_cache /home/jess.cal/caches/  
```

## Dependencies

[Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is required to build an environment with required workflow dependencies.

This bioinformatic pipeline requires Nextflow:
```
conda install -c bioconda nextflow
```
or download and add the nextflow executable to a location in your user $PATH variable:
```
curl -fsSL get.nextflow.io | bash
mv nextflow ~/bin/
```
Nextflow requires Java v8.0+, so check that it is installed:
```
java -version
```
The OS-independent conda environment activated upon running WasteFlow is specified in the
```wasteflow_env.yml``` file of the project directory and is built when 
```-profile conda``` is included in the command line. Nextflow will save
the environment to the project directory by default. Alternatively, the 
necessary conda environment can be saved to a different shared location 
accesible to compute nodes by adding ```--conda_cache /path/to/new/location/```.

## Installation

To copy the program into a directory of your choice, from desired directory run:
```
git clone https://github.com/j3551ca/WasteFlow.git
cd WasteFlow
nextflow run main.nf -profile conda --data_dir /path/to/input/data
```

## Quick Start

```
# Run all modules from the WasteFlow directory after cloning repository:
nextflow run main.nf -profile conda --data_dir /path/to/input/data

# Run all modules outside of WasteFlow
nextflow run j3551ca/WasteFLow -profile conda --data_dir /path/to/input/data

# Run all modules + collate all *mutation_tables.csv, construct mutation measurements, and write cumulative file
nextflow run main.nf -profile conda --data_dir /path/to/input/data --table_search_string "/path/to/dir/holding/*mutation_table.csv" --sum_dir /path/to/write/cumlative_file  
```

## Input

The pipeline requires the following files:

- Reference genome for guided read alignment [./data_dir/sequences.fasta]
- Paired-end sequencing reads. WasteFlow will accept *.fastq.gz, *.fq.gz, *.fastq, *.fq [./data_dir/*.fq]. 
Ensure the absolute path of the directory containing data to be analyzed is used, otherwise MultiQC will throw an error.
- Primer scheme bed file [./data_dir/primer[vers].bed]

## Output

## Workflow

WasteFlow was designed to allow multiple workflows for various pathogens. For example, reads can be processed with trim-galore or fastp, alignments can be generated with bwa or minimap2, variants can be called with iVar or Freebayes, and lineage calls made by Freyja can be performed on one directory or multiple. The modularity of Nextflow facilitates the addition of alternate packages. WasteFlow has currently been tested on SARS-CoV-2. Common workflows are outlined below: 

--combine_reps This allows multiple sequencing replicates to be combined into one file containing all forward reads and one file containing all reverse reads. 

--annotate_snps This will produce a table of annotated mutations for each sample. Currently variants are called by Freebayes and annotated with SnpEff. 

--primerless_reads 
--primer_pairs 

The default workflow with only data_dir and output_dir specified:

```mermaid
flowchart TB
    subgraph "INPUT"
    v0["*_{R{1,2}_001,R{1,2},{1,2}}*
    *.{fastq,fq,fastq.gz,fq.gz}"]
    v2["ref.fasta"]
    v18["ref.fasta"]
    end
    subgraph TREATMENT
    v3([clean_fastp])
    v5([qc_reads])
    v9([align_minimap2])
    v10([primer_trim])
    v12([qc_align])
    v15([var_call_freebayes])
    v16([annotate_snpeff])
    v1(( ))
    v4["*{R1,R2}.trim.fastq.gz,*failed_reads.txt,
    *fastp.html,*fastp.json"]
    v7["*{R1,R2}.trim.fastq.gz"]
    v11["*{R1,R2}.trim.fastq.gz"]
    v14(( ))
    v19(( ))
    end
    subgraph "OUTPUT"
    v6["multiqc_report.html"]
    v13["*.{depths,covstats}"]
    v17["*.ann.vcf"]
    v25["freyja_lineage_summary.tsv"]
    v27["*_WasteFlow.log"]
    end
    subgraph EFFLUENT
    v20([var_call_freyja])
    v22([lineage_freyja])
    v24([summarize_freyja])
    v26([barcode_version])
    v21["*.{txt,tsv}"]
    v23["*.{tsv,yml}"]
    end
    v0 --> v1
    v2 --> v7
    v2 --> v11
    v2 --> v14
    v1 --> v3
    v3 --> v4
    v3 --> v7
    v4 --> v5
    v5 --> v6
    v7 --> v9
    v9 --> v10
    v10 --> v11
    v10 --> v14
    v10 --> v19
    v11 --> v12
    v12 --> v13
    v14 --> v15
    v15 --> v16
    v16 --> v17
    v18 --> v19
    v19 --> v20
    v20 --> v21
    v21 --> v22
    v22 --> v23
    v23 --> v24
    v24 --> v25
    v26 --> v27
```
## References

