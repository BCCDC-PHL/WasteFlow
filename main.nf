#!/usr/bin/env nextflow

//enable domain-specific-language 2
nextflow.enable.dsl=2

/**
---------------------------------------------------------------------------------
function definition
---------------------------------------------------------------------------------
*/

def helpMe() {
  log.info """

Overview:
Nextflow pipeline for detection of pathogen strains in wastewater.


Usage:
nextflow run main.nf -profile conda [OPTIONS]

Mandatory arguments:
 --dir                    User's directory that contains input paired-end sequence reads (fastq files).
                          WasteFlow accepts gzip compressed or uncompressed files.
 --out_dir                User-specified directory to output WasteFlow results to. Default is data_dir/results.                         
 --ref                    Reference genome used to align reads to during guided assembly


Optional arguments:
 --trim_galore            Include this switch to use trim-galore to trim adapters/ low qual bases [default: fastp]
 --bwa                    Include this switch to use bwa-mem to align reads to reference [minimap2]
 --conda_cache            User-defined location to save conda env [/path/to/conda/env/location/]
 --adapters               Additional adapters to include during trimming with fastp
 --primers                Bed file containing primer scheme for trimming with iVar
 --ivar_flags             Additional options to pass to iVar during primer/ quality trimming
 --freebayes              Include this switch to use freebayes to call variants [Freyja::iVar]
 --freeb_flags            Additional options to pass freebayes during variant calling
 --sum_dir                User-defined location to save cumulative mutation table 
 --table_search_string    User-defined pathway/search string used to collect past mutation tables  
 --summarize              Summarize lineage calls across all samples
 --version                Current WasteFlow version number
 --help                   This usage statement
        """
}

def version() {
  log.info """
  WasteFlow version: ${workflow.manifest.version}
  """
}

//displays help upon request
if (params.help) {
  helpMe()
  exit 0 //stop running
}

//version upon request
if (params.version) {
  version()
  exit 0
}

/**
---------------------------------------------------------------------------------
program introduction
---------------------------------------------------------------------------------
*/

// this prints program header with mandatory input
log.info """


░██╗░░░░░░░██╗░█████╗░░██████╗████████╗███████╗███████╗██╗░░░░░░█████╗░░██╗░░░░░░░██╗
░██║░░██╗░░██║██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔════╝██║░░░░░██╔══██╗░██║░░██╗░░██║
░╚██╗████╗██╔╝███████║╚█████╗░░░░██║░░░█████╗░░█████╗░░██║░░░░░██║░░██║░╚██╗████╗██╔╝
░░████╔═████║░██╔══██║░╚═══██╗░░░██║░░░██╔══╝░░██╔══╝░░██║░░░░░██║░░██║░░████╔═████║░
░░╚██╔╝░╚██╔╝░██║░░██║██████╔╝░░░██║░░░███████╗██║░░░░░███████╗╚█████╔╝░░╚██╔╝░╚██╔╝░
░░░╚═╝░░░╚═╝░░╚═╝░░╚═╝╚═════╝░░░░╚═╝░░░╚══════╝╚═╝░░░░░╚══════╝░╚════╝░░░░╚═╝░░░╚═╝░░

\n===================================================================================
data directory: ${params.data_dir}
reference genome: ${params.ref}
results directory: ${params.out_dir}
previous mutation table search string: ${params.table_search_string != null ? params.table_search_string : 'working with current data' }
"""

/**
---------------------------------------------------------------------------------
workflow
---------------------------------------------------------------------------------
*/


include { INFLUENT } from "${projectDir}/modules/influent.nf"
include { TREATMENT } from "${projectDir}/modules/treatment.nf"
include { EFFLUENT } from "${projectDir}/modules/effluent.nf"

workflow {

INFLUENT | TREATMENT | EFFLUENT
 
}

/**
---------------------------------------------------------------------------------
optional notification of completion
---------------------------------------------------------------------------------
*/

/**
workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: '', subject: 'Flu Tree Analysis Complete on Sabin', body: msg)
}



process clean {
  tag "Removing Nextflow work directory?"
    shell:
        "rm -rfv work"
}
*/

