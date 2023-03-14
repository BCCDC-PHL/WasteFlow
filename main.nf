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
 --trim-galore            Include this switch to use trim-galore to trim adapters/ low qual bases [default: fastp]
 --bwa                    Include this switch to use bwa-mem to align reads to reference [minimap2]
 --conda_cache            User-defined location to save conda env [/path/to/conda/env/location/]
 --adapters               Additional adapters to include during trimming with fastp
 --primers                Bed file containing primer scheme for trimming with iVar
 --ivar_flags             Additional options to pass to iVar during primer/ quality trimming
 --freebayes              Include this switch to use freebayes to call variants [Freyja::iVar]
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

"""

/**
---------------------------------------------------------------------------------
process definition
---------------------------------------------------------------------------------
*/

process empty_samps {

tag "Removing empty fastq files from analysis"

input:
tuple val(sample_id), path(reads)

output:
tuple val(sample_id), path(reads)

"""

if [ \$(grep -c "@" ${reads[0]}) -lt 1 ] || [ \$(grep -c "@" ${reads[1]}) -lt 1 ]
then
  mv *${sample_id}* ${params.data_dir}excluded/
fi

"""
}

process clean_trim_galore {

  tag "Processing reads from ${sample_id} with Trim Galore (Cutadapt)"
  publishDir "${params.out_dir}/cutadapt_trimmed_reads/", mode: 'copy'

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("*_val_1.fq.gz"), file("*_val_2.fq.gz"), file("*fastqc.zip")//only works once fastqc.zip for all readpairs has been made

  when:
  params.trim_galore == true

  """
  trim_galore --fastqc \
      --paired ${reads[0]} ${reads[1]}

  """
}

process clean_fastp {
  //conda "/path/to/your/yaml"

  tag "Processing reads from ${sample_id} with fastp"
  publishDir "${params.out_dir}/fastp_trimmed_reads/", mode: 'copy'

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("*R1.trim.fastq.gz"), file("*R2.trim.fastq.gz"), \
  file("*failed_reads.txt"), file("*fastp.html"), file("*fastp.json")

  """
  fastp -w ${task.cpus} \
  -i ${reads[0]} \
  -I ${reads[1]} \
  -o ${sample_id}_R1.trim.fastq.gz \
  -O ${sample_id}_R2.trim.fastq.gz \
  --detect_adapter_for_pe \
  --failed_out ${sample_id}_failed_reads.txt \
  --html ${sample_id}_fastp.html \
  --json ${sample_id}_fastp.json

  """
}


process qc_reads {
  
  tag "Generating QC summary for processed reads with MultiQC"
  publishDir "${params.out_dir}/QC_reports/reads", mode: 'copy'

  input:
  file("*")

  output:
  tuple file("multiqc_report.html"), file("multiqc_data")

  """
  multiqc ${params.out_dir}/.

  """
}


process align_minimap2 {
  
  tag "Mapping reads from ${sample_id} to ${params.ref} with Minimap2"
  publishDir "${params.out_dir}/minimap2_alignments", mode: 'copy'

  input:
  tuple val(sample_id), path(trim_for), path(trim_rev), file(ref)

  output:
  tuple val(sample_id), file("*.sort.bam"), file("*.sort.bam.bai")

  """
  minimap2 \
  -ax sr \
  -t ${task.cpus} \
  ${ref} \
  ${trim_for} \
  ${trim_rev} \
  | samtools sort -@ ${task.cpus} \
  | samtools view -F4 -b -o  ${sample_id}.sort.bam -
  samtools index ${sample_id}.sort.bam

"""
}

process align_bwa {
  tag "Mapping reads from ${sample_id} to ${params.ref} with BWA-MEM"
  publishDir "${params.out_dir}/bwa_alignments", mode: 'copy'

  input:
  tuple val(sample_id), path(trim_for), path(trim_rev), file(ref)

  output:
  tuple val(sample_id), file("*.sort.bam"), file("*.sort.bam.bai")

  when:
  params.bwa == true

  """
  bwa mem -t ${task.cpus} \
  ${ref} \
  ${trim_for} \
  ${trim_rev} \
  | samtools sort -o ${sample_id}.sort.bam \
  | samtools view -F4 -b -o  ${sample_id}.sort.bam -
  samtools index ${sample_id}.sort.bam

"""
}

process primer_trim {
  tag "Trimming primers from reads in ${sample_id} with iVar"
  publishDir "${params.out_dir}/ivar_primer_removed_alignments", mode: 'copy'

  input:
  tuple val(sample_id), path(sort_bam), path(sort_bam_bai)

  output:
  tuple val(sample_id), file("*.sort.bam"), file("*.sort.bam.bai")

  """
  ivar trim \
  -i ${sort_bam} \
  -b ${params.primers} \
  ${params.ivar_flags} \
  -p primer_trim

  samtools sort -o ${sample_id}.primertrim.sort.bam primer_trim.bam
  samtools index ${sample_id}.primertrim.sort.bam
  """
}

process qc_align {
  
  tag "Generating alignment summary metrics for ${sample_id} with SAMtools"
  publishDir "${params.out_dir}/QC_reports/alignments", mode: 'copy'

  input:
  tuple val(sample_id), path(trim_sort_bam), path(trim_sort_bai)

  output:
  tuple file("*.mpileup"), file("*idxstats")

  """
  samtools mpileup \
  ${trim_sort_bam} \
  -o ${sample_id}.mpileup

  samtools idxstats \
  ${trim_sort_bam} \
  > ${sample_id}.idxstats
  """
}

process read_depths {
  
  tag "Extracting read depth coverage in ${sample_id} with SAMtools"
  publishDir "${params.out_dir}/samtools_depths", mode: 'copy'

  input:
  tuple val(sample_id), file(trim_sort_bam), path(trim_sort_bai), file(ref)

  output:
  tuple val(sample_id), file("*_depths.txt")

  """
  samtools mpileup \
  -aa \
  -A \
  -d 600000 \
  -Q 20 \
  -q 0 \
  -B \
  -f ${ref} \
  ${trim_sort_bam} \
  | tee >(cut -f1-4 > ${sample_id}_depths.txt)
  """
}

process var_call_freebayes {

  tag "Calling variants for ${sample_id} with freebayes"
  publishDir "${params.out_dir}/freebayes_variant_calls", mode: 'copy'

  input:
  tuple val(sample_id), path(trim_sort_bam), path(trim_sort_bai), path(ref)

  output:
  tuple val(sample_id), file("*_split.vcf")

  when:
  params.freebayes == true

  """
  freebayes \
  -f ${ref} \
  ${trim_sort_bam} \
  --vcf ${sample_id}.vcf 
  
  #split multiallelic sites for freyja lineage
  bgzip -f ${sample_id}.vcf
  tabix -f -p vcf ${sample_id}.vcf.gz
  bcftools norm -m- ${sample_id}.vcf.gz > \
  ${sample_id}_split.vcf 
  
  #hard-coded - needs to be adusted to take from ref fa header
  if [[ \$(grep -c "MN908947.3" ${sample_id}_split.vcf) -eq 1 ]];
  then
  echo "" > ${sample_id}_split.vcf
  fi


 # sed -i "s/0.5,0.5/0.5/g" ${sample_id}.vcf
 # sed -i "s/0,0/0/g" ${sample_id}.vcf
 # sed -i "s/0.5,0/0.5/g" ${sample_id}.vcf
 # sed -i "s/0,0.5/0.5/g" ${sample_id}.vcf

  #temp sol'n as freyja lineage call gets tripped by AF present in 50:50
  
  """
}
//2372Pool7dilx1_S8_L001: Freebayes produces a file, but there are no variants.. 
//this causes the pipeline to fail. For some reason, the depth file is produced even tho there is 0 coverage..
//so this sample is not being filtered out, as the pipeline would usually do. 

process var_call_freyja {
  
  tag "Calling variants and extracting read depth from ${sample_id} with Freyja (iVar)"
  publishDir "${params.out_dir}/freyja_variant_calls_depths", mode: 'copy'

  input:
  tuple val(sample_id), path(trim_sort_bam),path(trim_sort_bai), path(ref)

  output:
  tuple val(sample_id), file("*.txt"), file("*.tsv")

  """
  freyja variants \
  ${trim_sort_bam} \
  --variants ${sample_id}_freyja_vcf \
  --depths ${sample_id}_freyja_depths.txt \
  --ref ${ref}

  if [[ \$(grep -c "MN908947.3" ${sample_id}_freyja_vcf.tsv) -eq 0 ]];
  then
  echo "" > ${sample_id}_freyja_vcf.tsv
  fi

  """

}

process lineage_freyja {
    
  tag "Calculating relative viral lineage abundances from ${sample_id} with Freyja"
  publishDir "${params.out_dir}/freyja_individual_lineage_summaries", mode: 'copy'

  input:
  tuple val(sample_id), path(depths), path(vcf)

  output:
  tuple val(sample_id), file("*.tsv")

  """
  freyja demix\
  ${vcf} \
  ${depths} \
  --output ${sample_id}_freyja_lineage_summary.tsv
  """

}

process summarize_freyja {
  
  tag "Summarizing relative viral lineage abundances across all samples with Freyja"
  publishDir "${params.out_dir}/freyja_overall_lineage_summary", mode: 'copy'

  input:
  file("*")

  output:
 // tuple file("*.tsv"), file("*.pdf")
  file("*.tsv")
 
 when:
  params.summarize==true

  """
  freyja aggregate \
  ${params.out_dir}/freyja_individual_lineage_summaries/ \
  --output freyja_lineage_summary.tsv

#  freyja plot \
#  freyja_lineage_summary.tsv \
#  ${params.freyja_plot} \
#  --output freyja_lineage_summary.pdf
  """

}

/*
process voc_snv {
  tag "Determining frequency of lineage-specific SNVs in samples"
  publishDir "${params.out_dir}", mode: 'copy'
  //DIY script to come
  input:

  output:

  """

  """

}
/**
---------------------------------------------------------------------------------
workflow
---------------------------------------------------------------------------------
*/

workflow {
  reads_ch = Channel 
  .fromFilePairs(params.pe_reads)
  //.filter { !it[2].empty } //files null from here already
  //.collect() //reads_ch originally when omitting empty samps process

  ref_ch = Channel
  .fromPath(params.ref, checkIfExists:true)

  //reads_ch = empty_samps(read_ch)
  
  if (params.trim_galore) {
    trim_ch=clean_trim_galore(reads_ch)
    trim_ch.collect() | qc_reads
      } else {
    trim_ch=clean_fastp(reads_ch)
    trim_ch.collect() | qc_reads
  }

  //trim_ch.map{ it[0,1,2] }.view()
//align_minimap2(trim_ch.map{ it[0,1,2] }.combine(ref_ch))

if (params.bwa) {
    aln_ch = align_bwa(trim_ch
      .map{ it[0,1,2] }
      .combine(ref_ch))
  } else {
    aln_ch = align_minimap2(trim_ch
      .map{ it[0,1,2] }
      .combine(ref_ch))
  }

//trim_aln_ch = aln_ch
trim_aln_ch = primer_trim(aln_ch)
trim_aln_ch | qc_align
//aln_ch | qc_align

if (params.freebayes){
 trim_aln_ch
  .combine(ref_ch)| var_call_freebayes

  aln_depth_ch = trim_aln_ch

  aln_depth_ch
  .combine(ref_ch)|read_depths

  lineage_freyja(read_depths.out
    .filter { it[1].size()>0 }
    .join(var_call_freebayes
    .out
    .filter { it[1].size()>0 }))
} else {  
  trim_aln_ch
  .combine(ref_ch)| var_call_freyja
  depth_ch=var_call_freyja
  .out
  .filter { it[2].size()>0 && it[1].size()>0 }
  lineage_freyja(depth_ch)
}

lineage_ch = lineage_freyja.out

if (params.summarize){
  lineage_ch.collect() | summarize_freyja

}
//lineage_freyja(trim_aln_ch
//.combine(ref_ch)
//.combine(var_call.out))





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

