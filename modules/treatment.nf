#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process clean_trim_galore {

  tag "Processing reads from ${sample_id} with Trim Galore (Cutadapt)"
  publishDir "${params.out_dir}/cutadapt_trimmed_reads/", mode: 'copy'

  input:
  tuple val(sample_id), path(fwd_reads), path(rev_reads)

  output:
  tuple val(sample_id), path("*_val_1.fq.gz"), file("*_val_2.fq.gz"), file("*fastqc.zip")//only works once fastqc.zip for all readpairs has been made

  when:
  params.trim_galore == true

  """
  trim_galore --fastqc \
      --paired ${fwd_reads} ${rev_reads}

  """
}

process clean_fastp {
  //conda "/path/to/your/yaml"

  tag "Processing reads from ${sample_id} with fastp"
  publishDir "${params.out_dir}/fastp_trimmed_reads/", mode: 'copy'

  input:
  tuple val(sample_id), path(fwd_reads), path(rev_reads)

  output:
  tuple val(sample_id), path("*R1.trim.fastq.gz"), file("*R2.trim.fastq.gz"), \
  file("*failed_reads.txt"), file("*fastp.html"), file("*fastp.json")

  """
  fastp -w ${task.cpus} \
  -i ${fwd_reads} \
  -I ${rev_reads} \
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

/*
process annotate_bcf {

  tag "Annotating vcf file of ${sample_id} with BCFtools"
  publishDir "${params.out_dir}/bcft_vcf_annotation", mode: 'copy'
  
  input:
  tuple val(sample_id), file(free_vcf), file(ref)
  
  output:
  tuple val(sample_id), file("*.ann.vcf")
  
  """
bcftools csq --fasta-ref ${ref} --gff-annot ${params.gff3} --phase s --output ${sample_id}.ann.vcf --output-type v ${free_vcf} 
  
  """

}
**/
process var_call_freebayes {

  tag "Calling variants for ${sample_id} with freebayes"
  publishDir "${params.out_dir}/freebayes_variant_calls", mode: 'copy'

  input:
  tuple val(sample_id), path(trim_sort_bam), path(trim_sort_bai), path(ref)

  output:
  tuple val(sample_id), file("*_split.vcf")

  //when:
  //params.freebayes == true

  """
  freebayes \
  -f ${ref} \
  ${params.freeb_flags} \
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

//snpEff download -v MN908947.3 --> goes to /home/jess.cal/.conda/envs/snpeff_env/share/snpeff-5.1-2/./data/
process annotate_snpeff {

  tag "Annotating vcf file of ${sample_id} with SnpEff"
  publishDir "${params.out_dir}/snpeff_vcf_annotation", mode: 'copy'

  input:
  tuple val(sample_id), file(free_vcf)
  
  output:
  tuple val(sample_id), file("*.ann.vcf")
  
  """
  snpEff -noLog\
  -hgvs1LetterAa\
  MN908947.3\
  ${free_vcf} >\
  ${sample_id}.ann.vcf
  """
}


workflow TREATMENT {

  take:
  reads_ch

  main:
  ref_ch = Channel
  .fromPath(params.ref, 
  checkIfExists:true)
  
  if (params.trim_galore) {
      trim_ch=clean_trim_galore(reads_ch)
      trim_ch.collect() | qc_reads
      } else {
      trim_ch=clean_fastp(reads_ch)
      trim_ch.collect() | qc_reads
    }
    
  if (params.bwa) {
     aln_ch = align_bwa(trim_ch
     .map{ it[0,1,2] }
     .combine(ref_ch))
    } else {
    aln_ch = align_minimap2(trim_ch
    .map{ it[0,1,2] }
    .combine(ref_ch))
  }
  
  trim_aln_ch = primer_trim(aln_ch)
  trim_aln_ch | qc_align
  
  //generate mutation table input: 
  trim_aln_ch
  .combine(ref_ch)| var_call_freebayes 
  
  freeb_vcf_ch = var_call_freebayes.out
  
  annotate_snpeff(freeb_vcf_ch)
  
  ann_vcf_ch = annotate_snpeff.out
  
  emit:
  trim_aln_ch
  freeb_vcf_ch
  ann_vcf_ch
}
