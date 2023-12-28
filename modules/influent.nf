#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
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
**/
process merge_reps {

tag "Merging reads from replicates of ${sample_id}"
publishDir "${params.out_dir}/combined_replicate_reads", mode: 'copy'

input:
tuple val(sample_id), path(fwd_reps), path(rev_reps)

output:
tuple val(sample_id), path("*_cmbnd_R1.fq"), path("*_cmbnd_R2.fq")


"""
  fwds=\$(echo "${fwd_reps}" | sed -r 's/,/ /g' | sed -r 's/\\[/ /g' | sed -r 's/\\]/ /g')
  revs=\$(echo "${rev_reps}" | sed -r 's/,/ /g' | sed -r 's/\\[/ /g' | sed -r 's/\\]/ /g')
  
  zcat \${fwds} > ${sample_id}_cmbnd_R1.fq
  zcat \${revs} > ${sample_id}_cmbnd_R2.fq

"""

}

workflow INFLUENT {

 if(params.combine_reps){
 
   fwd_ch = Channel
  .fromPath(params.pe_reads)
  .filter{ it.name =~ /_R1_/  }
  .map{ tuple( it.baseName.split('-[0-9]{1}-|-COVIDWW-')[0], it ) }
  .groupTuple()
   
  rev_ch = Channel
  .fromPath(params.pe_reads)
  .filter{ it.name =~ /_R2_/  }
  .map{ tuple( it.baseName.split('-[0-9]{1}-|-COVIDWW-')[0], it ) }
  .groupTuple()
  
  merged_ch = fwd_ch
  .join(rev_ch)
  
  merge_reps(merged_ch)
  
  read_tuple_ch = merge_reps
  .out
  
  }else{
  
  read_tuple_ch = Channel
  .fromFilePairs(params.pe_reads)
  .map{tuple(it[0], it[1][0], it[1][1])}
  
  }
  
  emit:
  read_tuple_ch
  
}