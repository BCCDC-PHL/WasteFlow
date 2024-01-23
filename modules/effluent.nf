#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
  truncate -s 0 ${sample_id}_freyja_vcf.tsv
  fi

  """

}

process lineage_freyja {
    
  tag "Calculating relative viral lineage abundances from ${sample_id} with Freyja"
  publishDir "${params.out_dir}/freyja_individual_lineage_summaries", mode: 'copy'
  //errorStrategy { sample_id.toLowerCase() =~ /neg/ ? 'ignore' : 'terminate' }
  errorStrategy 'ignore'

  input:
  tuple val(sample_id), path(depths), path(vcf)

  output:
  tuple val(sample_id), file("*.tsv"), file("*.yml") 

  """
  freyja demix\
  --depthcutoff ${params.demixdepth} \
  --output ${sample_id}_freyja_lineage_summary.tsv \
  ${vcf} \
  ${depths}
  
  """
  
}

process bootstrap_freyja {

  tag "Bootstrapping lineage prevalences from ${sample_id} with Freyja"
  publishDir "${params.out_dir}/freyja_individual_bootstrapped_lineages", mode: 'copy'
  //errorStrategy { sample_id.toLowerCase() =~ /neg/ ? 'ignore' : 'terminate' }  
  errorStrategy 'ignore'

  input:
  tuple val(sample_id), path(depths), path(vcf)
  
  output:
  tuple val(sample_id), file("*_lineages.csv"), file("*_summarized.csv"), file("*_lineages.yml")
  
  """
  freyja boot \
  ${vcf} \
  ${depths} \
  --nt ${task.cpus} \
  --nb ${params.bootnum} \
  --output_base ${sample_id}_boot \
  --depthcutoff ${params.demixdepth}

  """

}

process barcode_version {

  tag "Logging WasteFlow pipeline metadata"
  publishDir "${params.out_dir}/pipeline_reports", mode: 'copy'

  output:
  file("*.log")
  
  """
  current_time=\$(date +%Y%m%d-%H%M%S)
  echo \$(freyja demix --version) > \${current_time}_WasteFlow.log
  echo -e "Pipeline name: ${workflow.manifest.name}  
  Pipeline version: ${workflow.manifest.version}  
  Nextflow version: ${nextflow.version}  
  Execution start time: ${workflow.start}
  Executed command: ${workflow.commandLine} 
  Data directory: ${params.data_dir}
  Results directory: ${params.out_dir}
  Project directory: ${workflow.projectDir}
  Launch directory: ${workflow.launchDir} 
  Session ID: ${workflow.sessionId} " >> \${current_time}_WasteFlow.log
  
  """

}

process summarize_freyja {
  
  tag "Summarizing relative viral lineage abundances across all samples with Freyja"
  publishDir "${params.out_dir}/freyja_overall_lineage_summary", mode: 'copy'

  input:
  file("*")

  output:
  file("*.tsv")

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

process vcf2table {
  tag "Extract and clean mutations from ${sample_id} annotated VCF file"
  publishDir "${params.out_dir}/mutation_table", mode: 'copy'

  input:
    tuple val(sample_id), path(ann_vcf)

  output:
    path("${sample_id}_mutation_table.csv")

  """
  SnpSift \
    extractFields -s "," -e "NA" \
    $ann_vcf \
    CHROM POS REF ALT AO DP TYPE "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].HGVS_P" > ${sample_id}_table.tsv

  Rscript ${params.bin}/process_table.R \
    --bin=${params.bin} \
    --mutations=${sample_id}_table.tsv \
    --sample=${sample_id} > ${sample_id}_mutation_table.csv
  """

}

process collectTables {
  tag "Generate mutation frequency and collapse into regional estimates"

  input:
    path MUTATION_TABLE_AGGREGATE

  """
  Rscript ${params.bin}/mutation_watchlist.R \
    --concat_mutations=$MUTATION_TABLE_AGGREGATE \
    --outdir="${params.mut_dir}"
  """
}


workflow EFFLUENT {

  take:
  trim_aln_ch
  freeb_vcf_ch
  ann_vcf_ch
  ref_ch
  
  main:
  
  if (params.rerun_lins != null){
  
  past_depth_vcf_ch = Channel
  .fromPath( params.rerun_lins,
   checkIfExists: true )
  .map{ tuple( it.baseName.split('_freyja')[0], it ) }
  .groupTuple(sort: true)
  .map{tuple(it[0], it[1][0], it[1][1])}
  .filter { it[2].size()>0 }
   
   if (params.freebayes){
    aln_depth_ch = trim_aln_ch
    aln_depth_ch
    .combine(ref_ch)|read_depths
    
    to_date_ch = read_depths.out
    .filter { it[1].size()>0 }
    .join(freeb_vcf_ch
    .filter { it[1].size()>0 })
    .concat(past_depth_vcf_ch)
   }else{
   trim_aln_ch
   .combine(ref_ch)| var_call_freyja
   depth_ch=var_call_freyja
   .out
   .filter { it[2].size()>0 }
   
   to_date_ch = depth_ch
   .concat(past_depth_vcf_ch)}

   to_date_ch | lineage_freyja
  
  }else if (params.freebayes){
  
    aln_depth_ch = trim_aln_ch
  
    aln_depth_ch
    .combine(ref_ch)|read_depths
  
    lineage_freyja(read_depths.out
      .filter { it[1].size()>0 }
      .join(freeb_vcf_ch
      .filter { it[1].size()>0 }))
      
      if (params.boot){
        bootstrap_freyja(read_depths.out
      .filter { it[1].size()>0 }
      .join(freeb_vcf_ch
      .filter { it[1].size()>0 }))   
      }
      
  }else {
    
    trim_aln_ch
    .combine(ref_ch)| var_call_freyja
    
    depth_ch=var_call_freyja
    .out
    .filter { it[2].size()>0 }
    
    lineage_freyja(depth_ch)
    
    if (params.boot){
        bootstrap_freyja(depth_ch)
      }
  
  }
  
  lineage_ch = lineage_freyja.out
  
  lineage_ch.collect() | summarize_freyja

  barcode_version()
  
  if (params.annotate_snps){
  // Extract mutations from VCF and clean entries
  vcf2table( ann_vcf_ch )
  }

  // Collect and clean all previous mutation table
  if ( params.rerun_mut != null ) {
      // Run with current and previous data
      previous_ch = Channel
                          .fromPath( params.rerun_mut,
                                     checkIfExists: true )

      previous_ch
                .concat( vcf2table( ann_vcf_ch ) )
                .collectFile( name: "mutation_tables.csv",
                              keepHeader: true,
                              sort: { it.simpleName } )
                .set { mutation_table_ch }

      // Clean and collapse mutations
      collectTables( mutation_table_ch )
  }
  
  emit:
  summarize_freyja.out
}
