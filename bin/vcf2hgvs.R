##################################################
## Mutation table functions 
## Author: Chad Fibke (chad.fibke@bccdc.ca)
## Last edited: May 25th, 2023
##################################################

##################################################
### extract_hgsv
# Constructs human genome variant standard formats from VCF entried based on mutation type, ref, alt, and position
#
# note, traditional HGVS are here:
# https://varnomen.hgvs.org/recommendations/DNA/variant/substitution/
# note will ignore designation of g. prefix, and assume all mutations
# are reported from a VCF/are genomic
#
## Input:
#   - pos: 1-based position of the first base in the reference feild (from VCF)
#   - ref: Reference base(s) at site of mutation (from VCF)
#   - alt: Alternate base(s) at site of mutation (from VCF)
#   - type: type of mutation (only: "snp", "mnp", "del", "ins", "complex") 
## Output:
#   - character string with following format:
#     - snp/mnp: {POS}{REF}>{ALT} (example 123A>G )
#     - del: {Position(s)_deleted}del (example 123_127del )
#     - ins: {positions_flanking}ins{inserted_sequence} (example 123_124insATCG )
#     - complex: {POS}{REF}>{ALT} (example 100ATG>GGG )
#
## Examples
#
# SNPS/MNP/COMPLEX
# extract_hgvs(pos = 100, ref = "T", alt = "G", type = "snp")       # "100T>G"
# extract_hgvs(pos = 100, ref = "GT", alt = "AA", type = "mnp")     # "100GT>AA"
# extract_hgvs(pos = 100, ref = "G", alt = "ATT", type = "complex") # "100G>ATTcomplex"
#
# Deletions
# extract_hgvs(pos = 100, ref = "TTA", alt = "T", type = "del")     # "101_102del"
# extract_hgvs(pos = 100, ref = "TTA", alt = "TT", type = "del")    # "102del"
# extract_hgvs(pos = 100, ref = "TA", alt = "A", type = "del")      # "100del"
# extract_hgvs(pos = 100, ref = "TAAA", alt = "TA", type = "del")   # "102_103del"
#
# Insertions
# extract_hgvs(pos = 100, ref = "A", alt = "ATGC", type = "ins")    # "100_101insTGC"
# extract_hgvs(pos = 100, ref = "T", alt = "AAACT", type = "ins")   # expect 99_100insTGC
##################################################

extract_hgvs <- function(pos, ref, alt, type){
   
  
  if (type %in% c("snp", "mnp") ){
    ##############################
    # SNP/MNP#####################
    ##############################
    # Assumes: mnps have been split onto seperate lines
  
    if (!nchar(ref) == nchar(alt)) {
      stop("Unexpected snp/mnp format")
    }
    
    hgvc <- paste0( pos, ref, ">", alt)
    
  } else if ( type == "del") {
    ##############################
    # Deletions###################
    ##############################
    # Assumes: deletions are left aligned (ref ="TTA", alt = "T"),
    # where the reference first base is the alt first base. However,
    # this can work if right aligned (very rare) (ref ="TA", alt = "A")
    
    if (nchar(ref) <= nchar(alt)) {
      stop("Unexpected del format")
    }
    
    # Identify deleted sequence and store in regex
    del <- str_remove(ref, alt)
    del_pattern <- paste0("(?=(",del, "))")
    
    # Identify location of remaining nucleotides
    index_kept <- str_locate(ref, alt)
    
    # Find all start and length of matches
    all_matches <- gregexpr(del_pattern, ref, perl=TRUE)
    
    starts <- as.vector( attr(all_matches[[1]],'capture.start') )
    lengths <- as.vector( attr(all_matches[[1]],'capture.length') )
    ends <- starts + lengths - 1
    
    pos_matches <- cbind(starts, ends)
    
    # 1 exact match
    if( nrow(pos_matches) == 1 ){
      
      start <- pos_matches[1,1][[1]]
      end <- pos_matches[1,2][[1]]
    }
    
    # multiple exact matches
    if( nrow(pos_matches) > 1 ){
      
      # collect first deleted match not containing remaining ref sequence
      # to account for repeates (TAAA -> TA, where the deleted AA has
      # multiple matches including a match in the reference portion )
      re_index <- min(which(pos_matches[,1] > index_kept[1,2]))
      start <- pos_matches[re_index,1][[1]]
      end <- pos_matches[re_index,2][[1]]
    }
    
    # Account for single-multi-base deletion format
    if (start == end){
      hgvc <- paste0( pos + (start - 1), "del")
    }
    
    if (start != end){
      hgvc <- paste0( pos + (start - 1), "_", pos + (end - 1), "del")
    }
    
  } else if ( type == "ins" ) {
    ##############################
    # Insertions###################
    ##############################
    # Assumes: insertions are left aligned (ref ="T", alt = "TGTC"),
    # where the reference first base is the alt first base. However,
    # this can work if right aligned (very rare) (ref ="T", alt = "AAACT")
    
    if (nchar(ref) >= nchar(alt)) {
      stop("Unexpected ins format")
    }
    
    # Identify inserted sequence
    insert <- str_remove(alt, ref)
    insert_pattern <- paste0("(?=(",insert, "))")
    
    # Identify location of ref within alt
    index_kept <- str_locate(alt, ref)
    
    match <- gregexpr(insert_pattern, alt, perl=TRUE)
    
    # Find start base left/upstream of the inserted seqence
    start <- as.vector( attr(match[[1]],'capture.start') )[1] - 1
    end <- start + 1
    
    hgvc <- paste0( pos + (start - 1), "_", pos + (end - 1), "ins", insert)
    
  } else if (type %in% c("complex") ){
    ##############################
    # Complex###################
    ##############################
    # No assumptions are made and ref/alt are taken as is
    
    hgvc <- paste0( pos, ref, ">", alt, "complex")
  } else {
    
    hgvc <- NA_character_
  }
  
  return(hgvc)
}
