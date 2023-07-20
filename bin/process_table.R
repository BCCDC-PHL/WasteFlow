
suppressPackageStartupMessages( library(optparse) )

###############################
message("Checking input given to script and assigning to variables...")
###############################
# command line options provided by user
option_list = list(
  make_option(c("-m", "--mutations"),
              type="character",
              default=NULL,
              help="Snpeff annotated VCF with headers removed",
              metavar="character"),
  make_option(c("-s", "--sample"),
              type="character",
              default=NULL,
              help="Name of sample/library",
              metavar="character"),
  make_option(c("-b", "--bin"),
              type="character",
              default=NULL,
              help="Pathway to functions",
              metavar="character") )

# Converting the provided commands to usable R objects

opt_parser <- OptionParser( option_list = option_list )
opt <- parse_args( opt_parser )

if ( is.null(opt$mutations) | is.null(opt$sample) ){
  stop( "--mutations, --sample are required")
}

message("Running script with the following options:")
message( sprintf("--mutations : %s\n--sample : %s",
                 opt$mutations, opt$sample) )


###############################
message("Loading libraries...")
###############################
suppressPackageStartupMessages( library(readr) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(data.table) )
source( sprintf( "%s/vcf2hgvs.R", opt$bin))

###############################
message("Read and clean mutation table...")
###############################

dt_mutation <- fread( file = opt$mutations,
                      sep = "\t",
                      header = TRUE,
                      na.strings = c("", "NA", "NULL"))

# Fix names
setnames(dt_mutation,
         old = c("ANN[0].ALLELE", "ANN[0].EFFECT", "ANN[0].GENE", "ANN[0].HGVS_P"),
         new = c("ALLELE", "EFFECT", "GENE", "ANNO"))


# Validate content
if ( nrow(dt_mutation) > 0 ) {
  
  ## Positions should be within know SARS-CoV-2 coordinate range 
  pos_min <- min(dt_mutation$POS)
  pos_max <- max(dt_mutation$POS)
  
  if ( any( c(pos_min, pos_max) %between% c(1,29903) == FALSE ) ) {
    stop("Positions not within 1 and 29903")
  }
  
  
  ## ref/alt content should be comprised of A T C and G
  ref_string <- unique( unlist(strsplit(unique(dt_mutation$REF), split = "")) )
  alt_string <- unique( unlist(strsplit(unique(dt_mutation$ALT), split = "")) )
  
  check_ref <- any( ref_string %in% c("A", "C", "G", "T") == FALSE )
  check_alt <- any( alt_string %in% c("A", "C", "G", "T") == FALSE )
  
  if ( check_ref | check_alt ) {
    stop("Ref or alt contain characters outside of A T C and G")
  }
  
  
  ## Normalization should index (seperate multiple alleles and
  # corresponding annotations to seperate rows
  
  if ( nrow(dt_mutation[ ALT != ALLELE ]) > 0 ) {
    stop("Normalization of VCF has failed")
  } 
  
}

# Generate/append variables
# Append sample name to variants acounting for missing/sample type
if ( nrow(dt_mutation) == 0 ) {
  
  dt_mutation <- rbind(dt_mutation, list(SAMPLE = str_replace( opt$sample, "(NEG|POS)", "\\1-") ), fill = TRUE)
  
} else {
  
  dt_mutation[ , SAMPLE := str_replace( opt$sample, "(NEG|POS)", "\\1-") ]
  
}

#[ Seperate name into informative info (time and location of sample)]
#[ Derive collection date (assuming yearmonthday)]
dt_mutation[ , c("LOCATION", "TIME") := tstrsplit(SAMPLE, split = "-", fixed = TRUE, keep=1:2) ][
  , `:=` (COLLECTION_DATE = sprintf("%s-%s-%s",
                                    substr(TIME, start = 1, stop = 4),
                                    substr(TIME, start = 5, stop = 6),
                                    substr(TIME, start = 7, stop = 8)))]



# Clean variables
# The SARS-CoV-2 NCBI reference (https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3)
# states the following 1-based boundaries:
# 5' UTR  1..265
# orf1ab  266..21555
# s       21563..25384
# ORF3a   25393..26220
# E       26245..26472
# M       26523..27191
# ORF6    27202..27387
# ORF7a   27394..27759
# ORF8    27894..28259
# N       28274..29533
# ORF10   29558..29674
# 3' UTR 29675..29903

#[ Assert Gene name based on genomic coordinate to account for intergenic regions ]
#[ Remove Protein prefix if p. annotation is present ]
#[ Strip annotation from synonymous variants (example D783D is not informative) ]
#[ If annotation is missing, construct hgsv format from mutation, position, and REF/ALT ]
#[ Concatenate Gene and annotation ]
#[ If redundant mutations present, sum counts ]
#[ If different mutations lead to same annotation, select most frequent mutation ]
#[ Replace AO with summed AO if <= depth, or replace with depth if SAO > depth ]
#[ Select vars of interest ]

dt_mutation_clean <- dt_mutation[ ! is.na(POS), GENE := fcase( POS %between% c(1, 265), "UTR5",
                                                               POS %between% c(266, 21555), "ORF1ab",
                                                               POS %between% c(21563, 25384), "S",
                                                               POS %between% c(25393, 26220), "ORF3a",
                                                               POS %between% c(26245, 26472), "E",
                                                               POS %between% c(26523, 27191), "M",
                                                               POS %between% c(27202, 27387), "ORF6",
                                                               POS %between% c(27394, 27759), "ORF7a",
                                                               POS %between% c(27894, 28259), "ORF8",
                                                               POS %between% c(28274, 29533), "N",
                                                               POS %between% c(29558, 29674), "ORF10",
                                                               POS %between% c(29675, 29903), "UTR3",
                                                               default = "INTER")][
                                   !is.na(ANNO), ANNO := str_remove(ANNO, pattern = "^p\\.") ][
                                     EFFECT == "synonymous_variant", ANNO :=  NA ][
                                       is.na(ANNO) & !is.na(POS) , ANNO := mapply(extract_hgvs, POS, REF, ALT, TYPE) ][
                                         !is.na(ANNO) , ANNO := sprintf("%s-%s", GENE, ANNO) ][
                                           , count := .N, by = ANNO ][
                                             count > 1, SAO := sum(AO), by = ANNO ][
                                               order(-AO), .SD[1], by = ANNO][
                                                 !is.na(SAO), AO := ifelse( SAO > DP, DP, SAO)][
                                                   , .(SAMPLE, POS, REF, ALT, ANNO, AO, DP, LOCATION, COLLECTION_DATE) ]


#Write to stdout
cat( format_csv( dt_mutation_clean ) )  
