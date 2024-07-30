
suppressPackageStartupMessages(library(optparse))

###############################
message("Checking input given to script and assigning to variables...")
###############################
# command line options provided by user
option_list = list(
  make_option(c("-c", "--concat_mutations"),
              type="character",
              default=NULL,
              help="All cleaned mutation tables concatenated into 1 file",
              metavar="character"),
  make_option(c("-o", "--outdir"),
              type="character",
              default=NULL,
              help="Output dir to write report to",
              metavar="character")
  )

# Converting the provided commands to usable R objects

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if ( is.null(opt$concat_mutations) | is.null(opt$outdir) ){
  stop("--concat_mutations and --outdir are required")
}

message("Running script with the following options:")
message(sprintf("--concat_mutations : %s\n--outdir : %s",
                opt$concat_mutations, opt$outdir))


###############################
message("Loading libraries...")
###############################
suppressPackageStartupMessages( library(lubridate) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(data.table) )


###############################
message("Read and clean mutation table...")
###############################

setDTthreads(threads = 60)

dt_concat_mutations <- fread(file = opt$concat_mutations,
                             sep = ",",
                             header = TRUE,
                             na.strings = c("", "NA", "NULL"),
                             colClasses=c("SAMPLE" = "character",
                                          "POS" = "integer",
                                          "REF" = "character",
                                          "ALT" = "character",
                                          "ANNO" = "character",
                                          "AO" = "integer",
                                          "DP" = "integer",
                                          "LOCATION"= "character",
                                          "COLLECTION_DATE"= "Date")
)


# Clean up dataset
#[ Remove control samples ]
#[ Reformat collection date to repersent week of sample, and collapse into regions ]
#[ select vars of interest ]
dt_clean <- dt_concat_mutations[ str_detect(LOCATION, pattern = "NEG|POS", negate = TRUE) ][
  , `:=` ( COLLECTION_DATE = floor_date(ymd(COLLECTION_DATE), unit = "week"),
           REGION = fcase(LOCATION %chin% c("Annacis", "NWL", "Iona", "LionsGate", "Lulu"), "Lower Mainland",
                          LOCATION %chin% c("Kamloops", "Kelowna", "Penticton"), "Interior",
                          LOCATION %chin% c("PrinceGeorge"), "Northern",
                          LOCATION %chin% c("CRD", "CVRD", "Nanaimo"), "Island",
                          str_detect(LOCATION, "YVR"), "YVR" ) ) ][
    , .(REGION, LOCATION, COLLECTION_DATE, ANNO, AO, DP) ]


# Summarize by REGION and LOCATION and merge back
#[ Calculate the variant allele frequency as the sum of alignments supporting allele / sum of depth by place/time/mutation ]

dt_location <- dt_clean[ , .( VAF = round( (sum(AO) / sum(DP)), digits= 4 ) ), by = .(LOCATION, COLLECTION_DATE, ANNO) ]

dt_region <- dt_clean[ , .( VAF = round( (sum(AO) / sum(DP)), digits= 4 ), LOCATION = REGION ), by = .(REGION, COLLECTION_DATE, ANNO) ][
  , REGION := NULL ]

# Combine
dt_final <- rbind(dt_location, dt_region)[ , setorder(.SD, COLLECTION_DATE, ANNO, LOCATION)]


# Write to file
if (!file.exists(opt$outdir)) {
  dir.create(opt$outdir)
}

save( dt_final, file = sprintf("%s/MT_WW_dashboard_input.RData", opt$outdir) )
