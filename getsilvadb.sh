#!/usr/bin/env bash
# This script downloads a desired version of SILVA and then generates 
# a usearch11 UDB database file as well as converting the taxonomy in 
# the fasta headers to be compatible with usearch and SINTAX for classification

# This BASH script is based on the template from https://github.com/kasperskytte/bash_template
# License is MIT, which means the author takes no responsibilities, but you can use it for anything you want

#exit when a command fails (use "|| :" after a command to allow it to fail)
set -o errexit

#exit when a pipe fails
set -o pipefail

#disallow undeclared variables
set -o nounset

#disallow clobbering (overwriting) of files
set -o noclobber

#print exactly what gets executed (useful for debugging)
#set -o xtrace

#default error message if bad usage
usageError() {
  echo "Invalid usage: $1" 1>&2
  echo ""
  eval "bash $0 -h"
}

#default settings
#use all logical cores except 2 unless adjusted by user
MAX_THREADS=${MAX_THREADS:-$(($(nproc)-2))}
VERSION="1.0"
output="refdatabases/"

#fetch and check options provided by user
#flags for required options, checked after getopts loop
r_flag=0
while getopts ":r:o:t:hv" opt; do
case ${opt} in
  h )
    echo "This script downloads a desired release version of the SILVA database and makes it ready for AutoTax."
    echo "Version: $VERSION"
    echo "Options:"
    echo "  -h    Display this help text and exit."
    echo "  -r    (required) The desired SILVA release version, fx \"138.1\"."
    echo "  -o    Output folder. (Default: refdatabases/)"
    echo "  -t    Max number of threads to use. (Default: all available except 2)"
    echo "  -v    Print version and exit."
    exit 1
    ;;
  r )
    SILVArelease="$OPTARG"
    r_flag=1
    ;;
  o )
    output="$OPTARG"
    ;;
  t )
    MAX_THREADS=$OPTARG
    ;;
  v )
    echo "Version: $VERSION"
    exit 0
    ;;
  \? )
    usageError "Invalid Option: -$OPTARG"
    exit 1
    ;;
  : )
    usageError "Option -$OPTARG requires an argument"
    exit 1
    ;;
esac
done
shift $((OPTIND -1)) #reset option pointer

#check all required options
if [ $r_flag -eq 0 ]
then
	usageError "option -r is required"
	exit 1
fi

#function to add timestamps to progress messages
scriptMessage() {
  #check user arguments
  if [ ! $# -eq 1 ]
  then
    echo "Error: function must be passed exactly 1 argument" >&2
    exit 1
  fi
  echo " *** [$(date '+%Y-%m-%d %H:%M:%S')] script message: $1"
}

##### START OF SCRIPT #####
echo "#################################################"
echo "Script version:         $VERSION"
echo "SILVA release version:  $SILVArelease"
echo "Output folder:          $output"
echo "Max. threads to use:    $MAX_THREADS"
echo "#################################################"
echo

mkdir -p "$output"
pushd "$output" 1> /dev/null

scriptMessage "Downloading SILVA release version ${SILVArelease} (NR99 FASTA sequences and metadata)"
SILVAfastafile="SILVA_${SILVArelease}_SSURef_NR99_tax_silva"
SILVAmetadatafile="SILVA_${SILVArelease}_SSURef.full_metadata"
wget -q https://www.arb-silva.de/fileadmin/silva_databases/release_${SILVArelease}/Exports/${SILVAfastafile}.fasta.gz
wget -q https://www.arb-silva.de/fileadmin/silva_databases/release_${SILVArelease}/Exports/full_metadata/${SILVAmetadatafile}.gz

scriptMessage "Unpacking files"
gunzip ${SILVAfastafile}.fasta.gz ${SILVAmetadatafile}.gz

scriptMessage "Extracting typestrain accession ID's from the database and converting taxonomy strings in FASTA headers to SINTAX format"
R --slave --args "$SILVAfastafile" "$SILVAmetadatafile" "$MAX_THREADS" << 'rscript'
#!/usr/local/bin/Rscript

#extract passed args from shell script
args <- commandArgs(trailingOnly = TRUE)
SILVAfastafile <- args[[1]]
SILVAmetadatafile <- args[[2]]
MAX_THREADS <- args[[3]]

#load R packages
suppressPackageStartupMessages({
    require("data.table")
    require("Biostrings")
    require("stringi")
    require("tidyr")
    require("doParallel")
})

data.table::setDTthreads(MAX_THREADS)

#read metadata file and extract typestrains accession ID's (those with [T] flags)
metadata <- fread(SILVAmetadatafile, select = c("acc", "flags"))
typestrains_accIDs <- metadata[stri_detect_regex(flags, "\\[t\\]", opts_regex = stri_opts_regex(case_insensitive = TRUE)), acc]
writeLines(typestrains_accIDs, "typestrains_accessionIDs.txt")

#read SILVA sequences and extract sequence headers
SILVA <- readBStringSet(paste0(SILVAfastafile, ".fasta"))
SILVA_names <- names(SILVA)

#uncomment to filter eukaryotes
#SILVA <- SILVA[!grepl("eukaryota", tolower(SILVA_names)))]

#extract sequence ID's and taxonomy strings from sequence headers into a 2-col table
names <- do.call(rbind, stri_split_fixed(SILVA_names, " ", n = 2))
namesDT <- data.table(ID = names[,1], tax = names[,2])

#extract taxonomy into separate columns
taxCols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sintax <- separate(
  namesDT,
  col = "tax",
  into = taxCols,
  sep = ";")

#cleanup
sintax[
  ,
  (taxCols) := lapply(.SD, function(x) {
    x <- stringi::stri_replace_all_regex(x, ",|;", ".")
    x[stringi::stri_detect_regex(tolower(x), "uncultured|unknown|unidentified|incertae sedis|metagenome|\\bbacterium\\b|\\bpossible\\b")] <- NA
    return(x)
    }),
  .SDcols = taxCols]
sintax[, Kingdom := ifelse(is.na(Kingdom), NA, paste0("d:", Kingdom))]
sintax[, Phylum := ifelse(is.na(Phylum), NA, paste0("p:", Phylum))]
sintax[, Class := ifelse(is.na(Class), NA, paste0("c:", Class))]
sintax[, Order := ifelse(is.na(Order), NA, paste0("o:", Order))]
sintax[, Family := ifelse(is.na(Family), NA, paste0("f:", Family))]
sintax[, Genus := ifelse(is.na(Genus), NA, paste0("g:", Genus))]
sintax[, Species := ifelse(is.na(Species), NA, paste0("s:", Species))]

#stitch together a string in SINTAX format
sintax_header <- paste0(sintax[,ID], ";tax=", unite(sintax[,-1], col = "tax", sep = ",", na.rm = TRUE)[,tax], ";")

#rename sequences and write out
names(SILVA) <- sintax_header
writeXStringSet(SILVA, paste0(SILVAfastafile, "_sintax.fasta"))
rscript

scriptMessage "Making UDB database from FASTA file"
usearch11 -makeudb_usearch ${SILVAfastafile}.fasta \
  -output ${SILVAfastafile}.udb

scriptMessage "Making UDB database from FASTA file with SINTAX compatible headers"
usearch11 -makeudb_usearch ${SILVAfastafile}_sintax.fasta \
  -output ${SILVAfastafile}_sintax.udb

scriptMessage "Extracting typestrain sequences from the database"
usearch11 -fastx_getseqs ${SILVAfastafile}.fasta \
  -labels typestrains_accessionIDs.txt \
  -fastaout ${SILVAfastafile}_typestrains.fasta \
  -label_substr_match \
  -threads ${MAX_THREADS}

scriptMessage "Making UDB database from typestrains FASTA file"
usearch11 -makeudb_usearch ${SILVAfastafile}_typestrains.fasta \
  -output ${SILVAfastafile}_typestrains.udb
##### END OF SCRIPT #####

#print elapsed time since script was invoked
duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
scriptMessage "Done in: $duration!"
popd

echo "As the ARB database file has a date in it's name, it's impossible to guess the filename of future releases of SILVA"
echo "Therefore head to: https://www.arb-silva.de/no_cache/download/archive/release_${SILVArelease}/ARB_files/"
echo "and download the SILVA_${SILVArelease}_SSURef_NR99_dd_mm_yy_opt.arb.gz file manually"
echo