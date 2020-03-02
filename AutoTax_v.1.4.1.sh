#!/bin/bash
VERSION="1.4.0"
##### requirements (tested with) #####
# SILVA database files (both SSURef and typestrains) in arb format (use latest)
# awk+grep+cat (included in most linux distributions)
# GNU parallel (version 20161222)
# usearch (version 11)
# SINA (version 1.6.0)
# R (version 3.6)
#   Note: The required R packages will install automatically, but Bioconductor can cause trouble
#   if it desides to update base R packages like MASS, mgcv, lattice and others.
#   If so run R as root and run: install.packages("BiocManager"); BiocManager::install()

#################################
############# setup #############
#################################

#Load Sina and R (for in-house use only)

module load Sina/1.6.0-rc.1-foss-2018a
module load R-bundle-Bioconductor/3.8-foss-2018a-R-3.5.0

#paths to executables

export sina=$(which sina)
export usearch=$(which usearch11)
export R=$(which R)
export Rscript=$(which Rscript)

# Set paths to the SILVA nr99 database and the typestrain database extracted from SILVA nr99 (see article supplementary for details). 
# .udb files have to be created first from fasta files using the example below:
# usearch11 -makeudb_usearch SILVA_132_SSURef_Nr99_tax_silva.fasta -output SILVA_132_SSURef_Nr99_tax_silva.udb

silva_db="/srv/PHN/users/md/Reference_databases/arb/SILVA_132_SSURef_NR99_13_12_17_opt.arb"
silva_udb="/srv/PHN/users/md/Reference_databases/SILVA_132_SSURef_Nr99_tax_silva.udb"
typestrains_udb="/srv/PHN/users/md/Reference_databases/SILVA_132_SSURef_Nr99_tax_silva_typestrains.udb"

#Define the de novo taxonomy prefix.
#Results will be in the format "prefix_g_123" for a de novo Genus based on ESV number 123

denovo_prefix="denovo"

##########################################################################
########## Check that everything is okay before running AutoTax ##########
##########################################################################

# Script must be run in bash, not shell
if [ ! -n "$BASH" ]
    then
    echo "This script must be run with BASH (bash), not Shell (sh)" 1>&2
    exit 1
fi

#set appropriate error handling
set -o errexit -o pipefail -o noclobber -o nounset


#adds a header to echo, for a better console output overview
echoWithHeader() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')]: $1"
}

userError() {
  self=`basename "$0"`
  echo "Invalid usage: $1" 1>&2
  echo ""
  echo "Run 'bash $self -h' for help"
}

#check if a folder is present and empty
checkFolder() {
  if [ -d $1 ]
      then
          echoWithHeader "A directory named '$1' already exists and is needed for this script to run. Do you want to clear the contents and continue (y/n)?"
          read ANSWER
          if [ $ANSWER = "y" ]
              then 
              echoWithHeader "Clearing the folder '$1'..."
              rm -rf $1/*
              else 
              if [ $ANSWER = "n" ]
                  then
                  echoWithHeader "Exiting script."
                  echo ""
                  exit 0
                  else
                  echoWithHeader "Exiting script."
                  echo ""
                  exit 1
              fi
          fi
      else
      	mkdir $1
  fi
}

#fetch and check options provided by user
while getopts ":hi:d:t:" opt; do
  case ${opt} in
    h )
      echo "Pipeline for extracting Exact Sequence Variants (ESV's) from full length 16S rRNA gene DNA sequences and generating de novo taxonomy"
      echo "Version: $VERSION"
      echo "Options:"
      echo "  -h    Display this help text."
      echo "  -i    Input FASTA file with full length DNA sequences to process (required)."
      echo "  -d    FASTA file with previously processed ESV sequences."
      echo "          ESV's generated from the input sequences will then be appended to this and de novo taxonomy is rerun."
      echo "  -t    Maximum number of threads to use. Default is all available cores except 2."
      exit 1
      ;;
    i )
      DATA=$OPTARG
      ;;
    d )
      ESVDB=$OPTARG
      ;;
    t )
      MAX_THREADS=$OPTARG
      ;;
    \? )
      userError "Invalid Option: -$OPTARG"
      exit 1
      ;;
    : )
      userError "Option -$OPTARG requires an argument"
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

#exit if no data provided
if [ -z "${DATA:-}" ]; then
  userError "No input data provided"
  exit 1
fi

#set threads to a default value if not provided by the user
MAX_THREADS=${MAX_THREADS:-$((`nproc`-2))}

#continue if no errors
echoWithHeader "Running fSSU de novo taxonomy pipeline (max threads: $MAX_THREADS)..."
checkFolder temp
checkFolder output
echoWithHeader "Checking for required R packages and installing if missing..."
#Run R and check for installed packages, install if needed
$R --slave << 'checkRpkgs'
  suppressPackageStartupMessages({
    #Biostrings (and BiocManager which is used to install Biostrings)
    if(!require("Biostrings")) {
      if(!require("BiocManager")) {
        install.packages("BiocManager")
      }
      BiocManager::install("Biostrings", update = FALSE, ask = FALSE)
    }

    #doParallel
    if(!require("doParallel"))
      install.packages("doParallel")
    
    #stringr
    if(!require("stringr"))
      install.packages("stringr")
      
    #stringi
    if(!require("stringi"))
      install.packages("stringi")

    #data.table
    if(!require("data.table"))
      install.packages("data.table")

    #tidyr
    if(!require("tidyr"))
      install.packages("tidyr")

    #dplyr
    if(!require("dplyr"))
      install.packages("dplyr")
    })
checkRpkgs

###################################
########## Generate ESVs ##########
###################################

# Define function to generate ESVs
generateESVs() {
  ##############
  # note: must set threads to 1 in the below usearch commands or the ordering of sequences
  # will be scrambled and inconsistent between individual runs. It doesn't use more than one 
  # thread anyways, so no speedup will be gained.
  
  # Orient all fSSU correctly based on the SILVA database
  echoWithHeader "  - Orienting sequences..."
  $usearch -orient $1 -db $silva_udb -fastaout temp/fSSUs_oriented.fa -threads 1 -quiet
  
  # Dereplicate sequences
  echoWithHeader "  - Dereplicating sequences..."
  $usearch -fastx_uniques temp/fSSUs_oriented.fa -fastaout temp/uniques_wsize.fa -sizeout -minuniquesize 1 -strand plus -relabel preESV -threads 1 -quiet
  
  # Denoise with UNOISE3 accepting sequences seen only twice (see article supplementary for why this is acceptable)
  echoWithHeader "  - Denoising sequences using UNOISE3"
  $usearch -unoise3 temp/uniques_wsize.fa -zotus temp/preESVs.fa -minsize 2
  #cp temp/uniques_wsize.fa temp/preESVs.fa

  # Remove shorter ESVs that also part of longer ESVs
  #To ensure reproducibility, the output ESVs will be ordered first by decending size (times the unique ESV has been observed) and when sizes are identical by preESV ID (ascending). 
  #Rename with new ID's to "ESV(ID).(length)" fx: "ESV1.1413"
  echoWithHeader "  - Finding the longest representative sequence of identical sequences, then reorder and rename..."
  $R --slave --args "$MAX_THREADS" << 'findLongestSortESVsBySizeAndID'
    #extract passed args from shell script
    args <- commandArgs(trailingOnly = TRUE)
    nProc <- as.integer(args[[1]])
    
    #load R packages
    suppressPackageStartupMessages({
      require("Biostrings")
      require("doParallel")
    })
    
    #read sequences
    ESVs <- Biostrings::readBStringSet("temp/preESVs.fa")
    seqs_chr <- as.character(ESVs)
    
    #Find the longest of representative ESVs and remove shorter, but otherwise identical ESVs
    registerDoParallel(cores = nProc)
    removeIDs <- foreach(
      i = seq_along(seqs_chr),
      .combine = c,
      .inorder = TRUE
    ) %dopar% {
      if(any(stringi::stri_detect_fixed(str = seqs_chr[-c(1:i)], pattern = seqs_chr[i])))
        return(i)
    }
    stopImplicitCluster()
    if(!is.null(removeIDs))
      ESVs <- ESVs[-removeIDs]
    
    #rename the sequences and write out
    names(ESVs) <- paste0("ESV", 1:length(ESVs), ".", lengths(ESVs))
    Biostrings::writeXStringSet(ESVs, file = "temp/ESVs.fa")
findLongestSortESVsBySizeAndID
}

# Run the generateESV function

echoWithHeader "Generating ESVs..."
generateESVs $DATA

#if expanding an existing ESV database (-d option), identify redundant new ESVs compared to the previous ESV database and merge the two before continuing

if [ -n "${ESVDB:-}" ]; then
  echoWithHeader "Finding new unique ESVs and adding them to the existing database..."
  cp temp/ESVs.fa output/newESVs.fa
  $R --slave --args "$ESVDB" "$MAX_THREADS" << 'addnewESVs'
    #extract passed args from shell script
    args <- commandArgs(trailingOnly = TRUE)
    ESVDB <- args[[1]]
    nProc <- as.integer(args[[2]])
    
    #load R packages
    suppressPackageStartupMessages({
      require("Biostrings")
      require("doParallel")
      require("stringi")
      require("data.table")
    })
    
    doParallel::registerDoParallel(cores = nProc)
    querySeqs <- Biostrings::readBStringSet(ESVDB)
    targetSeqs <- Biostrings::readBStringSet("temp/ESVs.fa")
    tlengths <- Biostrings::width(targetSeqs)
    
	#Match existing ESVs to the new ESVs to find sequences of equal or shorter length
	res <- foreach::foreach(
	  query = as.character(querySeqs),
	  name = names(querySeqs),
	  .inorder = TRUE,
	  .combine = "rbind") %dopar% {
	    qlength <- Biostrings::width(query)
	    targetSeqs <- targetSeqs[which(tlengths >= qlength)]
	    hitIDs <- stringi::stri_detect_fixed(pattern = query, str = targetSeqs)
	    nHits <- sum(hitIDs)
	    if(nHits == 0) {
	      return(NULL)
	    } else {
	      targetHits <- targetSeqs[hitIDs]
	      data.table::data.table(queryName = rep(name, times = nHits),
	                             queryLength = rep(qlength, times = nHits),
	                             querySeq = rep(as.character(query), times = nHits),
	                             targetName = names(targetHits),
	                             targetLength = Biostrings::width(targetHits),
	                             targetSeq = as.character(targetHits))
	    }
	  }
	doParallel::stopImplicitCluster()

	newESVs <- targetSeqs
	if(length(res) != 0) {
    newESVs <- newESVs[-which(names(targetSeqs) %in% res[["targetName"]])]
    #replace shorter sequences in existing ESV database with new longer ones
    replaceSeqs <- res[targetLength > queryLength, ]
    if(nrow(replaceSeqs) > 0) {
      #there may be more than one match of shorter sequences to longer ones, 
      #if so use the first of the longer ESV (lowest number, thus higher coverage)
      replaceSeqs <- replaceSeqs[, .SD[1], by = queryName]
      #make a column with the exact differences between the identical sequences
      replaceSeqs[,diff := stringi::stri_replace_all_fixed(str = targetSeq, 
                                                           pattern = querySeq, 
                                                           replacement = "..."), by = queryName]
      replaceIDs <- which(names(querySeqs) %in% replaceSeqs[["queryName"]])
      
      #replace shorter with longer
      querySeqs[replaceIDs] <- replaceSeqs[["targetSeq"]]
      #and adjust length annotation in the names of the new ESVs
      names(querySeqs)[replaceIDs] <- paste0(gsub("\\..*$", ".", names(querySeqs)[replaceIDs]), 
                                             Biostrings::width(querySeqs[replaceIDs]))
      
      #write information about the replacements to a log file
      replacementLog <- paste0(replaceSeqs[,queryName], 
                               " has been replaced with ",
                               replaceSeqs[,targetName], 
                               ", difference: ",
                               replaceSeqs[,diff])
      warning(paste0("  - ",
                     length(replaceIDs), 
                     if(length(replaceIDs) > 1) 
                       " ESVs have" 
                     else 
                       " ESV has",
                     " been replaced by a longer ESV, see the logfile \"./output/replacedESVs.log\" for details"),
              call. = FALSE)
      writeLines(replacementLog,
                 "./output/replacedESVs.log")
    } else if(nrow(replaceSeqs) == 0)
      warning("No sequences have been replaced", call. = FALSE)
  }
  #add the new sequences to the database with new names continuing ID numbering
	ESVs <- c(querySeqs, newESVs)
	names(ESVs) <- paste0("ESV", 1:length(ESVs), ".", lengths(ESVs))
	Biostrings::writeXStringSet(ESVs, "temp/ESVs.fa")
addnewESVs
fi

#############################################################################
########## Align and trim ESVs based on the global SILVA alignment ##########
#############################################################################

#Define function to align and trim sequences based on the global SILVA alignment using SINA

sinaAlign() {
  # Preparation
  local DATA=$1
  local OUTPUTID=$2
  local DB=$3
  local SINA_THREADS=${4:-$MAX_THREADS}
  
  $sina -i ${DATA} -o temp/${OUTPUTID}_aligned.fa -r $DB \
  --threads $SINA_THREADS --log-file temp/${OUTPUTID}_log.txt

  #trim sequences and strip alignment gaps
  echoWithHeader "  - Trimming, formatting, and sorting data..."
  awk '!/^>/ {$0=substr($0, 1048, 41788)}1' temp/${OUTPUTID}_aligned.fa > temp/${OUTPUTID}_trimmed.fa
  $usearch -quiet -fasta_stripgaps temp/${OUTPUTID}_trimmed.fa -fastaout temp/tmp.fa \
    && mv temp/tmp.fa temp/${OUTPUTID}_trimmed.fa
  
  #sort sequences and stats files by ESV ID using R
  $R --slave --args temp/${OUTPUTID}_trimmed.fa << 'sortSINAoutput'
	#extract passed args from shell script
	args <- commandArgs(trailingOnly = TRUE)

  #load R packages
  suppressPackageStartupMessages({
    require("Biostrings")
  })
  
  #reorder FASTA
  algn.fa <- Biostrings::readBStringSet(args[[1]])
  algn.fa <- algn.fa[order(as.integer(gsub("[^0-9+$]|\\..*$", "", names(algn.fa))))]
  Biostrings::writeXStringSet(x = algn.fa, 
                              filepath = paste0(tools::file_path_sans_ext(args[[1]]),
                                                "_sorted.", 
                                                tools::file_ext(args[[1]])))
sortSINAoutput
}

#Run the sinaAlign function
echoWithHeader "Aligning ESVs with SILVA database using SINA..."
sinaAlign temp/ESVs.fa ESVs_SILVA $silva_db $MAX_THREADS

#########################################################################
########## Map ESVs against the SILVA and typestrain databases ##########
#########################################################################

# Define functions to map ESVs against the SILVA and type strain database. For SILVA we find the top hit, 
# whereas for the type strain database we find all references with >=98.7% identity 
searchTaxDB() {
  IN=$1
  DB=$2
  OUT=$3
  $usearch -usearch_global $IN -db $DB -maxaccepts 0 -maxrejects 0 -top_hit_only -strand plus -id 0 -blast6out $OUT -threads $MAX_THREADS
}

searchTaxDB_typestrain() {
  IN=$1
  DB=$2
  OUT=$3
  $usearch -usearch_global $IN -db $DB -maxaccepts 0 -maxrejects 0 -strand plus -id 0.987 -blast6out $OUT -threads $MAX_THREADS
}

#SILVA
echoWithHeader "Finding taxonomy of best hit in SILVA database..."
searchTaxDB temp/ESVs_SILVA_trimmed_sorted.fa $silva_udb temp/tax_SILVA.txt

#Typestrains
echoWithHeader "Finding the taxonomy of species within the 98.7% threshold in the typestrains database..."
searchTaxDB_typestrain temp/ESVs_SILVA_trimmed_sorted.fa $typestrains_udb temp/tax_typestrains.txt

#######################################################################################
########## Perform clustering for the de novo taxonomy based on trimmed ESVs ##########
#######################################################################################

#assign with identity thresholds based on Yarza et al, 2014 using cluster_smallmem (no multithread support) to preserve order of input sequences.
echoWithHeader "Generating de novo taxonomy..."
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.987 -maxrejects 0 -uc temp/SILVA_ESV-S.txt -centroids temp/SILVA_ESV-S_centroids.fa -sortedby other
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.945 -maxrejects 0 -uc temp/SILVA_S-G.txt -centroids temp/SILVA_S-G_centroids.fa -sortedby other
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.865 -maxrejects 0 -uc temp/SILVA_G-F.txt -centroids temp/SILVA_G-F_centroids.fa -sortedby other
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.82 -maxrejects 0 -uc temp/SILVA_F-O.txt -centroids temp/SILVA_F-O_centroids.fa -sortedby other
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.785 -maxrejects 0 -uc temp/SILVA_O-C.txt -centroids temp/SILVA_O-C_centroids.fa -sortedby other
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.75 -maxrejects 0 -uc temp/SILVA_C-P.txt -centroids temp/SILVA_C-P_centroids.fa -sortedby other


#############################################################################
########## Format, merge, and ultimately generate AutoTax taxonomy ##########
#############################################################################

echoWithHeader "Merging and reformatting taxonomy..."

# Export the taxonomy assignment Rscript, so that it can be rerun independently if needed.
cat << 'generatedenovotax' > temp/Rscript.R

# 
#!/usr/bin/Rscript
#load R packages
suppressPackageStartupMessages({
  require("stringr")
  require("tidyr")
  require("dplyr")
  require("data.table")
})

## Define function to read usearch mapping files in blast6out format, curate the taxonomy, and export it in a tab-delimited format ##

read_clean_tax <- function(input) {
  tax <- read.delim(input,
                    sep = "\t",
                    header = FALSE,
                    quote = "\"",
                    fill = TRUE,
                    check.names = FALSE,
                    stringsAsFactors = FALSE)
  
  #order by ESV ID
  tax <- tax[,c(1,3,2)][order(as.integer(gsub("[^0-9+$]|\\..*$", "", tax[[1]]))),]
  colnames(tax) <- c("ESV", "idty", "tax")
  
  #remove database ID's (keep everything after first " ", and remove ";" from the end if any)
  tax[["tax"]] <- gsub("^[^ ]* *|;$", "", tax[["tax"]])
  
  #split tax into individual taxonomic levels
  tax <- suppressWarnings(
    tidyr::separate(tax,
                    col = "tax",
                    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                    sep = ";"))
  
  #replace NA's with empty string and remove all entries containing any of:
  #"uncultured", "unknown", "unidentified", "incertae sedis", "metagenome", "bacterium" (whole word, so fx methanobacterium will NOT be removed),  and "possible" (whole word).
  tax[,-c(1,2)] <- lapply(tax[,-c(1,2)], function(x) {
    x[grepl("uncultured|unknown|unidentified|incertae sedis|metagenome|\\bbacterium\\b|\\bpossible\\b", tolower(x))] <- ""
    x <- stringr::str_replace_all(x,
                                  c("candidatus" = "Ca",
                                    "Candidatus" = "Ca",
                                    " " = "_",
                                    #keep only letters, numbers and ".", "-", and "_"
                                    "[^[:alnum:]_\\.\\-]" = ""
                                  ))
    return(x)
  })
  
  #remove entries below identity threshold for each taxonomic rank
  tax[which(tax$idty < 98.7), "Species"] <- ""
  tax[which(tax$idty < 94.5), "Genus"] <- ""
  tax[which(tax$idty < 86.5), "Family"] <- ""
  tax[which(tax$idty < 82.0), "Order"] <- ""
  tax[which(tax$idty < 78.5), "Class"] <- ""
  tax[which(tax$idty < 75.0), "Phylum"] <- ""
  
  tax[is.na(tax)] <- ""
  invisible(tax)
}

write_tax <- function(tax, file) {
  data.table::fwrite(tax,
                     file,
                     quote = TRUE,
                     sep = ",")
}

# Define function to read and sort data from the denovo clustering in UCLUST-format  ##
read_sort_mappings <- function(path, colnames) {
  x <- data.table::fread(path,
                         sep = "\t",
                         fill = TRUE,
                         check.names = FALSE,
                         stringsAsFactors = FALSE)
  #keep only rows with H (hits) and S (singletons) in V1, keep only V9+V10
  x <- x[V1 %in% c("H", "S"),.(V9, V10)]
  #if * in V10, replace with V9
  x <- x[,V10 := ifelse(V10 == "*", V9, V10)]
  #order by ESV ID
  x <- x[order(as.integer(gsub("[^0-9+$]|\\..*$", "", V9))),]
  colnames(x) <- colnames
  invisible(x)
}


# Fix typestrain taxonomy
## Read typestrains tax, only Genus and Species
ESV_typestrain_tax <- select(read_clean_tax("./temp/tax_typestrains.txt"), ESV, idty, Genus, Species)

#remove strain information from Species names (keep only first two words)
ESV_typestrain_tax$Species[which(ESV_typestrain_tax$Species != "")] <- 
  sapply(stringr::str_split(ESV_typestrain_tax$Species[which(ESV_typestrain_tax$Species != "")], "_"), 
         function(x) {
           paste0(x[1:2], collapse = "_")
         })

#remove ESVs that hit more than one Species using data.table
  ESV_typestrain_tax <- as.data.table(ESV_typestrain_tax)

  ESV_typestrain_tax <- ESV_typestrain_tax[,hits:=uniqueN(Species),by=ESV][hits==1,][,hits:=NULL][,idty:=NULL][Species!="",]

  ESV_typestrain_tax <- unique(ESV_typestrain_tax)

#write out
write_tax(tax = ESV_typestrain_tax, file = "./output/tax_typestrains.csv")

##### SILVA
#read SILVA tax, without Species
ESV_SILVA_tax <- select(read_clean_tax(input = "./temp/tax_SILVA.txt"), -Species)

#write out
write_tax(tax = ESV_SILVA_tax, file = "./output/tax_SILVA.csv")

##### merge typestrains+SILVA taxonomy by ESV and Genus #####
ESV_slv_typestr_tax <- left_join(ESV_SILVA_tax[,-2], ESV_typestrain_tax, by = c("ESV", "Genus"))
ESV_slv_typestr_tax[is.na(ESV_slv_typestr_tax)] <- ""

#write out
write_tax(ESV_slv_typestr_tax,
          file = "./output/tax_slv_typestr.csv")

##### denovo taxonomy #####
ESV_S <- read_sort_mappings("./temp/SILVA_ESV-S.txt", c("ESV", "Species"))
S_G <- read_sort_mappings("./temp/SILVA_S-G.txt", c("Species", "Genus"))
G_F <- read_sort_mappings("./temp/SILVA_G-F.txt", c("Genus", "Family"))
F_O <- read_sort_mappings("./temp/SILVA_F-O.txt", c("Family", "Order"))
O_C <- read_sort_mappings("./temp/SILVA_O-C.txt", c("Order", "Class"))
C_P <- read_sort_mappings("./temp/SILVA_C-P.txt", c("Class", "Phylum"))

#merge each taxonomic level according to the mapping results
denovo_tax <- left_join(ESV_S, S_G, by = "Species")
denovo_tax <- left_join(denovo_tax, G_F, by = "Genus")
denovo_tax <- left_join(denovo_tax, F_O, by = "Family")
denovo_tax <- left_join(denovo_tax, O_C, by = "Order")
denovo_tax <- left_join(denovo_tax, C_P, by = "Class")

#reorder columns and remove length from ESV ID's (".xxxx")
denovo_tax <- denovo_tax[,c("ESV", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
denovo_tax[,2:7] <- lapply(denovo_tax[,2:7], gsub, pattern = "\\..*$", replacement = "")

#generate denovo names per taxonomic level based on ESV ID
#get denovo prefix string passed from shell script, use "denovo" if blank
cmdArgs <- commandArgs(trailingOnly = TRUE)
if(length(cmdArgs) > 0) {
  prefix <- as.character(cmdArgs[[1]])
} else
  prefix <- "denovo"

denovo_tax[["Species"]] <- gsub("^[^0-9]+", paste0(prefix, "_s_"), denovo_tax[["Species"]])
denovo_tax[["Genus"]] <- gsub("^[^0-9]+", paste0(prefix, "_g_"), denovo_tax[["Genus"]])
denovo_tax[["Family"]] <- gsub("^[^0-9]+", paste0(prefix, "_f_"), denovo_tax[["Family"]])
denovo_tax[["Order"]] <- gsub("^[^0-9]+", paste0(prefix, "_o_"), denovo_tax[["Order"]])
denovo_tax[["Class"]] <- gsub("^[^0-9]+", paste0(prefix, "_c_"), denovo_tax[["Class"]])
denovo_tax[["Phylum"]] <- gsub("^[^0-9]+", paste0(prefix, "_p_"), denovo_tax[["Phylum"]])

#write out
write_tax(denovo_tax,
          file = "./output/tax_denovo.csv")

#merge SILVA+typestrains+denovo taxonomy
#merge by ESV
merged_tax <- left_join(x = ESV_slv_typestr_tax,
                        y = denovo_tax,
                        by = "ESV",
                        suffix = c("", ".denovo"))

#fill out empty entries in typestrains+SILVA with denovo taxonomy
merged_tax <- data.table(merged_tax)
merged_tax[which(Species %in% c(NA, "")), Species:=Species.denovo]
merged_tax[which(Genus %in% c(NA, "")), Genus:=Genus.denovo]
merged_tax[which(Family %in% c(NA, "")), Family:=Family.denovo]
merged_tax[which(Order %in% c(NA, "")), Order:=Order.denovo]
merged_tax[which(Class %in% c(NA, "")), Class:=Class.denovo]
merged_tax[which(Phylum %in% c(NA, "")), Phylum:=Phylum.denovo]

#order by ESV ID
merged_tax <- merged_tax[order(as.integer(gsub("[^0-9+$]|\\..*$", "", ESV))), 1:8]

##### search and replace according to a replacement file #####
replacement_file <- list.files(".",
                               pattern = "replacement",
                               ignore.case = TRUE, 
                               recursive = FALSE,
                               include.dirs = FALSE,
                               full.names = TRUE)
if(length(replacement_file) == 1 && file.exists(replacement_file)) {
  replacements <- fread(replacement_file, data.table = FALSE)
  if(any(duplicated(replacements[[1]])))
    stop("All taxa in the replacements file must be unique", call. = FALSE)
  
  #per taxonomic level, search by first column in replacements, and replace by second column
  merged_tax[] <- lapply(merged_tax, function(col) {
    unlist(lapply(col, function(taxa) {
      if(any(taxa == replacements[[1]])) {
        taxa <- replacements[which(replacements[[1]] == taxa),][[2]]
      } else
        taxa <- taxa
    }), use.names = FALSE)
  })
} else if (length(replacement_file) > 1) {
  message("More than one replacement file found, skipping...")
} else
  message("No replacement file found, skipping...")

##### fix taxa with more than one parent #####
#generate a log file first
polyTaxa <- list()
polyTaxa$Species <- merged_tax[, .(ESV = first(ESV), nParents = uniqueN(Genus), Genus = paste0(sort(unique(Genus)), collapse = ", "), new = first(Genus)), by = Species][nParents > 1,]
polyTaxa$Genus <- merged_tax[, .(ESV = first(ESV), nParents = uniqueN(Family), Family = paste0(sort(unique(Family)), collapse = ", "), new = first(Family)), by = Genus][nParents > 1,]
polyTaxa$Family <- merged_tax[, .(ESV = first(ESV), nParents = uniqueN(Order), Order = paste0(sort(unique(Order)), collapse = ", "), new = first(Order)), by = Family][nParents > 1,]
polyTaxa$Order <- merged_tax[, .(ESV = first(ESV), nParents = uniqueN(Class), Class = paste0(sort(unique(Class)), collapse = ", "), new = first(Class)), by = Order][nParents > 1,]
polyTaxa$Class <- merged_tax[, .(ESV = first(ESV), nParents = uniqueN(Phylum), Phylum = paste0(sort(unique(Phylum)), collapse = ", "), new = first(Phylum)), by = Class][nParents > 1,]
polyTaxa$Phylum <- merged_tax[, .(ESV = first(ESV), nParents = uniqueN(Kingdom), Kingdom = paste0(sort(unique(Kingdom)), collapse = ", "), new = first(Kingdom)), by = Phylum][nParents > 1,]

polyTaxaLog <- unlist(sapply(polyTaxa, function(x) {
  if(nrow(x) > 1) {
    paste0(colnames(x)[1],
           " ",
           x[[1]], 
           " has ", 
           x[[3]],
           " parents: \"",
           x[[4]],
           "\", and has been assigned the ",
           colnames(x)[4],
           " of ",
           x[[2]],
           ": ",
           x[[5]])
  }
}), use.names = FALSE)
nTaxa <- length(polyTaxaLog)

#issue a warning if one or more taxa had more than one parent, write out logfile
if(nTaxa > 0) {
  warning(paste0(nTaxa, 
                 " taxa had more than one parent, see the logfile \"./output/polyphyletics.log\" for details"), 
          call. = FALSE)
  writeLines(polyTaxaLog,
             "./output/polyphyletics.log")
}

#fix them
merged_tax[, Genus := first(Genus), by = Species]
merged_tax[, Family := first(Family), by = Genus]
merged_tax[, Order := first(Order), by = Family]
merged_tax[, Class := first(Class), by = Order]
merged_tax[, Phylum := first(Phylum), by = Class]
merged_tax[, Kingdom := first(Kingdom), by = Phylum]

#write out
write_tax(merged_tax,
          file = "./output/tax_complete.csv")

##### export as different formats #####
## export ESVs with SINTAX taxonomy in headers 
ESVs.fa <- Biostrings::readBStringSet("temp/ESVs.fa")
sintax <- left_join(data.frame(ESV = names(ESVs.fa), stringsAsFactors = FALSE), 
                    merged_tax,
                    by = "ESV")
sintax[["Kingdom"]] <- paste0("d:", sintax[["Kingdom"]])
sintax[["Phylum"]] <- paste0("p:", sintax[["Phylum"]])
sintax[["Class"]] <- paste0("c:", sintax[["Class"]])
sintax[["Order"]] <- paste0("o:", sintax[["Order"]])
sintax[["Family"]] <- paste0("f:", sintax[["Family"]])
sintax[["Genus"]] <- paste0("g:", sintax[["Genus"]])
sintax[["Species"]] <- paste0("s:", sintax[["Species"]])
sintax <- sintax[,c("ESV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
sintax_header <- paste0(sintax[["ESV"]], ";tax=", apply(sintax[,-1], 1, paste0, collapse = ","), ";")
names(ESVs.fa) <- sintax_header
Biostrings::writeXStringSet(ESVs.fa, "output/ESVs_w_sintax.fa")

## export taxonomy in QIIME format
qiime_tax <- merged_tax
qiime_tax[["Kingdom"]] <- paste0("k__", qiime_tax[["Kingdom"]])
qiime_tax[["Phylum"]] <- paste0("p__", qiime_tax[["Phylum"]])
qiime_tax[["Class"]] <- paste0("c__", qiime_tax[["Class"]])
qiime_tax[["Order"]] <- paste0("o__", qiime_tax[["Order"]])
qiime_tax[["Family"]] <- paste0("f__", qiime_tax[["Family"]])
qiime_tax[["Genus"]] <- paste0("g__", qiime_tax[["Genus"]])
qiime_tax[["Species"]] <- paste0("s__", qiime_tax[["Species"]])
qiime_tax <- qiime_tax[,c("ESV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
writeLines(paste0(qiime_tax[["ESV"]], "\t", apply(qiime_tax[,-1], 1, paste0, collapse = "; ")), "output/tax_complete_qiime.txt")
generatedenovotax
$Rscript --vanilla temp/Rscript.R $denovo_prefix
#########################################################################################
#Done, clean up
##############
cp temp/ESVs.fa output/ESVs.fa

#uncomment the below to remove temporary files
#echoWithHeader "Removing temporary files and cleaning up..."
#rm -rf temp/
duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echoWithHeader "Done in: $duration! Results are in the ./output/ folder, enjoy!"
exit 0
