#!/bin/bash
VERSION="1.2.5"
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
#paths to executables
export sina=$(which sina)
export usearch=$(which usearch11)
export R=$(which R)
export Rscript=$(which Rscript)

# Paths to SILVA nr99 database, and typestrains database extracted from there. 
# .udb files have to be created first from fasta files using for example:
# $usearch -makeudb_usearch refdatabases/SILVA_132_SSURef_Nr99_tax_silva.fasta -output $silva_udb
silva_db="refdatabases/SILVA_132_SSURef_NR99_13_12_17_opt.arb"
silva_udb="refdatabases/SILVA_132_SSURef_Nr99_tax_silva.udb"
typestrains_db="refdatabases/SILVA132-typestrains.arb"
typestrains_udb="refdatabases/SILVA_132_SSURef_Nr99_typestrains.udb"

#de novo taxonomy prefix. Results will be in the format "prefix_g_123" for a de novo Genus based on ESV number 123
denovo_prefix="midas"

##################################
########## end of setup ##########
##################################
# Script must be run in bash, not shell
if [ ! -n "$BASH" ]
    then
    echo "This script must be run with BASH (bash), not Shell (sh)" 1>&2
    exit 1
fi

#set appropriate error handling
set -o errexit -o pipefail -o noclobber -o nounset

################################## functions ##################################
#They probably only work here in this script!

#adds a header to echo, for a better console output overview
echoWithHeader() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')]: $1"
}

#error handling
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

#generate ESVs
generateESVs() {
  ##############
  # note: must set threads to 1 in the below usearch commands or the ordering of sequences
  # will be scrambled and inconsistent between individual runs. It doesn't use more than one 
  # thread anyways, so no speedup will be gained.
  echoWithHeader "  - Finding unique sequences occuring at least 2 times..."
  $usearch -fastx_uniques $1 -quiet -fastaout temp/preESV_wsize.fa -sizeout -minuniquesize 2 -strand both -relabel preESV -threads 1
  
  # Orient to the same strand
  echoWithHeader "  - Orienting sequences..."
  $usearch -quiet -orient temp/preESV_wsize.fa -db $silva_udb -fastaout temp/preESV_wsize_oriented.fa -tabbedout temp/preESV_wsize_oriented.txt -threads 1
  
  echoWithHeader "  - Finding the longest representative sequence of identical sequences, then reorder and rename..."
  #The output centroids will be ordered by size (coverage), but sequences with identical size
  #will be ordered randomly between runs. The below R script first orders by size (descending) 
  #and then by ESV ID (ascending) when sizes are identical. 
  #Rename with new ID's to "ESV(ID).(length)" fx: "ESV1.1413"
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
    ESVs <- Biostrings::readBStringSet("temp/preESV_wsize_oriented.fa")
    seqs_chr <- as.character(ESVs)
    
    #Find the longest of representative ESVs and remove shorter, but otherwise identical ESVs
    registerDoParallel(cores = nProc)
    removeIDs <- foreach(
      i = seq_along(seqs_chr),
      .combine = c
    ) %dopar% {
      if(any(stringi::stri_detect_fixed(str = seqs_chr[-i], pattern = seqs_chr[i])))
        return(i)
    }
    stopImplicitCluster()
    if(!is.null(removeIDs))
      ESVs <- ESVs[-removeIDs]
    
    #extract names
    headernames <- names(ESVs)
    #make a data frame with the FASTA header split in ID and size
    names <- data.frame(name = headernames,
                        ID = as.numeric(gsub("[^0-9*$]", "", gsub(";size=.*$", "", headernames))),
                        size = as.numeric(gsub(".*;size=", "", gsub(";$", "", headernames))),
                        length = lengths(ESVs),
                        check.names = FALSE,
                        stringsAsFactors = FALSE)
    #reorder the data frame first by size (descending), then by ID (ascending)
    names_ordered <- names[with(names, order(-size, ID)), ]
    #reorder the sequences
    ESVs <- ESVs[names_ordered[["name"]]]
    #rename the sequences and write out
    names(ESVs) <- paste0("ESV", 1:length(ESVs), ".", names_ordered[["length"]])
    Biostrings::writeXStringSet(ESVs, file = "temp/ESVs.fa")
findLongestSortESVsBySizeAndID
}

#align and trim sequences with SINA
sinaAlign() {
  # Preparation
  local DATA=$1
  local OUTPUTID=$2
  local DB=$3
  local SINA_THREADS=${4:-$MAX_THREADS}
  
  $sina -i ${DATA} -o temp/${OUTPUTID}_aligned.fa -r $DB \
  --threads $SINA_THREADS --log-file temp/${OUTPUTID}_log.txt

  echoWithHeader "  - Trimming, formatting, and sorting data..."
  #trim sequences and strip alignment gaps
  awk '!/^>/ {$0=substr($0, 1048, 41788)}1' temp/${OUTPUTID}_aligned.fa > temp/${OUTPUTID}_trimmed.fa
  $usearch -quiet -fasta_stripgaps temp/${OUTPUTID}_trimmed.fa -fastaout temp/tmp.fa \
    && mv temp/tmp.fa temp/${OUTPUTID}_trimmed.fa
  
  #sort sequences and stats files by ESV ID using R
  #careful with the order of arguments passed on to R
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

#search/assign tax DB
searchTaxDB() {
  IN=$1
  DB=$2
  OUT=$3
  $usearch -quiet -usearch_global $IN -db $DB -maxaccepts 1 -maxrejects 0 -strand plus -id 0 -blast6out $OUT -threads $MAX_THREADS
}

################################## end of functions ##################################
################################## start of pipeline ##################################
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

#########################################################################################
#Generate ESVs
##############
echoWithHeader "Generating ESVs..."
generateESVs $DATA

#########################################################################################
#if -d is provided, identify redundant ESVs compared to the ESV database
#and merge the two before continuing
##############
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
	  nReplacements <- nrow(replaceSeqs)
	  if(nReplacements > 0) {
	    replaceSeqs[,diff := gsub(querySeq, "...", targetSeq), by = queryName]
	    replaceIDs <- which(names(querySeqs) %in% replaceSeqs[["queryName"]])
	    querySeqs[replaceIDs] <- replaceSeqs[["targetSeq"]]
	    #adjust length annotation in the names of the new ESVs
	    names(querySeqs)[replaceIDs] <- paste0(gsub("\\..*$", ".", names(querySeqs)[replaceIDs]), 
	                                           Biostrings::width(querySeqs[replaceIDs]))
	    
	    #write information about the replacements to a log file
	    replacementLog <- paste0(replaceSeqs[,queryName], 
	                             " has been replaced with ",
	                             replaceSeqs[,targetName], 
	                             ", difference: ",
	                             replaceSeqs[,diff])
	    warning(paste0(nReplacements, 
	                   if(nReplacements > 1) 
	                     " ESVs have" 
	                   else 
	                     " ESV has",
	                   " been replaced by a longer ESV, see the logfile \"./output/replacedESVs.log\" for details"))
	    writeLines(replacementLog,
	               "./output/replacedESVs.log")
	  }
	}

	#find new unique sequences that are not in the existing database
	#and add to the database with new names continuing ID numbering
	lastID <- as.integer(gsub("[^0-9+$]|\\..*$", "", names(querySeqs[length(querySeqs)])))
	names(newESVs) <- paste0("ESV", 1:length(newESVs)+lastID, ".", gsub("^.*\\.", "", names(newESVs)))
	ESVs <- c(querySeqs, newESVs)
	Biostrings::writeXStringSet(ESVs, "temp/ESVs.fa")
addnewESVs
fi

#########################################################################################
#Align ESVs using SINA with SILVA and SILVA typestrains databases, then trim and sort
##############
#typestrains
echoWithHeader "Aligning ESVs with typestrains database using SINA..."
sinaAlign temp/ESVs.fa ESVs_typestrains $typestrains_db $MAX_THREADS

#SILVA
echoWithHeader "Aligning ESVs with SILVA database using SINA..."
sinaAlign temp/ESVs.fa ESVs_SILVA $silva_db $MAX_THREADS

#########################################################################################
#Assign taxonomy of best hit
##############
#typestrains
echoWithHeader "Finding taxonomy of best hit in typestrains database..."
searchTaxDB temp/ESVs_typestrains_trimmed_sorted.fa $typestrains_udb temp/tax_typestrains.txt

#SILVA
echoWithHeader "Finding taxonomy of best hit in SILVA database..."
searchTaxDB temp/ESVs_SILVA_trimmed_sorted.fa $silva_udb temp/tax_SILVA.txt

#########################################################################################
#De novo taxonomy 
##############
#assign with identity thresholds based on Yarza et al, 2014
#using cluster_smallmem (no multithread support) and not cluster_fast to preserve order
#of input sequences, cluster_fast runs on 1 thread anyways even if set to more than 1
echoWithHeader "Generating de novo taxonomy..."
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.987 -maxrejects 0 -uc temp/SILVA_ESV-S.txt -sortedby other
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.945 -maxrejects 0 -uc temp/SILVA_S-G.txt -sortedby other
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.865 -maxrejects 0 -uc temp/SILVA_G-F.txt -sortedby other
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.82 -maxrejects 0 -uc temp/SILVA_F-O.txt -sortedby other
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.785 -maxrejects 0 -uc temp/SILVA_O-C.txt -sortedby other
$usearch -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.75 -maxrejects 0 -uc temp/SILVA_C-P.txt -sortedby other

#########################################################################################
echoWithHeader "Merging and reformatting taxonomy..."
#run R script to format, merge, and ultimately generate de novo taxonomy
cat << 'generatedenovotax' > temp/Rscript.R
#!/usr/bin/Rscript
#load R packages
suppressPackageStartupMessages({
  require("stringr")
  require("tidyr")
  require("dplyr")
  require("data.table")
})

##### FUNCTIONS #####
## read/write taxonomy ##
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
  #"uncultured"
  #"unknown"
  #"unidentified"
  #"incertae sedis"
  #"metagenome"
  #"bacterium" (whole word, so fx methanobacterium will NOT be removed)
  #"possible" (whole word)
  tax[,-c(1,2)] <- lapply(tax[,-c(1,2)], function(x) {
    x[grepl("uncultured|unknown|unidentified|incertae sedis|metagenome|\\bbacterium\\b|\\bpossible\\b", tolower(x))] <- ""
    x <- stringr::str_replace_all(x,
                                  c("candidatus" = "Candidatus",
                                    " " = "_",
                                    #keep only letters, numbers and ".", "-", and "_"
                                    "[^[:alnum:]_\\.\\-]" = ""
                                  ))
    return(x)
  })
  
  #remove entries below identity threshold per level
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

## read and sort mappings ##
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
  x <- x[order(as.integer(gsub("[^0-9+$]", "", V9))),]
  colnames(x) <- colnames
  invisible(x)
}

##### Fix taxonomy #####
##### typestrains
#read typestrains tax, only Genus and Species
ESV_typestrain_tax <- select(read_clean_tax("temp/tax_typestrains.txt"), ESV, idty, Genus, Species)

#remove strain information from Species names (keep only first two words)
ESV_typestrain_tax$Species[which(ESV_typestrain_tax$Species != "")] <- 
  sapply(stringr::str_split(ESV_typestrain_tax$Species[which(ESV_typestrain_tax$Species != "")], "_"), 
         function(x) {
           paste0(x[1:2], collapse = "_")
         })

#write out
write_tax(tax = ESV_typestrain_tax,
          file = "./output/tax_typestrains.csv")

##### SILVA
#read SILVA tax, without Species
ESV_SILVA_tax <- select(read_clean_tax(input = "./temp/tax_SILVA.txt"), -Species)

#write out
write_tax(tax = ESV_SILVA_tax,
          file = "./output/tax_SILVA.csv")

##### merge typestrains+SILVA taxonomy by ESV and Genus #####
ESV_slv_typestr_tax <- left_join(ESV_SILVA_tax[,-2], ESV_typestrain_tax[,-2], by = c("ESV", "Genus"))
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
#make data.tables
ESV_slv_typestr_tax <- data.table(ESV_slv_typestr_tax, key = "ESV")
denovo_tax <- data.table(denovo_tax, key = "ESV")

#merge by ESV
merged_tax <- ESV_slv_typestr_tax[denovo_tax]

#fill out empty entries in typestrains+SILVA with denovo taxonomy
merged_tax[which(Species %in% c(NA, "")), Species:=i.Species]
merged_tax[which(Genus %in% c(NA, "")), Genus:=i.Genus]
merged_tax[which(Family %in% c(NA, "")), Family:=i.Family]
merged_tax[which(Order %in% c(NA, "")), Order:=i.Order]
merged_tax[which(Class %in% c(NA, "")), Class:=i.Class]
merged_tax[which(Phylum %in% c(NA, "")), Phylum:=i.Phylum]

#order by ESV ID
merged_tax <- merged_tax[order(as.integer(gsub("[^0-9+$]", "", ESV))), 1:8]

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

#issue a warning if one or more taxa have more than one parent, write out logfile
if(nTaxa > 0) {
  warning(paste0(nTaxa, 
                 " taxa ",
                 if(nTaxa > 1)
                   "have" 
                 else 
                   "has",
                 " more than one parent, see the logfile \"./output/polyphyletics.log\" for details"), 
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
mv temp/ESVs.fa output/ESVs.fa

#uncomment the below to remove temporary files
#echoWithHeader "Removing temporary files and cleaning up..."
#rm -rf temp/
duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echoWithHeader "Done in: $duration! Results are in the ./output/ folder, enjoy!"
exit 0