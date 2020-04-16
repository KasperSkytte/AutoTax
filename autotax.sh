#!/bin/bash
export VERSION="1.5"

#################################
############# setup #############
#################################
# Set paths to the SILVA nr99 database and the typestrain database extracted from SILVA nr99 (see article supplementary for details). 
# .udb files have to be created first from fasta files using the example below:
# usearch11 -makeudb_usearch SILVA_132_SSURef_Nr99_tax_silva.fasta -output SILVA_132_SSURef_Nr99_tax_silva.udb
export silva_db="refdatabases/SILVA_138_SSURef_NR99_11_11_19_opt.arb"
export silva_udb="refdatabases/SILVA_138_SSURef_NR99_tax_silva.udb"
export typestrains_db="refdatabases/SILVA_138_SSURef_NR99_11_11_19_opt_typestrains.arb"
export typestrains_udb="refdatabases/SILVA_138_SSURef_NR99_tax_silva_typestrains.udb"

#de novo taxonomy prefix. Results will be in the format "prefix_g_123" for a de novo Genus based on ESV number 123
export denovo_prefix="denovo"

#set threads to a default value if not provided by the user
export MAX_THREADS=${MAX_THREADS:-$((`nproc`-2))}

##################################
########## end of setup ##########
##################################

################################## functions ##################################

#error handling
userError() {
  local self=`basename "$0"`
  echo "Invalid usage: $1" 1>&2
  echo ""
  echo "Run 'bash $self -h' for help"
}

checkRPkgs() {
  echoWithHeader "Checking for required R packages and installing if missing..."
  #Run R and check for installed packages, install if needed
  R --slave << 'checkRpkgs'
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
}

#adds a header to echo, for a better console output overview
echoWithHeader() {
  #check user arguments
  if [ ! $# -eq 1 ]
  then
    echo "Error: function must be passed exactly 1 argument" >&2
    exit 1
  fi
  echo "[$(date '+%Y-%m-%d %H:%M:%S')]: $1"
}

#check if script is run with BASH
checkBASH() {
  #check user arguments
  if [ ! $# -eq 0 ]
  then
    echo "Error: function must not be passed any arguments" >&2
    exit 1
  fi
  if [ ! -n "$BASH" ]
    then
    echo "Error: script must be run with BASH (bash)!" 1>&2
    exit 1
  fi
}

#check if a folder is present and empty
checkFolder() {
  #check user arguments
  if [ ! $# -eq 1 ]
  then
    echo "Error: function must be passed exactly 1 argument" >&2
    exit 1
  fi
  if [ -d $1 ]
      then
          echoWithHeader "A directory named '$1' already exists and is needed for this script to run. Please backup or delete the folder."
          echoWithHeader "Exiting script."
          exit 1
      else
        mkdir -p $1
  fi
}

checkInputData() {
  #exit if no data provided
  if [ -z "${DATA:-}" ]; then
    userError "No input data provided"
    exit 1
  fi
}

orient() {
  #check user arguments
  local OPTIND
  while getopts ":i:d:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      d )
        local database=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  echoWithHeader "  - Orienting sequences..."
  #note: threads must be set to 1 to make sure ordering is the same between runs
  usearch11 -orient $input -db $database -fastaout $output -threads 1 -quiet
}

derep() {
  #check user arguments
  local OPTIND
  while getopts ":i:d:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  echoWithHeader "  - Dereplicating sequences..."
  #note: threads must be set to 1 to make sure ordering is the same between runs
  usearch11 -fastx_uniques $input -fastaout $output -sizeout -minuniquesize 1 -strand plus -relabel preESV -threads 1 -quiet
}

denoise() {
  #check user arguments
  local OPTIND
  while getopts ":i:d:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  # Denoise with UNOISE3 accepting sequences seen only twice (see article supplementary for why this is acceptable)
  echoWithHeader "  - Denoising sequences using UNOISE3"
  usearch11 -unoise3 $input -zotus $output -minsize 2
  #cp temp/uniques_wsize.fa temp/preESVs.fa
}

findLongest() {
  #check user arguments
  local OPTIND
  while getopts ":i:t:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      t )
        local MAX_THREADS=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  # Remove shorter ESVs that also part of longer ESVs
  #To ensure reproducibility, the output ESVs will be ordered first by decending size (times the unique ESV has been observed) and when sizes are identical by preESV ID (ascending). 
  echoWithHeader "  - Finding the longest representative sequence of identical sequences, then reorder and rename..."
  #Rename with new ID's to "ESV(ID).(length)" fx: "ESV1.1413"
  R --slave --args "$input" "$output" "$MAX_THREADS" << 'findLongestSortESVsBySizeAndID'
    #extract passed args from shell script
    args <- commandArgs(trailingOnly = TRUE)
    input <- as.character(args[[1]])
    output <- as.character(args[[2]])
    nProc <- as.integer(args[[3]])
    
    #load R packages
    suppressPackageStartupMessages({
      require("Biostrings")
      require("doParallel")
    })
    
    #read sequences
    ESVs <- Biostrings::readBStringSet(input)
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
    Biostrings::writeXStringSet(ESVs, file = output)
findLongestSortESVsBySizeAndID
}

addESVs() {
  #check user arguments
  local OPTIND
  while getopts ":i:t:d:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      t )
        local MAX_THREADS=$OPTARG
        ;;
      d )
        local database=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  echoWithHeader "Finding new unique ESVs and adding them to the existing database..."
  cp $input output/allNewESVs.fa
  R --slave --args "$input" "$database" "$output" "$MAX_THREADS" << 'addnewESVs'
    #extract passed args from shell script
    args <- commandArgs(trailingOnly = TRUE)
    input <- args[[1]]
    ESVDB <- args[[2]]
    output <- args[[3]]
    nProc <- as.integer(args[[4]])
    
    #load R packages
    suppressPackageStartupMessages({
      require("Biostrings")
      require("doParallel")
      require("stringi")
      require("data.table")
    })
    
    doParallel::registerDoParallel(cores = nProc)
    querySeqs <- Biostrings::readBStringSet(ESVDB)
    targetSeqs <- Biostrings::readBStringSet(input)
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
  Biostrings::writeXStringSet(ESVs, output)
addnewESVs
}

#Define function to align and trim sequences based on the global SILVA alignment using SINA
sinaAlign() {
  #check user arguments
  local OPTIND
  while getopts ":i:o:d:t:l:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      o )
        local output=$OPTARG
        ;;
      d )
        local database=$OPTARG
        ;;
      t )
        local MAX_THREADS=$OPTARG
        ;;
      l )
        local logfile=$OPTARG
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
  echoWithHeader "Aligning ESVs with SILVA database using SINA..."
  sina -i $input -o $output -r $database --threads $MAX_THREADS --log-file $logfile
}

trimStripAlignment() {
  #check user arguments
  local OPTIND
  while getopts ":i:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  echoWithHeader "  - Trimming, formatting, and sorting data..."
  #trim sequences and strip alignment gaps
  awk '!/^>/ {$0=substr($0, 1048, 41788)}1' $input > $output
  usearch11 -quiet -fasta_stripgaps $output -fastaout tmp.fa \
    && mv tmp.fa $output
}

sortESVs() {
  #check user arguments
  local OPTIND
  while getopts ":i:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  #sort sequences and stats files by ESV ID using R
  R --slave --args $input $output << 'sortSINAoutput'
  #extract passed args from shell script
  args <- commandArgs(trailingOnly = TRUE)
  input <- args[[1]]
  output <- args[[2]]

  #load R packages
  suppressPackageStartupMessages({
    require("Biostrings")
  })
  
  #reorder FASTA
  algn.fa <- Biostrings::readBStringSet(input)
  algn.fa <- algn.fa[order(as.integer(gsub("[^0-9+$]|\\..*$", "", names(algn.fa))))]
  Biostrings::writeXStringSet(x = algn.fa, 
                              filepath = output)
sortSINAoutput
}

# Define functions to map ESVs against the SILVA and type strain database. For SILVA we find the top hit, 
# whereas for the type strain database we find all references with >=98.7% identity 
searchTaxDB() {
  #check user arguments
  local OPTIND
  while getopts ":i:t:d:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      t )
        local MAX_THREADS=$OPTARG
        ;;
      d )
        local database=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  echoWithHeader "Finding taxonomy of best hit in SILVA database..."
  usearch11 -usearch_global $input -db $database -maxaccepts 0 -maxrejects 0 -top_hit_only -strand plus -id 0 -blast6out $output -threads $MAX_THREADS
}

searchTaxDB_typestrain() {
  #check user arguments
  local OPTIND
  while getopts ":i:t:d:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      t )
        local MAX_THREADS=$OPTARG
        ;;
      d )
        local database=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  echoWithHeader "Finding the taxonomy of species within the 98.7% threshold in the typestrains database..."
  usearch11 -usearch_global $input -db $database -maxaccepts 0 -maxrejects 0 -strand plus -id 0.987 -blast6out $output -threads $MAX_THREADS
}

#assign with identity thresholds based on Yarza et al, 2014 using cluster_smallmem (no multithread support) to preserve order of input sequences.
clusterSpecies() {
  #check user arguments
  local OPTIND
  while getopts ":i:c:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      c )
        local centroids=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  usearch11 -quiet -cluster_smallmem $input -id 0.987 -maxrejects 0 -uc $output -centroids $centroids -sortedby other
}

clusterGenus() {
  #check user arguments
  local OPTIND
  while getopts ":i:c:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      c )
        local centroids=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  usearch11 -quiet -cluster_smallmem $input -id 0.945 -maxrejects 0 -uc $output -centroids $centroids -sortedby other
}

clusterFamily() {
  #check user arguments
  local OPTIND
  while getopts ":i:c:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      c )
        local centroids=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  usearch11 -quiet -cluster_smallmem $input -id 0.865 -maxrejects 0 -uc $output -centroids $centroids -sortedby other
}

clusterOrder() {
  #check user arguments
  local OPTIND
  while getopts ":i:c:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      c )
        local centroids=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  usearch11 -quiet -cluster_smallmem $input -id 0.82 -maxrejects 0 -uc $output -centroids $centroids -sortedby other
}

clusterClass() {
  #check user arguments
  local OPTIND
  while getopts ":i:c:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      c )
        local centroids=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  usearch11 -quiet -cluster_smallmem $input -id 0.785 -maxrejects 0 -uc $output -centroids $centroids -sortedby other
}

clusterPhylum() {
  #check user arguments
  local OPTIND
  while getopts ":i:c:o:" opt; do
    case ${opt} in
      i )
        local input=$OPTARG
        ;;
      c )
        local centroids=$OPTARG
        ;;
      o )
        local output=$OPTARG
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
  usearch11 -quiet -cluster_smallmem $input -id 0.75 -maxrejects 0 -uc $output -centroids $centroids -sortedby other
}

Rstuff() {
  cat << 'generatedenovotax' > temp/Rscript.R 
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

    ESV_typestrain_tax <- ESV_typestrain_tax[,hits:=uniqueN(Species),by=ESV][hits==1][,hits:=NULL][,idty:=NULL][Species!=""]

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
  echoWithHeader "Generating de novo taxonomy..."
  Rscript --vanilla temp/Rscript.R $denovo_prefix
cp temp/ESVs.fa output/ESVs.fa
}

echoDuration() {
  duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
  echoWithHeader "Done in: $duration! Results are in the ./output/ folder, enjoy!"
}

autotax() {
  #set appropriate error handling
  set -o errexit -o pipefail -o nounset #-o noclobber
  checkBASH
  checkInputData
  checkFolder temp
  checkFolder output
  checkRPkgs
  orient -i $DATA -d $silva_udb -o temp/fSSUs_oriented.fa
  derep -i temp/fSSUs_oriented.fa -o temp/uniques_wsize.fa
  denoise -i temp/uniques_wsize.fa -o temp/preESVs.fa
  findLongest -i temp/preESVs.fa -o temp/ESVs.fa
  #if -d is provided, identify redundant ESVs compared to the ESV database
  #and merge the two before continuing
  if [ -n "${ESVDB:-}" ]; then
    addESVs -i temp/ESVs.fa -d $ESVDB -o temp/ESVs.fa -t $MAX_THREADS
  fi
  sinaAlign -i temp/ESVs.fa -o temp/ESVs_SILVA_aln.fa -d $silva_db -t $MAX_THREADS -l temp/sinaAlign_log.txt
  trimStripAlignment -i temp/ESVs_SILVA_aln.fa -o temp/ESVs_SILVA_aln_trimmed.fa
  sortESVs -i temp/ESVs_SILVA_aln_trimmed.fa -o temp/ESVs_SILVA_aln_trimmed_sorted.fa
  searchTaxDB -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -d $silva_udb -o temp/tax_SILVA.txt -t $MAX_THREADS
  searchTaxDB_typestrain -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -d $typestrains_udb -o temp/tax_typestrains.txt -t $MAX_THREADS
  clusterSpecies -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_ESV-S.txt -c temp/SILVA_ESV-S_centroids.fa
  clusterGenus -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_S-G.txt -c temp/SILVA_S-G_centroids.fa
  clusterFamily -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_G-F.txt -c temp/SILVA_G-F_centroids.fa
  clusterOrder -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_F-O.txt -c temp/SILVA_F-O_centroids.fa
  clusterClass -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_O-C.txt -c temp/SILVA_O-C_centroids.fa
  clusterPhylum -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_C-P.txt -c temp/SILVA_C-P_centroids.fa
  Rstuff
  echoDuration
}

#only run autotax if script is not sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]
then
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
  autotax
  if [ $? -gt 0 ]
  then
    exit 1
  fi
fi
