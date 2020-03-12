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
# Set paths to the SILVA nr99 database and the typestrain database extracted from SILVA nr99 (see article supplementary for details). 
# .udb files have to be created first from fasta files using the example below:
# usearch11 -makeudb_usearch SILVA_132_SSURef_Nr99_tax_silva.fasta -output SILVA_132_SSURef_Nr99_tax_silva.udb
silva_db="refdatabases/SILVA_138_SSURef_NR99_11_11_19_opt.arb"
silva_udb="refdatabases/SILVA_138_SSURef_NR99_tax_silva.udb"
typestrains_db="refdatabases/SILVA_138_SSURef_NR99_11_11_19_opt_typestrains.arb"
typestrains_udb="refdatabases/SILVA_138_SSURef_NR99_tax_silva_typestrains.udb"

#de novo taxonomy prefix. Results will be in the format "prefix_g_123" for a de novo Genus based on ESV number 123
denovo_prefix="denovo"

#set threads to a default value if not provided by the user
MAX_THREADS=${MAX_THREADS:-$((`nproc`-2))}

##################################
########## end of setup ##########
##################################

################################## functions ##################################

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

# Define function to generate ESVs
generateESVs() {
  ##############
  # note: must set threads to 1 in the below usearch commands or the ordering of sequences
  # will be scrambled and inconsistent between individual runs. It doesn't use more than one 
  # thread anyways, so no speedup will be gained.
  
  # Orient all fSSU correctly based on the SILVA database
  echoWithHeader "  - Orienting sequences..."
  usearch11 -orient $1 -db $silva_udb -fastaout temp/fSSUs_oriented.fa -threads 1 -quiet
  
  # Dereplicate sequences
  echoWithHeader "  - Dereplicating sequences..."
  usearch11 -fastx_uniques temp/fSSUs_oriented.fa -fastaout temp/uniques_wsize.fa -sizeout -minuniquesize 1 -strand plus -relabel preESV -threads 1 -quiet
  
  # Denoise with UNOISE3 accepting sequences seen only twice (see article supplementary for why this is acceptable)
  echoWithHeader "  - Denoising sequences using UNOISE3"
  usearch11 -unoise3 temp/uniques_wsize.fa -zotus temp/preESVs.fa -minsize 2
  #cp temp/uniques_wsize.fa temp/preESVs.fa

  # Remove shorter ESVs that also part of longer ESVs
  #To ensure reproducibility, the output ESVs will be ordered first by decending size (times the unique ESV has been observed) and when sizes are identical by preESV ID (ascending). 
  echoWithHeader "  - Finding the longest representative sequence of identical sequences, then reorder and rename..."
  #Rename with new ID's to "ESV(ID).(length)" fx: "ESV1.1413"
  R --slave --args "$MAX_THREADS" << 'findLongestSortESVsBySizeAndID'
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

#Define function to align and trim sequences based on the global SILVA alignment using SINA

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