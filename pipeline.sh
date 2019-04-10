#!/bin/bash
VERSION="1.2"

##### requirements (tested with) #####
# SILVA database files (SSURef+typestrains) in arb format (use latest)
# awk+grep+cat (included in most linux distributions)
# usearch (version10)
# SINA (version 1.6.0)
# R version (version 3.5.x)
#   Note: The required R packages will install automatically, but Bioconductor can cause trouble
#   if it desides to update base R packages like MASS, mgcv, lattice and others.
#   If so run R as root and run: install.packages("BiocManager"); BiocManager::install()

#################################
############# setup #############
#################################

DATA=$1
MAX_THREADS=${2:-$((`nproc`-2))}
#paths to executables
export sina=$(which sina)
export usearch=$(which usearch10)
export R=$(which R)
export Rscript=$(which Rscript)

# Paths to SILVA nr99 database, and typestrains database extracted from there. 
# .udb files have to be created first from fasta files using for example:
# $usearch10 -makeudb_usearch refdatabases/SILVA_132_SSURef_Nr99_tax_silva.fasta -output $silva_udb
silva_db="refdatabases/SILVA_132_SSURef_NR99_13_12_17_opt.arb"
silva_udb="refdatabases/SILVA_132_SSURef_Nr99_tax_silva.udb"
typestrains_db="refdatabases/SILVA132-typestrains.arb"
typestrains_udb="refdatabases/SILVA_132_SSURef_Nr99_typestrains.udb"

#de novo taxonomy prefix. Results in fx "prefix_g_123" for a de novo Genus based on ESV 123
denovo_prefix="midas"

##################################
########## end of setup ##########
##################################

#adds a header to echo, for a better console output overview
echoWithHeader() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')]: $1"
}

#check folders
if [ -d "temp" ]
    then
        echoWithHeader "A directory named 'temp' already exists and is needed for this script to run. Do you want to clear the contents and continue (y/n)?"
        read ANSWER
        if [ $ANSWER = "y" ]
            then 
            echoWithHeader "Clearing the folder 'temp'..."
            rm -rf temp/*
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
    	mkdir temp
fi

if [ -d "output" ]
	then
	  echoWithHeader "Clearing the folder 'output'..."
		rm -rf output/*
		else
		mkdir output
fi
echoWithHeader "Running fSSU taxonomy pipeline (max threads: $MAX_THREADS)..."

################################## functions ##################################
#They probably only work here in this script!
sina_align () {
  # Preparation
  local DATA=$1
  local OUTPUTID=$2
  local DB=$3
  local PT_SERVERS=${4:-1}
  local SINA_THREADS=${5:-$MAX_THREADS}
  
  $sina140 -i ${DATA} -o temp/${OUTPUTID}_aligned.fa -r $DB \
  --num-pts $PT_SERVERS --threads $SINA_THREADS \
  --log-file temp/${OUTPUTID}_log.txt

  echoWithHeader "Trimming, formatting, and sorting data..."
  #trim sequences and strip alignment gaps
  awk '!/^>/ {$0=substr($0, 1048, 41788)}1' temp/${OUTPUTID}_aligned.fa > temp/${OUTPUTID}_trimmed.fa
  $usearch10 -quiet -fasta_stripgaps temp/${OUTPUTID}_trimmed.fa -fastaout temp/tmp.fa \
    && mv temp/tmp.fa temp/${OUTPUTID}_trimmed.fa
  
  #sort sequences and stats files by ESV ID using R
  #careful with the order of arguments passed on to R
  R --slave --args temp/${OUTPUTID}_trimmed.fa << 'sortSINAoutput'
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

########################################################################################
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

    #stringr
    if(!require("stringr"))
      install.packages("stringr")

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
# note: must set threads to 1 in the below usearch10 commands or the ordering of sequences
# will be scrambled and inconsistent between individual runs. It doesn't use more than one 
# thread anyways, so no speedup will be gained.
echoWithHeader "Finding unique sequences occuring at least 2 times..."
$usearch10 -fastx_uniques $DATA -quiet -fastaout temp/preESV_wsize.fa -sizeout -minuniquesize 2 -strand both -relabel preESV -threads 1

# Orient to the same strand
echoWithHeader "Orienting sequences..."
$usearch10 -quiet -orient temp/preESV_wsize.fa -db $silva_udb -fastaout temp/preESV_wsize_oriented.fa -tabbedout temp/preESV_wsize_oriented.txt -threads 1

echoWithHeader "Clustering sequences and finding representative sequences (cluster centroids)..."
$usearch10 -quiet -cluster_fast temp/preESV_wsize_oriented.fa -sizein -sizeout -sort length -id 1 -maxrejects 0 -centroids temp/ESVs_wsize.fa -uc temp/preESVs_redundancy.uc -threads 1

#The output centroids will be ordered by size (coverage), but sequences with identical size
#will be ordered randomly between runs. The below R script first orders by size (descending) 
#and then by ESV ID (ascending) when sizes are identical. 
#Rename with new ID's to "ESV(ID).(length)" fx: "ESV1.1413"
R --slave << 'sortESVsBySizeAndID'
  #load R packages
  suppressPackageStartupMessages({
    require("Biostrings")
  })
  
  #read sequences
  ESVs <- Biostrings::readBStringSet("temp/ESVs_wsize.fa")
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
sortESVsBySizeAndID

#########################################################################################
#Align ESVs using SINA with SILVA and SILVA typestrains databases, then trim and sort
##############
#typestrains
echoWithHeader "Aligning with typestrains database using SINA..."
sina_align temp/ESVs.fa ESVs_typestrains $typestrains_db $((MAX_THREADS / 10)) $MAX_THREADS

#SILVA
echoWithHeader "Aligning with SILVA database using SINA..."
sina_align temp/ESVs.fa ESVs_SILVA $silva_db $((MAX_THREADS / 2)) $MAX_THREADS

#########################################################################################
#Assign LCA taxonomy 
##############
#typestrains
$usearch10 -quiet -usearch_global temp/ESVs_typestrains_trimmed_sorted.fa -db $typestrains_udb -maxaccepts 1 -maxrejects 0 -strand plus -id 0 -blast6out temp/tax_typestrains.txt -threads $MAX_THREADS
#SILVA
$usearch10 -quiet -usearch_global temp/ESVs_SILVA_trimmed_sorted.fa -db $silva_udb -maxaccepts 1 -maxrejects 0 -strand plus -id 0 -blast6out temp/tax_SILVA.txt -threads $MAX_THREADS

#########################################################################################
#Denovo taxonomy 
##############
#assign with identity thresholds based on Yarza et al, 2014
#using cluster_smallmem (no multithread support) and not cluster_fast to preserve order
#of input sequences, cluster_fast runs on 1 thread anyways even if set to more than 1
echoWithHeader "Clustering..."
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.987 -maxrejects 0 -uc temp/SILVA_ESV-S.txt -sortedby other
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.945 -maxrejects 0 -uc temp/SILVA_S-G.txt -sortedby other
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.865 -maxrejects 0 -uc temp/SILVA_G-F.txt -sortedby other
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.82 -maxrejects 0 -uc temp/SILVA_F-O.txt -sortedby other
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.785 -maxrejects 0 -uc temp/SILVA_O-C.txt -sortedby other
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.75 -maxrejects 0 -uc temp/SILVA_C-P.txt -sortedby other

#########################################################################################
echoWithHeader "Merging and reformatting taxonomy..."
#run R script to format, merge, and ultimately generate denovo taxonomy
$Rscript --vanilla Rscript.R $denovo_prefix

#########################################################################################
mv temp/ESVs.fa output/ESVs.fa

#uncomment the below to remove temporary files
#echoWithHeader "Removing temporary files..."
#rm -rf temp/
duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echoWithHeader "Done in: $duration! Results are in the ./output/ folder, enjoy!"
