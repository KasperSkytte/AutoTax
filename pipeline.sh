#!/bin/bash
# usage: ./pipeline.sh fSSU.fa 36

##### requirements (tested with) #####
# SILVA database files (SSURef+typestrains) in arb format (use latest)
# usearch (version10)
# SINA (version 1.4.0)
# R version (version 3.5.1)
#   Note: The required R packages will install automatically, but Bioconductor can cause trouble
#   if it desides to update base R packages like MASS, mgcv, lattice and others.
#   If so run R as root and run: install.packages("BiocManager"); BiocManager::install()

##### setup #####
DATA=$1
MAX_THREADS=${2:-36}
#path to sina executable
export sina=/home/kapper/.conda/envs/sina/bin/sina

#path to usearch executable
#export usearch10=/opt/bin/usearch10

# Provide the path the the SILVA nr99 arb-database and the typestrain database extracted from the SILVA nr99 arb-database.
silva_db="refdatabases/SILVA_132_SSURef_NR99_13_12_17_opt.arb"
typestrain_db="refdatabases/SILVA132-typestrains.arb"

##### end of setup #####
#adds a header to echo, for a better console output overview
echoWithHeader() {
  header() {
    echo "################################################################################"
  }
  header
  echo $1
  header
}
date

echoWithHeader "Running fSSU taxonomy pipeline..."
#check folders
if [ -d "temp" ]
    then
        echo ""
        echo "A directory named 'temp' already exists and is needed for this script to run. Do you want to clear the contents and continue (y/n)?"
        read ANSWER
        if [ $ANSWER = "y" ]
            then 
            echo "Clearing the folder 'temp'..."
            rm -rf temp/*
            else 
            if [ $ANSWER = "n" ]
                then
                echo "    Exiting script."
                echo ""
                exit 0
                else
                echo "    Exiting script."
                echo ""
                exit 1
            fi
        fi
    else
    	mkdir temp
fi

if [ -d "output" ]
	then
	  echo "Clearing the folder 'output'..."
		rm -rf output/*
		else
		mkdir output
fi

################################## functions ##################################
#They probably only work here in this script!
#Run SINA using GNU parallel
sina_align_par () {
  # Preparation
  local DATA=$1
  local OUTPUTID=$2
  local DB=$3
  local LCA_FIELDS=$4
  local JOBS=${5:-MAX_THREADS}

  echo "Running SINA in parallel using GNU parallel..."
  #determine optimal block size for each job
  FILE_SIZE=$(stat --printf="%s" ${DATA})
  BLOCK_SIZE=$((FILE_SIZE / (JOBS+1) / 1000))
  if [ ${BLOCK_SIZE} -lt 10 ] ; then
    BLOCK_SIZE=10
  fi
  if [ ${BLOCK_SIZE} -gt 2000000 ] ; then
    BLOCK_SIZE=2000000
  fi
  BLOCK_SIZE=${BLOCK_SIZE}k
  echo "Splitting data into blocks with size: $BLOCK_SIZE"
  
  #duplicate fasta headers (cause SINA is weird) and
  #split FASTA file in smaller files and run SINA in parallel
  if [ -d "parallel_temp" ]
    then
      echo "Clearing the folder 'parallel_temp'..."
      rm -rf parallel_temp/*
    else
    mkdir parallel_temp
  fi
  
  cat ${DATA} | awk -F' ' -vOFS=' ' '/^>/ {$1=$1 " " substr($1,2,30)} {print}' |\
    parallel --progress -j $JOBS --block $BLOCK_SIZE --recstart '>' \
    --pipe "cat > parallel_temp/temp_block_{#}.fa; $sina -i parallel_temp/temp_block_{#}.fa \
    -o parallel_temp/${OUTPUTID}_block_{#}_aligned.fa --meta-fmt csv -r $DB --search \
    --calc-idty --lca-fields $LCA_FIELDS --search-min-sim 0.5 --search-max-result 1 \
    --num-pts 1 --threads 1"

  echo ""
  echo "Concatenating, sorting, reformatting, and trimming data..."
  
  #combine results
  cat parallel_temp/${OUTPUTID}_block_*.fa > temp/${OUTPUTID}_aligned.fa
  head -n 1 parallel_temp/${OUTPUTID}_block_1_aligned.fa.csv > temp/${OUTPUTID}_aligned.csv
  find parallel_temp/ -name "${OUTPUTID}_block_*.fa.csv" | xargs -n 1 tail -n +2 >> temp/${OUTPUTID}_aligned.csv
  
  #trim sequences
  awk '!/^>/ {$0=substr($0, 1048, 41788)}1' temp/${OUTPUTID}_aligned.fa | \
    sed 's/ ESV[0-9]*//g' | tr -d "-" | tr -d "." > temp/${OUTPUTID}_trimmed.fa
  
  #remove duplicate names, and any []() from taxonomy
  sed -i -e 's/[][)(]//g' -e 's/ ESV[0-9]*//g' temp/${OUTPUTID}_aligned.csv

  #sort sequences and stats files by ESV ID using R
  #careful with the order of arguments passed on to R
  R --slave --args temp/${OUTPUTID}_trimmed.fa temp/${OUTPUTID}_aligned.csv ${LCA_FIELDS} << 'sortSINAoutput'
  #extract passed args from shell script
  args <- commandArgs(trailingOnly = TRUE)

  #load R packages
  suppressPackageStartupMessages({
    require("Biostrings")
  })
  
  #reorder FASTA
  algn.fa <- Biostrings::readBStringSet(args[[1]])
  algn.fa <- algn.fa[order(as.integer(gsub("[^0-9+$]", "", names(algn.fa))))]
  Biostrings::writeXStringSet(x = algn.fa, 
                              filepath = paste0(tools::file_path_sans_ext(args[[1]]),
                                                "_sorted.", 
                                                tools::file_ext(args[[1]])))
  
  #reorder stats
  #use base read.csv, data.table::fread() doesn't understand separators within quotes
  stats <- read.csv(args[[2]],
                    sep = ",",
                    header = TRUE,
                    quote = "\"",
                    fill = TRUE,
                    check.names = FALSE,
                    stringsAsFactors = FALSE)
  stats <- stats[order(as.integer(gsub("[^0-9+$]", "", stats[["full_name"]]))),c("full_name", "align_ident_slv", paste0("lca_", args[[3]]))]
  write.table(stats,
              paste0(tools::file_path_sans_ext(args[[2]]), 
                           "_sorted.", 
                           tools::file_ext(args[[2]])),
              sep = ",",
              na = "",
              append = FALSE,
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
sortSINAoutput

  echo "Removing temporary files and killing all running ARB pt servers running in the background..."
  rm -rf parallel_temp
  pgrep -u $(whoami) arb_pt_server
  pkill -u $(whoami) arb_pt_server
  #kill -9 $(ps -auxww | grep '[p]t_server' | awk '{print $2}')
}

sina_140_align () {
  # Preparation
  local DATA=$1
  local OUTPUTID=$2
  local DB=$3
  local LCA_FIELDS=$4
  local PT_SERVERS=${5:-1}
  local SINA_THREADS=${6:-$MAX_THREADS}
  
  cat ${DATA} | awk -F' ' -vOFS=' ' '/^>/ {$1=$1 " " substr($1,2,30)} {print}' |\
    $sina -o temp/${OUTPUTID}_aligned.fa --meta-fmt csv -r $DB --search \
    --calc-idty --lca-fields $LCA_FIELDS --search-min-sim 0.5 \
    --search-max-result 1 --num-pts $PT_SERVERS --threads $SINA_THREADS
  
  echo "Trimming, formatting, and sorting data"
  #trim sequences
  awk '!/^>/ {$0=substr($0, 1048, 41788)}1' temp/${OUTPUTID}_aligned.fa | \
    sed 's/ ESV[0-9]*//g' | tr -d "-" | tr -d "." > temp/${OUTPUTID}_trimmed.fa
  #remove duplicate names, and any []() from taxonomy
  sed -i -e 's/[][)(]//g' -e 's/ ESV[0-9]*//g' temp/${OUTPUTID}_aligned.csv
  
  #sort sequences and stats files by ESV ID using R
  #careful with the order of arguments passed on to R
  R --slave --args temp/${OUTPUTID}_trimmed.fa temp/${OUTPUTID}_aligned.csv ${LCA_FIELDS} << 'sortSINAoutput'
	#extract passed args from shell script
	args <- commandArgs(trailingOnly = TRUE)

  #load R packages
  suppressPackageStartupMessages({
    require("Biostrings")
  })
  
  #reorder FASTA
  algn.fa <- Biostrings::readBStringSet(args[[1]])
  algn.fa <- algn.fa[order(as.integer(gsub("[^0-9+$]", "", names(algn.fa))))]
  Biostrings::writeXStringSet(x = algn.fa, 
                              filepath = paste0(tools::file_path_sans_ext(args[[1]]),
                                                "_sorted.", 
                                                tools::file_ext(args[[1]])))
  
  #reorder stats
  #use base read.csv, data.table::fread() doesn't understand separators within quotes
  stats <- read.csv(args[[2]],
                    sep = ",",
                    header = TRUE,
                    quote = "\"",
                    fill = TRUE,
                    check.names = FALSE,
                    stringsAsFactors = FALSE)
  stats <- stats[order(as.integer(gsub("[^0-9+$]", "", stats[["full_name"]]))),c("full_name", "align_ident_slv", paste0("lca_", args[[3]]))]
  write.table(stats,
              paste0(tools::file_path_sans_ext(args[[2]]), 
                           "_sorted.", 
                           tools::file_ext(args[[2]])),
              sep = ",",
              na = "",
              append = FALSE,
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
sortSINAoutput
}

###################

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
# note: must set threads to 1 in the below usearch10 commands or the ordering of sequences
# will be scrambled and inconsistent between individual runs. It doesn't use more than one 
# thread anyways, so no speedup will be gained.
echoWithHeader "Finding unique sequences occuring at least 2 times..."
usearch10 -fastx_uniques $DATA -fastaout temp/preESV_wsize.fa -sizeout -minuniquesize 2 -strand both -relabel preESV -threads 1

### Orient to the same strand
# Create reference database, only needed once
#usearch10 -makeudb_usearch refdatabases/SILVA_132_SSURef_Nr99_tax_silva.fasta -output refdatabases/SILVA_132_SSURef_Nr99_tax_silva.udb

# Run the orientation check
echoWithHeader "Orienting sequences..."
usearch10 -orient temp/preESV_wsize.fa -db refdatabases/SILVA_132_SSURef_Nr99_tax_silva.udb -fastaout temp/preESV_wsize_oriented.fa -tabbedout temp/preESV_wsize_oriented.txt -threads 1

echoWithHeader "Clustering sequences and finding representative sequences (cluster centroids)..."
usearch10 -cluster_fast temp/preESV_wsize_oriented.fa -sizein -sizeout -sort length -id 1 -maxrejects 0 -centroids temp/ESVs_wsize.fa -uc temp/preESVs_redundancy.uc -threads 1

#The output centroids will be ordered by size (coverage), but sequences with identical size
#will be ordered randomly between runs. The below R script first orders by size (descending) 
#and then by ESV ID (ascending) when sizes are identical
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
                      ID = as.numeric(gsub("[^0-9*$]","", gsub(";size=.*$", "", headernames))),
                      size = as.numeric(gsub(".*;size=", "", gsub(";$", "", headernames))),
                      check.names = FALSE,
                      stringsAsFactors = FALSE)
  #reorder the data frame first by size (descending), then by ID (ascending)
  names_ordered <- names[with(names, order(-size, ID)), ]
  #reorder the sequences
  ESVs <- ESVs[names_ordered[["name"]]]
  #rename the sequences and write out
  names(ESVs) <- paste0("ESV", 1L:length(ESVs))
  Biostrings::writeXStringSet(ESVs, file = "temp/ESVs.fa")
sortESVsBySizeAndID

#########################################################################################
## DONT USE ##
#to add new sequences in the future use something like this
# Example command: ./PrepareESV.sh all_new_fSSUs.fa oldESV.fa
#usearch10 -search_global new_ESV.fa -db $2 -id 1 -strand plus -maxrejects 0 -notmatched new_unique_ESVs.fa -threads $MAX_THREADS
#cat $2 new_unique_ESVs.fa > updated_ESVs_temp.fa
#usearch10 -fastx_relabel updated_ESVs_temp.fa -prefix ESV -fastaout ./ESV/updated_ESVs.fa -threads $MAX_THREADS
#rm new_unique_ESVs.fa
#rm updated_ESVs_temp.fa
#rm ESVs.fa
## DONT USE ##

#########################################################################################
#typestrains
echoWithHeader "Aligning with typestrains database and assigning taxonomy using SINA..."
#with GNU parallel
sina_align_par temp/ESVs.fa ESVs_typestrains $typestrain_db tax_ltp_name $MAX_THREADS
#with SINA 1.4.0 with multithreading
#sina_140_align temp/ESVs.fa ESVs_typestrains $typestrain_db tax_ltp_name $((MAX_THREADS / 10)) $MAX_THREADS

#########################################################################################
#SILVA
#careful with too many jobs here, the SILVA DB causes SINA to use about
#5-6GB of working memory per job, so set accordingly! 16 jobs requires ~100GB of memory
echoWithHeader "Aligning with SILVA database and assigning taxonomy using SINA..."
#with GNU parallel
sina_align_par temp/ESVs.fa ESVs_SILVA $silva_db tax_slv 16
#sina_140_align temp/ESVs.fa ESVs_SILVA $silva_db tax_slv $((MAX_THREADS / 5)) $MAX_THREADS

#assign with identity thresholds based on Yarza et al, 2014
#using cluster_smallmem (no multithread support) and not cluster_fast to preserve order
#of input sequences, and cluster_fast runs on 1 thread anyways even if set to more than 1
echoWithHeader "Clustering..."
usearch10 -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.987 -maxrejects 0 -centroids temp/SILVA_DeNovoSpecies.fa -uc temp/SILVA_ESV-S.txt -sortedby other
usearch10 -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.945 -maxrejects 0 -centroids temp/SILVA_DeNovoGenus.fa -sortedby other
usearch10 -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.865 -maxrejects 0 -centroids temp/SILVA_DeNovoFamily.fa -sortedby other
usearch10 -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.82 -maxrejects 0 -centroids temp/SILVA_DeNovoOrder.fa -sortedby other
usearch10 -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.785 -maxrejects 0 -centroids temp/SILVA_DeNovoClass.fa -sortedby other
usearch10 -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.75 -maxrejects 0 -centroids temp/SILVA_DeNovoPhylum.fa -sortedby other

#uncomment the below to sort by ESV ID
# R --slave << EOF
# library(Biostrings)
# files <- c("temp/SILVA_DeNovoSpecies.fa",
#            "temp/SILVA_DeNovoGenus.fa",
#            "temp/SILVA_DeNovoFamily.fa",
#            "temp/SILVA_DeNovoOrder.fa",
#            "temp/SILVA_DeNovoClass.fa",
#            "temp/SILVA_DeNovoPhylum.fa")
# lapply(files, function(x) {
#   s <- Biostrings::readBStringSet(x)
#   s <- s[order(as.numeric(unlist(stringr::str_extract_all(names(s), "[:digit:]+$"), use.names = FALSE)))]
#   Biostrings::writeXStringSet(x = s, x)
# })
# EOF

#map per taxonomic level
echoWithHeader "Mapping each taxonomic level..."
usearch10 -usearch_global temp/SILVA_DeNovoSpecies.fa -db temp/SILVA_DeNovoGenus.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_S-G.txt -threads $MAX_THREADS
usearch10 -usearch_global temp/SILVA_DeNovoGenus.fa -db temp/SILVA_DeNovoFamily.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_G-F.txt -threads $MAX_THREADS
usearch10 -usearch_global temp/SILVA_DeNovoFamily.fa -db temp/SILVA_DeNovoOrder.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_F-O.txt -threads $MAX_THREADS
usearch10 -usearch_global temp/SILVA_DeNovoOrder.fa -db temp/SILVA_DeNovoClass.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_O-C.txt -threads $MAX_THREADS
usearch10 -usearch_global temp/SILVA_DeNovoClass.fa -db temp/SILVA_DeNovoPhylum.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_C-P.txt -threads $MAX_THREADS

#########################################################################################
echoWithHeader "Merging and reformatting taxonomy..."
#Use R to fix output
R --slave << 'generateTaxonomy'
  #load R packages
  suppressPackageStartupMessages({
    require("stringr")
    require("tidyr")
    require("dplyr")
    require("data.table")
  })
  
  ##### FUNCTIONS #####
  ## read/write taxonomy ##
  read_tax <- function(input, tax_field) {
    tax <- read.csv(input,
                    sep = ",",
                    header = TRUE,
                    quote = "\"",
                    fill = TRUE,
                    check.names = FALSE,
                    stringsAsFactors = FALSE)
    tax <- tax[,c("full_name", "align_ident_slv", tax_field)]
    colnames(tax)[3] <- "tax"
    tax[["tax"]] <- stringr::str_replace_all(tax[["tax"]], 
                                             c("candidatus " = "Candidatus_",
                                               "Candidatus " = "Candidatus_", 
                                               ";$" = ""))
    invisible(tax)
  }
  
  write_tax <- function(tax, file) {
    fwrite(tax,
           file,
           sep = ",")
  }
  
  ## read and sort mappings ##
  read_sort_mappings <- function(input) {
    x <- fread(input,
               sep = "\t",
               fill = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE)
    x <- x[,1:2]
    colnames(x) <- c("V1", "V2")
    x <- x[order(as.integer(gsub("[^0-9+$]", "", V1))),]
    # #write out
    # fwrite(S_G,
    #        file = "./output/SILVA_S-G.txt",
    #        sep = "\t",
    #        row.names = FALSE,
    #        col.names = FALSE)
    invisible(x)
  }
  
  ##### Fix taxonomy #####
  #typestrains
  #read typestrains stats
  ESV_typestrain_algn <- read_tax(input = "./temp/ESVs_typestrains_aligned.csv",
                                  tax_field = "lca_tax_ltp_name")
  
  #split tax into individual taxonomy variables
  ESV_typestrain_algn <- suppressWarnings(
    separate(ESV_typestrain_algn,
             col = "tax",
             into = c("Genus", "Species"),
             sep = " "))
  
  #remove entries below identity threshold per level
  ESV_typestrain_algn[which(ESV_typestrain_algn$align_ident_slv < 98.7), "Species"] <- ""
  ESV_typestrain_algn[which(ESV_typestrain_algn$align_ident_slv < 94.5), "Genus"] <- ""
  
  #replace NA's with empty string
  ESV_typestrain_algn[is.na(ESV_typestrain_algn)] <- ""
  
  #remove all entries containing "uncultured", "unknown", and "incertae sedis"
  ESV_typestrain_algn[,-c(1,2)] <- lapply(ESV_typestrain_algn[,-c(1,2)], function(x) {
    x[grepl("uncultured|unknown|incertae sedis", tolower(x))] <- ""
    return(x)
  })
  
  #write out
  write_tax(tax = ESV_typestrain_algn,
            file = "./output/tax_typestrains.csv")
  
  #SILVA
  #read SILVA stats
  ESV_SILVA_algn <- read_tax(input = "./temp/ESVs_SILVA_aligned.csv",
                             tax_field = "lca_tax_slv")
  
  #split tax into individual taxonomy variables
  ESV_SILVA_algn <- suppressWarnings(
    separate(ESV_SILVA_algn,
             col = "tax",
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
             sep = ";"))
  
  #remove entries below identity threshold per level
  ESV_SILVA_algn[which(ESV_SILVA_algn$align_ident_slv < 94.5), "Genus"] <- ""
  ESV_SILVA_algn[which(ESV_SILVA_algn$align_ident_slv < 86.5), "Family"] <- ""
  ESV_SILVA_algn[which(ESV_SILVA_algn$align_ident_slv < 82.0), "Order"] <- ""
  ESV_SILVA_algn[which(ESV_SILVA_algn$align_ident_slv < 78.5), "Class"] <- ""
  ESV_SILVA_algn[which(ESV_SILVA_algn$align_ident_slv < 75.0), "Phylum"] <- ""
  
  #replace NA's with empty string
  ESV_SILVA_algn[is.na(ESV_SILVA_algn)] <- ""
  
  #remove all entries containing "uncultured", "unknown", and "incertae sedis"
  ESV_SILVA_algn[,-c(1,2)] <- lapply(ESV_SILVA_algn[,-c(1,2)], function(x) {
    x[grepl("uncultured|unknown|incertae sedis", tolower(x))] <- ""
    return(x)
  })
  
  #write out
  write_tax(tax = ESV_SILVA_algn,
            file = "./output/tax_SILVA.csv")
  
  #merge typestrains+SILVA taxonomy by ESV and Genus
  taxonomy <- left_join(ESV_SILVA_algn[,-2], ESV_typestrain_algn[,-2], by = c("full_name", "Genus"))
  colnames(taxonomy)[1] <- "ESV"
  taxonomy[is.na(taxonomy)] <- ""
  write_tax(taxonomy,
            file = "./output/tax_slv_typestr.csv")
  
  ##### Fix mappings + denovo taxonomy #####
  ESV_S <- fread("./temp/SILVA_ESV-S.txt",
                 sep = "\t",
                 fill = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE)
  #keep only rows with H (hits) and S (singletons) in V1, keep only V9+V10
  ESV_S <- ESV_S[V1 %in% c("H", "S"),.(V9, V10)]
  #if * in V10, replace with V9
  ESV_S <- ESV_S[,V10 := ifelse(V10 == "*", V9, V10)]
  colnames(ESV_S) <- c("ESV", "Species")
  
  #order by ESV ID
  ESV_S <- ESV_S[order(as.integer(gsub("[^0-9+$]", "", ESV))),]
  
  #sort
  S_G <- read_sort_mappings("temp/SILVA_S-G.txt")
  colnames(S_G) <- c("Species", "Genus")
  G_F <- read_sort_mappings("temp/SILVA_G-F.txt")
  colnames(G_F) <- c("Genus", "Family")
  F_O <- read_sort_mappings("temp/SILVA_F-O.txt")
  colnames(F_O) <- c("Family", "Order")
  O_C <- read_sort_mappings("temp/SILVA_O-C.txt")
  colnames(O_C) <- c("Order", "Class")
  C_P <- read_sort_mappings("temp/SILVA_C-P.txt")
  colnames(C_P) <- c("Class", "Phylum")
  
  #merge each taxonomic level mapping
  denovo_midas <- left_join(ESV_S, S_G, by = "Species")
  denovo_midas <- left_join(denovo_midas, G_F, by = "Genus")
  denovo_midas <- left_join(denovo_midas, F_O, by = "Family")
  denovo_midas <- left_join(denovo_midas, O_C, by = "Order")
  denovo_midas <- left_join(denovo_midas, C_P, by = "Class")
  denovo_midas <- denovo_midas[,c("ESV", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
  
  #generate new names per taxonomic level
  denovo_midas[["Species"]] <- gsub("ESV", "midas_s_", denovo_midas[["Species"]])
  denovo_midas[["Genus"]] <- gsub("ESV", "midas_g_", denovo_midas[["Genus"]])
  denovo_midas[["Family"]] <- gsub("ESV", "midas_f_", denovo_midas[["Family"]])
  denovo_midas[["Order"]] <- gsub("ESV", "midas_o_", denovo_midas[["Order"]])
  denovo_midas[["Class"]] <- gsub("ESV", "midas_c_", denovo_midas[["Class"]])
  denovo_midas[["Phylum"]] <- gsub("ESV", "midas_p_", denovo_midas[["Phylum"]])
  
  #write out
  write_tax(denovo_midas,
            file = "./output/tax_midas_denovo.csv")
  
  #merge SILVA+typestrains+denovo MiDAS taxonomy
  #make data.tables
  taxonomy <- data.table(taxonomy, key = "ESV")
  denovo_midas <- data.table(denovo_midas, key = "ESV")
  
  #merge by ESV
  merged_tax <- taxonomy[denovo_midas]
  
  #fill out empty entries in typestrains+SILVA with denovo MiDAS taxonomy
  merged_tax[which(Species %in% c(NA, "")), Species:=i.Species]
  merged_tax[which(Genus %in% c(NA, "")), Genus:=i.Genus]
  merged_tax[which(Family %in% c(NA, "")), Family:=i.Family]
  merged_tax[which(Order %in% c(NA, "")), Order:=i.Order]
  merged_tax[which(Class %in% c(NA, "")), Class:=i.Class]
  merged_tax[which(Phylum %in% c(NA, "")), Phylum:=i.Phylum]
  
  #order by ESV ID
  merged_tax <- merged_tax[order(as.integer(gsub("[^0-9+$]", "", ESV))), 1:8]
  
  #write out
  write_tax(merged_tax,
            file = "./output/tax_complete.csv")
  
  ##### export ESVs with SINTAX taxonomy in headers #####
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
  
  ##### export taxonomy in QIIME format #####
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
generateTaxonomy

#########################################################################################
###### NEED TO REMOVE TEMP FILES ######
mv temp/ESVs.fa output/ESVs.fa

#echoWithHeader "Removing temporary files..."
#rm -rf temp/
duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
date
echoWithHeader "Done in: $duration! Results are in the output/ folder, enjoy!"
echo "$DATA: $duration" >> benchmark.txt
