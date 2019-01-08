#!/bin/bash
# usage: ./pipeline.sh fSSU.fa 36

##### requirements (tested with) #####
# SILVA database files (SSURef+typestrains) in arb format (use latest)
# usearch (version10)
# SINA (version 1.4.0)
# R version (version 3.5.x)
#   Note: The required R packages will install automatically, but Bioconductor can cause trouble
#   if it desides to update base R packages like MASS, mgcv, lattice and others.
#   If so run R as root and run: install.packages("BiocManager"); BiocManager::install()

##### setup #####
DATA=$1
MAX_THREADS=${2:-$((`nproc`-2))}
#path to sina executable
export sina140=$(which sina)
export usearch10=$(which usearch10)

# Provide the path the the SILVA nr99 arb-database and the typestrains database extracted from the SILVA nr99 arb-database.
silva_db="refdatabases/SILVA_132_SSURef_NR99_13_12_17_opt.arb"
silva_udb="refdatabases/SILVA_132_SSURef_Nr99_tax_silva.udb"
typestrains_db="refdatabases/SILVA132-typestrains.arb"
typestrains_udb="refdatabases/SILVA_132_SSURef_Nr99_typestrains.udb"

##### end of setup #####
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

### Orient to the same strand
# Create reference database, only needed once
#$usearch10 -makeudb_usearch refdatabases/SILVA_132_SSURef_Nr99_tax_silva.fasta -output refdatabases/SILVA_132_SSURef_Nr99_tax_silva.udb

# Run the orientation check
echoWithHeader "Orienting sequences..."
$usearch10 -quiet -orient temp/preESV_wsize.fa -db refdatabases/SILVA_132_SSURef_Nr99_tax_silva.udb -fastaout temp/preESV_wsize_oriented.fa -tabbedout temp/preESV_wsize_oriented.txt -threads 1

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
#of input sequences, and cluster_fast runs on 1 thread anyways even if set to more than 1
echoWithHeader "Clustering..."
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.987 -maxrejects 0 -centroids temp/SILVA_DeNovoSpecies.fa -uc temp/SILVA_ESV-S.txt -sortedby other
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.945 -maxrejects 0 -centroids temp/SILVA_DeNovoGenus.fa -sortedby other
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.865 -maxrejects 0 -centroids temp/SILVA_DeNovoFamily.fa -sortedby other
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.82 -maxrejects 0 -centroids temp/SILVA_DeNovoOrder.fa -sortedby other
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.785 -maxrejects 0 -centroids temp/SILVA_DeNovoClass.fa -sortedby other
$usearch10 -quiet -cluster_smallmem temp/ESVs_SILVA_trimmed_sorted.fa -id 0.75 -maxrejects 0 -centroids temp/SILVA_DeNovoPhylum.fa -sortedby other

#map per taxonomic level
echoWithHeader "Mapping each taxonomic level..."
$usearch10 -quiet -usearch_global temp/SILVA_DeNovoSpecies.fa -db temp/SILVA_DeNovoGenus.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_S-G.txt -threads $MAX_THREADS
$usearch10 -quiet -usearch_global temp/SILVA_DeNovoGenus.fa -db temp/SILVA_DeNovoFamily.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_G-F.txt -threads $MAX_THREADS
$usearch10 -quiet -usearch_global temp/SILVA_DeNovoFamily.fa -db temp/SILVA_DeNovoOrder.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_F-O.txt -threads $MAX_THREADS
$usearch10 -quiet -usearch_global temp/SILVA_DeNovoOrder.fa -db temp/SILVA_DeNovoClass.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_O-C.txt -threads $MAX_THREADS
$usearch10 -quiet -usearch_global temp/SILVA_DeNovoClass.fa -db temp/SILVA_DeNovoPhylum.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_C-P.txt -threads $MAX_THREADS

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
  read_sort_mappings <- function(input) {
    x <- data.table::fread(input,
                           sep = "\t",
                           fill = TRUE,
                           check.names = FALSE,
                           stringsAsFactors = FALSE)
    x <- x[,1:2]
    colnames(x) <- c("V1", "V2")
    x[] <- lapply(x, gsub, pattern = "\\..*$", replacement = "")
    x <- x[order(as.integer(gsub("[^0-9+$]", "", V1))),]
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
  ESV_S <- data.table::fread("./temp/SILVA_ESV-S.txt",
                             sep = "\t",
                             fill = TRUE,
                             check.names = FALSE,
                             stringsAsFactors = FALSE)
  #keep only rows with H (hits) and S (singletons) in V1, keep only V9+V10
  ESV_S <- ESV_S[V1 %in% c("H", "S"),.(V9, V10)]
  #if * in V10, replace with V9
  ESV_S <- ESV_S[,V10 := ifelse(V10 == "*", V9, V10)]
  colnames(ESV_S) <- c("ESV", "Species")
  #remove length from ESV ID's (".xxxx")
  ESV_S$Species <- gsub("\\..*$", "", ESV_S$Species)
  #order by ESV ID
  ESV_S <- ESV_S[order(as.integer(gsub("[^0-9+$]", "", ESV))),]
  
  #read and sort mappings
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
  
  #merge each taxonomic level according to the mapping results
  denovo_midas <- left_join(ESV_S, S_G, by = "Species")
  denovo_midas <- left_join(denovo_midas, G_F, by = "Genus")
  denovo_midas <- left_join(denovo_midas, F_O, by = "Family")
  denovo_midas <- left_join(denovo_midas, O_C, by = "Order")
  denovo_midas <- left_join(denovo_midas, C_P, by = "Class")
  #reorder columns
  denovo_midas <- denovo_midas[,c("ESV", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
  
  #generate denovo names per taxonomic level based on ESV ID
  denovo_midas[["Species"]] <- gsub("^[^0-9]+", "midas_s_", denovo_midas[["Species"]])
  denovo_midas[["Genus"]] <- gsub("^[^0-9]+", "midas_g_", denovo_midas[["Genus"]])
  denovo_midas[["Family"]] <- gsub("^[^0-9]+", "midas_f_", denovo_midas[["Family"]])
  denovo_midas[["Order"]] <- gsub("^[^0-9]+", "midas_o_", denovo_midas[["Order"]])
  denovo_midas[["Class"]] <- gsub("^[^0-9]+", "midas_c_", denovo_midas[["Class"]])
  denovo_midas[["Phylum"]] <- gsub("^[^0-9]+", "midas_p_", denovo_midas[["Phylum"]])
  
  #write out
  write_tax(denovo_midas,
            file = "./output/tax_midas_denovo.csv")
  
  #merge SILVA+typestrains+denovo MiDAS taxonomy
  #make data.tables
  ESV_slv_typestr_tax <- data.table(ESV_slv_typestr_tax, key = "ESV")
  denovo_midas <- data.table(denovo_midas, key = "ESV")
  
  #merge by ESV
  merged_tax <- ESV_slv_typestr_tax[denovo_midas]
  
  #fill out empty entries in typestrains+SILVA with denovo MiDAS taxonomy
  merged_tax[which(Species %in% c(NA, "")), Species:=i.Species]
  merged_tax[which(Genus %in% c(NA, "")), Genus:=i.Genus]
  merged_tax[which(Family %in% c(NA, "")), Family:=i.Family]
  merged_tax[which(Order %in% c(NA, "")), Order:=i.Order]
  merged_tax[which(Class %in% c(NA, "")), Class:=i.Class]
  merged_tax[which(Phylum %in% c(NA, "")), Phylum:=i.Phylum]
  
  #order by ESV ID
  merged_tax <- merged_tax[order(as.integer(gsub("[^0-9+$]", "", ESV))), 1:8]
  
  #fix taxa with more than one parent
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
               "./output/polyphyletics.log",)
  }
  
  #fix them
  merged_tax[, Genus := first(Genus), by = Species]
  merged_tax[, Family := first(Family), by = Genus]
  merged_tax[, Order := first(Order), by = Family]
  merged_tax[, Class := first(Class), by = Order]
  merged_tax[, Phylum := first(Phylum), by = Class]
  merged_tax[, Kingdom := first(Kingdom), by = Phylum]
  
  ##### search and replace according to a replacement file #####
  replacements <- fread("190107_replacements.txt", data.table = FALSE)
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
generateTaxonomy

#########################################################################################
###### NEED TO REMOVE TEMP FILES ######
mv temp/ESVs.fa output/ESVs.fa

#echoWithHeader "Removing temporary files..."
#rm -rf temp/
duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echoWithHeader "Done in: $duration! Results are in the ./output/ folder, enjoy!"
