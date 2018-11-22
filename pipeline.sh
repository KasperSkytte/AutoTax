#!/bin/bash
# usage: ./pipeline.sh fSSU.fa 36
DATA=$1
THREADS=${2:-36}

# Provide the path the the SILVA nr99 arb-database and the typestrain database extracted from the SILVA nr99 arb-database.
silva_db="refdatabases/SILVA_132_SSURef_NR99_13_12_17_opt.arb"
typestrain_db="refdatabases/SILVA132-typestrains.arb"

echoWithHeader() {
  header() {
    echo "################################################################################"
  }
  header
  echo $1
  header
}

echoWithHeader "Running fSSU taxonomy pipeline..."
date
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
  local JOBS=${5:-$THREADS}

  echo "Running SINA in parallel using GNU parallel..."
  FILE_SIZE=$(stat --printf="%s" ${DATA})
  BLOCK_SIZE=$((FILE_SIZE / JOBS))
  if [ ${BLOCK_SIZE} -lt 1000 ] ; then
    BLOCK_SIZE=1000
  fi
  if [ ${BLOCK_SIZE} -gt 2000000 ] ; then
    BLOCK_SIZE=2000000
  fi
  echo "Splitting data into blocks with size: $BLOCK_SIZE"

  mkdir parallel_temp
  cat ${DATA} | awk -F' ' -vOFS=' ' '/^>/ {$1=$1 " " substr($1,2,30)} {print}' |\
    parallel --progress -j $JOBS --block $BLOCK_SIZE --recstart '>' \
    --pipe "cat > parallel_temp/temp_block_{#}.fa; sina -i parallel_temp/temp_block_{#}.fa \
    --intype fasta -o parallel_temp/${OUTPUTID}_block_{#}_aligned.fa --meta-fmt \
    csv nearest nuc achieved_idty name -r $DB --search --search-db $DB \
    --lca-fields $LCA_FIELDS --search-min-sim 0.5 --search-max-result 1"

  echo ""
  echo "Concatenating, sorting, reformatting, and trimming data..."
  #sequences
  cat parallel_temp/${OUTPUTID}_block_*.fa > temp/${OUTPUTID}_aligned.fa
  sed -i 's/ ESV[0-9]*//g' temp/${OUTPUTID}_aligned.fa
  awk '!/^>/ {$0=substr($0, 1036, 41788)}1' temp/${OUTPUTID}_aligned.fa > temp/${OUTPUTID}_aligned_trimmed.fa
  rm temp/${OUTPUTID}_aligned.fa
  cat temp/${OUTPUTID}_aligned_trimmed.fa | tr -d "-" | tr -d "." > temp/${OUTPUTID}_aligned_trimmed2.fa
  rm temp/${OUTPUTID}_aligned_trimmed.fa

  #stats
  head -n 1 parallel_temp/${OUTPUTID}_block_1_aligned.fa.csv > temp/${OUTPUTID}_stats.csv
  find parallel_temp/ -name "${OUTPUTID}_block_*.csv" | xargs -n 1 tail -n +2 >> temp/${OUTPUTID}_stats.csv
  sed -i 's/ ESV[0-9]*//g' temp/${OUTPUTID}_stats.csv

  #sort in R
  #careful with the order of arguments passed on to R
  R --slave --args temp/${OUTPUTID}_aligned_trimmed2.fa temp/${OUTPUTID}_stats.csv << EOF
	#extract passed args from shell script
	files <- commandArgs(trailingOnly = TRUE)

	##### check packages, install if not available #####
	suppressPackageStartupMessages({
		#Biostrings (and BiocManager which is used to install Biostrings)
		if(!require("Biostrings")) {
		  if(!require("BiocManager")) {
		    install.packages("BiocManager")
		  }
		  BiocManager::install("Biostrings", 
		                       update = FALSE, 
		                       ask = FALSE)
		  require("Biostrings")
		}

		#stringr
		if(!require("stringr")) {
		  install.packages("stringr")
		  require("stringr")
		}

		#data.table
		if(!require("data.table")) {
		  install.packages("data.table")
		  require("data.table")
		}
	})

	#reorder FASTA
	algn.fa <- Biostrings::readBStringSet(files[[1]])
	algn.fa <- algn.fa[order(as.numeric(unlist(stringr::str_extract_all(names(algn.fa), "[:digit:]+$"), use.names = FALSE)))]
	Biostrings::writeXStringSet(x = algn.fa, filepath = files[[1]])

	#reorder stats
	stats <- fread(files[[2]],
	               sep = ",",
	               sep2 = ";",
	               fill = TRUE,
	               check.names = FALSE,
	               stringsAsFactors = FALSE)
	stats[, name := as.numeric(gsub("ESV", "", name))]
	stats <- stats[order(name),]
	stats[, name := paste0("ESV", name)]
	fwrite(stats,
	       file = files[[2]],
	       sep = ",",
	       row.names = FALSE,
	       col.names = TRUE)
EOF

  echo "Removing temporary files and killing running ARB PT servers..."
  rm -rf parallel_temp
  pgrep arb_pt_server
  pkill arb_pt_server
}

#########################################################################################
echoWithHeader "Finding unique sequences occuring at least 2 times..."
usearch10 -fastx_uniques $DATA -fastaout temp/all_ESVs.fa -sizeout -minuniquesize 2 -strand both -threads $THREADS

echoWithHeader "Clustering sequences and finding representative sequences (cluster centroids)..."
usearch10 -cluster_smallmem temp/all_ESVs.fa -id 1 -maxrejects 0 -centroids temp/ESV_centroids.fa -relabel ESV -sortedby other

### Orient to the same strand
# Create reference database from the rdp_train_set_16.fa
#usearch10 -makeudb_usearch refdatabases/SILVA_132_SSURef_Nr99_tax_silva.fasta -output refdatabases/SILVA_132_SSURef_Nr99_tax_silva.udb
# Run the orientation check
#usearch10 -orient temp/ESV_centroids.fa -db refdatabases/SILVA_132_SSURef_Nr99_tax_silva.udb -fastaout temp/ESV_centroids_oriented.fa -tabbedout temp/ESV_centroids_oriented.txt -threads $THREADS

#########################################################################################
# Example command: ./PrepareESV.sh all_new_fSSUs.fa oldESV.fa
#usearch10 -search_global new_ESV.fa -db $2 -id 1 -strand plus -maxrejects 0 -notmatched new_unique_ESVs.fa -threads $THREADS
#cat $2 new_unique_ESVs.fa > updated_ESVs_temp.fa
#usearch10 -fastx_relabel updated_ESVs_temp.fa -prefix ESV -fastaout ./ESV/updated_ESVs.fa -threads $THREADS
#rm new_unique_ESVs.fa
#rm updated_ESVs_temp.fa

#rm all_ESVs.fa

#########################################################################################
#typestrains
echoWithHeader "Aligning centroids with typestrains database using SINA..."
sina_align_par temp/ESV_centroids.fa ESVs_typestrains $typestrain_db tax_ltp_name $THREADS

#########################################################################################
#SILVA
echoWithHeader "Aligning centroids with SILVA database using SINA..."
sina_align_par temp/ESV_centroids.fa ESVs_SILVA $silva_db tax_slv 16

echoWithHeader "Clustering and generating denovo taxonomy..."
usearch10 -cluster_smallmem temp/ESVs_SILVA_aligned_trimmed2.fa -id 0.987 -maxrejects 0 -centroids temp/SILVA_DeNovoSpecies.fa -uc temp/SILVA_ESV-S.txt -sortedby other
usearch10 -cluster_smallmem temp/ESVs_SILVA_aligned_trimmed2.fa -id 0.945 -maxrejects 0 -centroids temp/SILVA_DeNovoGenus.fa -sortedby other
usearch10 -cluster_smallmem temp/ESVs_SILVA_aligned_trimmed2.fa -id 0.865 -maxrejects 0 -centroids temp/SILVA_DeNovoFamily.fa -sortedby other
usearch10 -cluster_smallmem temp/ESVs_SILVA_aligned_trimmed2.fa -id 0.82 -maxrejects 0 -centroids temp/SILVA_DeNovoOrder.fa -sortedby other
usearch10 -cluster_smallmem temp/ESVs_SILVA_aligned_trimmed2.fa -id 0.785 -maxrejects 0 -centroids temp/SILVA_DeNovoClass.fa -sortedby other
usearch10 -cluster_smallmem temp/ESVs_SILVA_aligned_trimmed2.fa -id 0.75 -maxrejects 0 -centroids temp/SILVA_DeNovoPhylum.fa -sortedby other

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

echoWithHeader "Mapping per taxonomic level..."
usearch10 -usearch_global temp/SILVA_DeNovoSpecies.fa -db temp/SILVA_DeNovoGenus.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_S-G.txt -threads $THREADS
usearch10 -usearch_global temp/SILVA_DeNovoGenus.fa -db temp/SILVA_DeNovoFamily.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_G-F.txt -threads $THREADS
usearch10 -usearch_global temp/SILVA_DeNovoFamily.fa -db temp/SILVA_DeNovoOrder.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_F-O.txt -threads $THREADS
usearch10 -usearch_global temp/SILVA_DeNovoOrder.fa -db temp/SILVA_DeNovoClass.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_O-C.txt -threads $THREADS
usearch10 -usearch_global temp/SILVA_DeNovoClass.fa -db temp/SILVA_DeNovoPhylum.fa -maxrejects 0 -strand plus -id 0.5 -blast6out temp/SILVA_C-P.txt -threads $THREADS


#########################################################################################
echoWithHeader "Merging taxonomy..."
#Use R to fix output
R --slave << EOF
  ##### check packages, install if not available #####
  suppressPackageStartupMessages({
    #data.table 
    if(!require("data.table")) {
      install.packages("data.table")
      require("data.table")
    }
    #tidyr
    if(!require("tidyr")) {
      install.packages("tidyr")
      require("tidyr")
    }
    #stringr
    if(!require("stringr")) {
      install.packages("stringr")
      require("stringr")
    }
    #stringr
    if(!require("dplyr")) {
      install.packages("dplyr")
      require("dplyr")
    }
  })
  
  ##### FUNCTIONS #####
  ## read/write taxonomy ##
  read_tax <- function(input, tax_field) {
    tax <- fread(input,
                 sep = ",",
                 sep2 = ";",
                 fill = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE)
    tax <- tax[, .(full_name, align_ident_slv, get(tax_field))]
    colnames(tax)[3] <- "tax"
    tax[,tax := stringr::str_replace_all(tax, 
                                         c("candidatus " = "candidatus_",
                                           "Candidatus " = "Candidatus_", 
                                           ";$" = "",
                                           "uncultured" = "",
                                           "unknown" = "",
                                           "incertae sedis" = ""))]
    invisible(tax)
  }
  
  write_tax <- function(tax, file) {
    fwrite(tax,
           file = file,
           sep = ",",
           row.names = FALSE,
           col.names = TRUE)
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
    x <- x[order(as.numeric(gsub("ESV", "", V1))),]
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
  ESV_typestrain_algn <- read_tax(input = "./temp/ESVs_typestrains_stats.csv",
                                  tax_field = "lca_tax_ltp_name")
  ESV_typestrain_algn <- suppressWarnings(
    separate(ESV_typestrain_algn,
             col = "tax",
             into = c("Genus", "Species"),
             sep = " "))
  ESV_typestrain_algn[which(align_ident_slv < 98.7), Species := ""]
  ESV_typestrain_algn[which(align_ident_slv < 94.5), Genus := ""]
  write_tax(tax = ESV_typestrain_algn,
            file = "./output/ESVs_typestrains_stats.csv")
  
  #SILVA
  ESV_SILVA_algn <- read_tax(input = "./temp/ESVs_SILVA_stats.csv",
                             tax_field = "lca_tax_slv")
  ESV_SILVA_algn <- suppressWarnings(
    separate(ESV_SILVA_algn,
             col = "tax",
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
             sep = ";"))
  ESV_SILVA_algn[which(align_ident_slv < 94.5), Genus := ""]
  ESV_SILVA_algn[which(align_ident_slv < 86.5), Family := ""]
  ESV_SILVA_algn[which(align_ident_slv < 82.0), Order := ""]
  ESV_SILVA_algn[which(align_ident_slv < 78.5), Class := ""]
  ESV_SILVA_algn[which(align_ident_slv < 75.0), Phylum := ""]
  write_tax(tax = ESV_SILVA_algn,
            file = "./output/ESVs_SILVA_stats.csv")
  
  #combine typestrains+SILVA
  taxonomy <- left_join(ESV_SILVA_algn[,-2], ESV_typestrain_algn[,-2], by = c("full_name", "Genus"))
  colnames(taxonomy)[1] <- "ESV"
  write_tax(taxonomy,
            file = "./output/tax_slv_typestr.csv")
  
  ##### Fix mappings + denovo taxonomy #####
  ESV_S <- fread("./temp/SILVA_ESV-S.txt",
                 sep = "\t",
                 fill = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE)
  #keep only rows with H and S in V1, keep only V9+V10
  ESV_S <- ESV_S[V1 %in% c("H", "S"),.(V9, V10)]
  #if * in V10, replace with V9
  ESV_S <- ESV_S[,V10 := ifelse(V10 == "*", V9, V10)]
  colnames(ESV_S) <- c("ESV", "Species")
  ESV_S <- ESV_S[,V3 := gsub("ESV", "", ESV)]
  ESV_S <- ESV_S[order(as.numeric(V3)),1:2]
  
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
  
  denovo_midas <- left_join(ESV_S, S_G, by = "Species")
  denovo_midas <- left_join(denovo_midas, G_F, by = "Genus")
  denovo_midas <- left_join(denovo_midas, F_O, by = "Family")
  denovo_midas <- left_join(denovo_midas, O_C, by = "Order")
  denovo_midas <- left_join(denovo_midas, C_P, by = "Class")
  denovo_midas <- denovo_midas[,c("ESV", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
  denovo_midas[["Species"]] <- gsub("ESV", "midas_s_", denovo_midas[["Species"]])
  denovo_midas[["Genus"]] <- gsub("ESV", "midas_g_", denovo_midas[["Genus"]])
  denovo_midas[["Family"]] <- gsub("ESV", "midas_f_", denovo_midas[["Family"]])
  denovo_midas[["Order"]] <- gsub("ESV", "midas_o_", denovo_midas[["Order"]])
  denovo_midas[["Class"]] <- gsub("ESV", "midas_c_", denovo_midas[["Class"]])
  denovo_midas[["Phylum"]] <- gsub("ESV", "midas_p_", denovo_midas[["Phylum"]])
  
  write_tax(denovo_midas,
            file = "./output/tax_midas_denovo.csv")
  
  taxonomy <- data.table(taxonomy, key = "ESV")
  denovo_midas <- data.table(denovo_midas, key = "ESV")
  merged_tax <- taxonomy[denovo_midas] #merge by ESV
  merged_tax[which(Species %in% c(NA, "")), Species:=i.Species]
  merged_tax[which(Genus %in% c(NA, "")), Genus:=i.Genus]
  merged_tax[which(Family %in% c(NA, "")), Family:=i.Family]
  merged_tax[which(Order %in% c(NA, "")), Order:=i.Order]
  merged_tax[which(Class %in% c(NA, "")), Class:=i.Class]
  merged_tax[which(Phylum %in% c(NA, "")), Phylum:=i.Phylum]
  merged_tax[which(Kingdom %in% c(NA, "")), Kingdom:=i.Kingdom]
  merged_tax <- merged_tax[order(as.numeric(gsub("ESV", "", ESV))), 1:8]
  
  write_tax(merged_tax,
            file = "./output/tax_complete.csv")
EOF

#########################################################################################
###### NEED TO REMOVE TEMP FILES ######

#echo "Removing temporary files..."
#rm -rf temp/
duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
date
echoWithHeader "Done in: $duration! Results are in the output/ folder, enjoy!"
echo "$DATA: $duration" >> benchmark.txt