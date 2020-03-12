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
                                  c("candidatus" = "Ca",
                                    "Candidatus" = "Ca",
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
  x <- x[order(as.integer(gsub("[^0-9+$]|\\..*$", "", V9))),]
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
