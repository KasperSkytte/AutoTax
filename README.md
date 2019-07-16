# AutoTax
*AutoTax* is a linux BASH script that automatically generates de novo taxonomy from full length 16S rRNA exact sequence variants (ESVs). This allows generation of eco-system specific de novo taxonomic databases based on any environmental sample(s). It does so by combining several different software tools, listed below, into a single BASH script that otherwise only requires a single FASTA file as input. For a more detailed description of *AutoTax*, please refer to the paper [Dueholm et al, 2019]( https://doi.org/10.1101/672873). *AutoTax* has only been tested on Ubuntu 18.04 LTS, but will probably run just fine on other Linux distributions as long as the required software listed below is installed.

## Usage
Adjust the variables in the SETUP chunk at the start of the `autotax.sh` script to match the paths to the database files and executables.
Then simply run the script with fx `bash autotax.sh -i myseqs.fa`. Make sure the file is executable with `chmod +x autotax.sh`.
Type `bash autotax.sh -h` to show available options.

## What the script does
In brief, the script performs the following steps:
 - Check user input, files and folders, and check for installed R packages, installing missing ones
 
 **Generate/identify ESVs**
 
 - Dereplicate the input sequences (both strands), and determine the coverage of each unique sequence (usearch)
 - Orient the sequences based on the SILVA taxonomic database (usearch)
 - Remove all sequences that match exactly (100% identity) with other, but longer sequences (R)
 - Sort the sequences based on coverage, and rename the sequences in order of occurence, in the format `ESVx.length`, e.g. `ESV123.1410` (R)
 - If desired, update an existing ESV database (FASTA file) by matching the generated ESVs to the database, replacing identical ESVs with longer sequences if any, and adding the new ones to the end of the FASTA file, renamed to continue numbering from the database (R)
 
 **Generate de novo taxonomy**
 
 - Perform a multiple sequence alignment of the ESVs with both the SILVA and SILVA typestrains databases using SINA, then trim, strip gaps, format, and sort based on ESV IDs (multithreading doesn't always preserve ordering) (SINA+awk+R)
 - Assign taxonomy to that of the best hit in both the SILVA and SILVA typestrains databases (usearch)
 - Cluster the ESVs at different identity thresholds each corresponding to a taxonomic level and use the ESV ID of the cluster centroids as a de novo placeholder name at each level (usearch, thresholds from [Yarza et al, 2014](https://www.nature.com/articles/nrmicro3330))
 - Reformat the output from the last 2 steps into 3 separate tables where each column contains the taxonomy at each taxonomic level (Kingdom->Species) of each ESV (R)  
 - Merge the 3 tables so that the de novo taxonomy fills in where the assigned taxonomy based on SILVA and SILVA typestrains are below the taxonomic thresholds (R)
 - Manually curate the taxonomy based on a replacement file if any (R)
 
 **Output the taxonomy in the following formats:**
 - ESVs in FASTA format with usearch SINTAX formatted taxonomy in the headers (R)
 - QIIME formatted table (R)
 - CSV files of the individual tables mentioned earlier as well as the combined, complete taxonomy for each ESV (R)


## Requirements
 - Linux OS including the standard tools `awk`, `grep`, and `cat`, which are included in most Linux distribution, and BASH.
 - GNU parallel (version 20161222)
 - usearch (version 10 or later)
 - SINA (version 1.6 or later)
 - R (version 3.5 or later) with the following packages installed (the script will try to install if missing):
   - Biostrings (from Bioconductor, be ready for trouble if you dont have administrative rights, try manually if it fails)
   - doParallel
   - stringr (and stringi)
   - data.table
   - tidyr
   - dplyr

Other than these software tools, SILVA and SILVA typestrains database files in both UDB and ARB format are needed. You can also use other databases, but the script is made to handle the finicky details of SILVA particularly. If you want to use other databases, you will need to adjust the script.

## Generating input full-length 16S sequences
*AutoTax* is made to take input sequences obtained from the method described in [Karst et al, 2018](https://www.nature.com/articles/nbt.4045). The sequences need to be processed first using the Perl scripts in the `/fSSU-pipelines` subfolder. 

It is also possible to use full-length 16S sequences obtained from other methods such as that described in [Callahan et al, 2019](https://doi.org/10.1101/392332) based on PacBio circular consensus sequencing as long as the the error-rate is near-zero.

## See also
In the future full ribosomal operon taxonomic databases may be possible using Nanopore Sequencing and Unique Molecular Tagging, see https://www.biorxiv.org/content/10.1101/645903v2.
