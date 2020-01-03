# About AutoTax
*AutoTax* is a linux BASH script that automatically generates de novo taxonomy from full length 16S rRNA exact sequence variants (ESVs). This allows generation of eco-system specific de novo taxonomic databases based on any environmental sample(s). It does so by combining several different software tools, listed below, into a single BASH script that otherwise only requires a single FASTA file as input. For a more detailed description of *AutoTax*, please refer to the paper [Dueholm et al, 2019](https://doi.org/10.1101/672873). *AutoTax* has only been tested on Ubuntu 18.04 LTS, but will probably run just fine on other Linux distributions as long as the required software listed below is installed.

Table of Contents
=================

   * [About AutoTax](#about-autotax)
   * [What the script does](#what-the-script-does)
   * [Installation and requirements](#installation-and-requirements)
      * [Software](#software)
      * [Database files](#database-files)
   * [Usage](#usage)
   * [Running AutoTax from a docker container (recommended)](#running-autotax-from-a-docker-container-recommended)
   * [Example data](#example-data)
   * [Generating input full-length 16S sequences](#generating-input-full-length-16s-sequences)
   * [See also](#see-also)
   * [Notes](#notes)


# What the script does
In brief, the script performs the following steps:
 - Check user input, files and folders, and check for installed R packages, installing missing ones
 
 **Generate/identify ESVs**
 
 - Orient the sequences based on the SILVA taxonomic database (usearch)
 - Dereplicate the input sequences (both strands), and determine the coverage of each unique sequence (usearch)
 - Denoise the dereplicated sequences using UNOISE3, with `minsize = 2` (usearch)
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

# Installation and requirements
As AutoTax is simply a BASH script that wraps and combines other software tools and their outputs, so there is no installation to do for the AutoTax script itself. Simply download the `autotax.sh` script by either:
```
wget https://raw.githubusercontent.com/KasperSkytte/AutoTax/master/autotax.sh
```

or clone the github repository by (make sure git is installed):
```
git clone https://github.com/KasperSkytte/AutoTax.git
cd AutoTax
```

<<<<<<< Updated upstream
Other than the standard linux tools `awk`, `grep`, and `cat` (which is included in most Linux distributions), AutoTax depends on a few other software tools, however, which need to be installed and be available in the [PATH variable](https://opensource.com/article/17/6/set-path-linux). The tools can be installed manually by refering to the documentation of the individual tools. It is recommended to run AutoTax through the docker container image based on Ubuntu linux 18.04, however, with everything pre-installed and tested (except database files), see [this section](#running-autotax-from-a-docker-container-recommended).
=======
Other than the standard linux tools `awk`, `grep`, and `cat` (which is included in most Linux distributions), AutoTax depends on a few other software tools, however, which need to be installed and be available in the [PATH variable](https://opensource.com/article/17/6/set-path-linux). The tools can be installed manually by refering to the documentation of the individual tools. It is recommended to run AutoTax through the docker container image based on Ubuntu linux 18.04, however, with everything pre-installed and tested (except database files), see the [docker section](#running-autotax-from-a-docker-container-recommended).
>>>>>>> Stashed changes

## Software
 - GNU parallel (version 20161222)
 - usearch (version 10 or later)
 - SINA (version 1.6 or later)
 - R (version 3.5 or later) with the following packages installed (the script will attempt to install if missing):
   - Biostrings (from Bioconductor, be ready for trouble if you dont have administrative rights, try manually if it fails)
   - doParallel
   - stringr (and stringi)
   - data.table
   - tidyr
   - dplyr

## Database files
Other than these software tools, SILVA and SILVA typestrains database files in both UDB and ARB format are needed. A zip file of all 4 files can be found on figshare [here](https://doi.org/10.6084/m9.figshare.9994568) (both SILVA release 132 and 138). **Make sure the paths to these files are set correctly in the [`autotax.sh`](https://github.com/KasperSkytte/AutoTax/blob/master/autotax.sh) script**. You can also use other databases, but the script is made to handle the finicky details of SILVA particularly. If you want to use other databases, you will need to adjust the script.

# Usage
Adjust the variables in the SETUP chunk at the start of the [`autotax.sh`](https://github.com/KasperSkytte/AutoTax/blob/master/autotax.sh) script to match the paths to the database files and executables.
Then simply run the script with fx `bash autotax.sh -i myseqs.fa`. Make sure the script is executable with `chmod +x autotax.sh`.
Type `bash autotax.sh -h` to show available options and version:
```
$ bash autotax.sh -h
Pipeline for extracting Exact Sequence Variants (ESV's) from full length 16S rRNA gene DNA sequences and generating de novo taxonomy
Version: 1.3.1
Options:
  -h    Display this help text.
  -i    Input FASTA file with full length DNA sequences to process (required).
  -d    FASTA file with previously processed ESV sequences.
          ESV's generated from the input sequences will then be appended to this and de novo taxonomy is rerun.
  -t    Maximum number of threads to use. Default is all available cores except 2.
```

# Running AutoTax from a docker container (recommended)
To run AutoTax through a docker container first install [Docker Engine - Community](https://docs.docker.com/install/linux/docker-ce/ubuntu/) as described there. A prebuilt image `autotax` based on Ubuntu Linux 18.04 can then be retrieved from [Docker Hub](https://hub.docker.com/) with all the required software and dependencies preinstalled (exact versions that are tested and guaranteed to work as intended):
```
sudo docker pull kasperskytte/autotax:1.0
```

Alternatively build the image manually by downloading the [Dockerfile](https://github.com/KasperSkytte/AutoTax/blob/master/docker/Dockerfile) directly from the github repository (may take 10-20 minutes):
```
git clone https://github.com/KasperSkytte/AutoTax.git
cd AutoTax
sudo docker build -t autotax:1.0 docker/
```

Then make a suitable working folder and copy the [`autotax.sh`](https://github.com/KasperSkytte/AutoTax/blob/master/autotax.sh) there together with your data and reference databases. As [usearch](http://drive5.com/usearch/) is non-free software it is not included in the image. You must buy it (the free 32-bit version is limited to 4GB memory and is doubtfully going to be sufficient) and place it in the same working folder. Please respect the [usearch software license](http://drive5.com/usearch/license64comm.html). Make sure the name of the file matches exactly that defined in the [`autotax.sh`](https://github.com/KasperSkytte/AutoTax/blob/master/autotax.sh) script (default is just "usearch11") or change the line `export usearch=$(which usearch11)` in the script accordingly. Now run AutoTax with the current working directory mounted inside the container as `/home/user/`:
```
sudo docker run -it --name autotax_run1 -v ${PWD}:/home/user/ autotax:1.0 bash autotax.sh -h
```

When running through the docker container all paths must relative to the working directory. Absolute paths (i.e. starts with `/`) won't work as the container file system is separate from the host file system. Furthermore, the output folders `temp` and `output` will be owned by root, so it's a good idea to change ownership afterwards with fx:
```
sudo chown -R $(id -u ${USER}):$(id -g ${USER}) temp/ output/
```

# Example data
To test if AutoTax is working as intended, you can use the [10.000 example sequences](https://github.com/KasperSkytte/AutoTax/blob/master/example_data/10k_fSSUs.fa) and the database files for SILVA 138 NR 99 provided on [figshare](https://doi.org/10.6084/m9.figshare.9994568), and check that the results match those in the [example_run](https://github.com/KasperSkytte/AutoTax/blob/master/example_run/output/) folder. The file `ESVs_w_sintax.fa` should be enough to check as it contains both the sequences and the taxonomy in one file. The `md5sum` of the file when using SILVA 138 NR99 should be `ad78b25419edc2cfe30c1bd28b3d6b18`. If not, feel free to post an issue.


# Generating input full-length 16S sequences
*AutoTax* is made to take input sequences obtained from the method described in [Karst et al, 2018](https://www.nature.com/articles/nbt.4045). The sequences need to be processed first using the Perl scripts in the `/fSSU-pipelines` subfolder. 

It is also possible to use full-length 16S sequences obtained from other methods such as that described in [Callahan et al, 2019](https://doi.org/10.1101/392332) based on PacBio circular consensus sequencing as long as the the error-rate is near-zero.

# See also
In the future full ribosomal operon taxonomic databases may be possible using Nanopore Sequencing and Unique Molecular Tagging, see https://www.biorxiv.org/content/10.1101/645903v2.

# Notes
If you want to use vsearch instead of usearch, feel free to adjust the script accordingly and test it. We would love to hear about the results if so, we have only used and tested usearch10 and usearch11. 
