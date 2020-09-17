# About AutoTax

[![Docker Pulls](https://img.shields.io/docker/pulls/kasperskytte/autotax)][hub]

*AutoTax* is a linux BASH script that automatically generates de novo taxonomy from full length 16S rRNA amplicon sequence variants (FL-ASVs). This allows generation of eco-system specific de novo taxonomic databases based on any environmental sample(s). It does so by combining several different software tools, listed below, into a single BASH script that otherwise only requires a single FASTA file as input. For a more detailed description of *AutoTax*, please refer to the paper [Dueholm et al, 2019](https://doi.org/10.1101/672873). *AutoTax* has only been tested on Ubuntu 18.04 LTS, but will probably run just fine on other Linux distributions as long as the required software listed below is installed.

Table of Contents
=================

   * [About AutoTax](#about-autotax)
   * [Table of Contents](#table-of-contents)
   * [What the script does](#what-the-script-does)
   * [Installation and requirements](#installation-and-requirements)
      * [Software](#software)
      * [Database files](#database-files)
   * [Usage](#usage)
   * [Running AutoTax from a docker container (recommended)](#running-autotax-from-a-docker-container-recommended)
      * [Important notes when running through docker container](#important-notes-when-running-through-docker-container)
   * [Unit tests](#unit-tests)
   * [Generating input full-length 16S sequences](#generating-input-full-length-16s-sequences)
   * [See also](#see-also)
   * [Notes](#notes)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc)

# What the script does
In brief, the script performs the following steps:
 - Check user input, files and folders, and check for installed R packages, installing missing ones
 
 **Generate/identify FL-ASVs**
 
 - Orient the sequences based on the SILVA taxonomic database (usearch)
 - Dereplicate the input sequences (both strands), and determine the coverage of each unique sequence (usearch)
 - Denoise the dereplicated sequences using UNOISE3, with `minsize = 2` (usearch)
 - Remove all sequences that match exactly (100% identity) with other, but longer sequences (R)
 - Sort the sequences based on coverage, and rename the sequences in order of occurence, in the format `FLASVx.length`, e.g. `FLASV123.1410` (R)
 - If desired, update an existing FL-ASV database (FASTA file) by matching the generated FL-ASVs to the database, replacing identical FL-ASVs with longer sequences if any, and adding the new ones to the end of the FASTA file, renamed to continue numbering from the database (R)
 
 **Generate de novo taxonomy**
 
 - Perform a multiple sequence alignment of the FL-ASVs with both the SILVA and SILVA typestrains databases using SINA, then trim, strip gaps, format, and sort based on FL-ASV IDs (multithreading doesn't always preserve ordering) (SINA+awk+R)
 - Assign taxonomy to that of the best hit in both the SILVA and SILVA typestrains databases (usearch)
 - Cluster the FL-ASVs at different identity thresholds each corresponding to a taxonomic level and use the FL-ASV ID of the cluster centroids as a de novo placeholder name at each level (usearch, thresholds from [Yarza et al, 2014](https://www.nature.com/articles/nrmicro3330))
 - Reformat the output from the last 2 steps into 3 separate tables where each column contains the taxonomy at each taxonomic level (Kingdom->Species) of each FL-ASV (R)  
 - Merge the 3 tables so that the de novo taxonomy fills in where the assigned taxonomy based on SILVA and SILVA typestrains are below the taxonomic thresholds (R)
 - Manually curate the taxonomy based on a replacement file if any (R)
 
 **Output the taxonomy in the following formats:**
 - FL-ASVs in FASTA format with usearch SINTAX formatted taxonomy in the headers (R)
 - QIIME formatted table (R)
 - CSV files of the individual tables mentioned earlier as well as the combined, complete taxonomy for each FL-ASV (R)

# Installation and requirements
As AutoTax is simply a BASH script that wraps and combines other software tools and their outputs, so there is no installation to do for the AutoTax script itself. Simply download the `autotax.bash` script by either:
```
wget https://raw.githubusercontent.com/KasperSkytte/AutoTax/master/autotax.bash
```

or clone the github repository by (make sure git is installed):
```
git clone https://github.com/KasperSkytte/AutoTax.git
cd AutoTax
```

Other than the standard linux tools `awk`, `grep`, and `cat` (which is included in most Linux distributions), AutoTax depends on a few other software tools, however, which need to be installed and be available in the [PATH variable](https://opensource.com/article/17/6/set-path-linux). The tools can be installed manually by refering to the documentation of the individual tools. It is recommended to run AutoTax through the docker container image based on Ubuntu linux 18.04, however, with everything pre-installed and tested (except database files), see [this section](#running-autotax-from-a-docker-container-recommended).

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
Other than these software tools, SILVA and SILVA typestrains database files in both UDB and ARB format are needed. A zip file of all 4 files can be found on figshare [here](https://doi.org/10.6084/m9.figshare.9994568) (both SILVA release 132 and 138). **Make sure the paths to these files are set correctly in the [`autotax.bash`](https://github.com/KasperSkytte/AutoTax/blob/master/autotax.bash) script**. You can also use other databases, but the script is made to handle the finicky details of SILVA particularly. If you want to use other databases, you will need to adjust the script.

For SILVA version 138 this can be done from a shell by the following commands:
```
wget https://ndownloader.figshare.com/files/22790396 -O SILVA138_NR99.zip
unzip SILVA138_NR99.zip -d refdatabases/
```

# Usage
Adjust the variables in the SETUP chunk at the start of the [`autotax.bash`](https://github.com/KasperSkytte/AutoTax/blob/master/autotax.bash) script to match the paths to the database files and executables. If you downloaded SILVA138 using the link above, you don't have to adjust anything if you create a folder named `refdatabases` and extract all the files into the folder. Make sure the script is executable with `chmod +x autotax.bash`.
Type `bash autotax.bash -h` to show available options and version:
```
$ bash autotax.bash -h
Pipeline for extracting Full-length 16S rRNA Amplicon Sequence Variants (FL-ASVs) from full length 16S rRNA gene DNA sequences and generating de novo taxonomy
Version: 1.5.4
Options:
  -h    Display this help text and exit.
  -i    Input FASTA file with full length DNA sequences to process (required).
  -c    Cluster the resulting FL-ASVs at 99% (before generating de novo taxonomy),
          do chimera filtering on the clusters, and then add them on top in the same way as when using -d.
  -d    FASTA file with previously processed FL-ASV sequences.
          FL-ASVs generated from the input sequences will then be appended to this and de novo taxonomy is rerun.
  -t    Maximum number of threads to use. Default is all available cores except 2.
  -b    Run all BATS unit tests to assure everything is working as intended (requires git).
  -v    Print version and exit.
```

Using the example data in `/test/example_data/` a usage example would be:
`bash autotax.bash -i test/example_data/10k_fSSUs.fa -t 20`.

The main output files can then be found in the `output/` folder and all intermediate files along the way in `temp/`.

# Running AutoTax from a docker container (recommended)
To run AutoTax through a docker container first install [Docker Engine - Community](https://docs.docker.com/install/linux/docker-ce/ubuntu/) as described there. A prebuilt image `autotax` based on Ubuntu Linux 18.04 can then be retrieved from [Docker Hub](https://hub.docker.com/) with all the required software and dependencies preinstalled (exact versions that are tested and guaranteed to work as intended):
```
sudo docker pull kasperskytte/autotax:latest
```

Alternatively build the image manually by downloading the [Dockerfile](https://github.com/KasperSkytte/AutoTax/blob/master/docker/Dockerfile) directly from the github repository (may take 10-20 minutes):
```
git clone https://github.com/KasperSkytte/AutoTax.git
cd AutoTax
sudo docker build -t kasperskytte/autotax:latest docker/
```

The image also contains the autotax github repository itself (most recent from master branch) in `/opt/autotax/`. Now run AutoTax with the current working directory mounted inside the container as `/autotax`:
```
sudo docker run -it --rm --name autotax -v ${PWD}:/autotax kasperskytte/autotax:latest -h
```

## Important notes when running through docker container
As [usearch](http://drive5.com/usearch/) is non-free software it is not included in the image. You must buy it or use the free 32-bit version (limited to 4GB memory and is doubtfully going to be sufficient, but you are welcome to try) and place the executable in the same folder that is mounted inside the container and name it `usearch11`. Please respect the [usearch software license](http://drive5.com/usearch/license64comm.html).

By default the [`autotax.bash`](https://github.com/KasperSkytte/AutoTax/blob/master/autotax.bash) script included in the image is executed, which assumes you have extracted the SILVA138 database (most recent as of the time of writing) into a folder named `refdatabases` in the current working directory as described in [Database files](#database-files). If you wish to use a different version you need to adjust the paths in the script itself, hence you must also copy the [`autotax.bash`](https://github.com/KasperSkytte/AutoTax/blob/master/autotax.bash) script into the current working folder, adjust the paths, and run that instead of that included in the image.

When running through the docker container all paths must relative to the working directory. Absolute paths (i.e. starts with `/`) won't work as the container file system is separate from the host file system. Furthermore, the output folders `temp` and `output` will be owned by root, so it's a good idea to change ownership afterwards with fx:
```
sudo chown -R $(id -u ${USER}):$(id -g ${USER}) temp/ output/
```

# Unit tests
AutoTax is being unit tested by the [Bash Automated Testing System](https://github.com/bats-core/bats-core). To run the tests, preferably before running with your own data, you can do so with the `autotax.bash -b` argument. This requires you to run from the root of a clone of the AutoTax git repository as several additional test files are needed. The test result is printed to the terminal as well as a log file `test_result.log`. If you want to run through docker, you can run the tests properly with the following command:

```
sudo docker run -it --rm --name autotax -v ${PWD}:/autotax kasperskytte/autotax:latest -b
```

The exact docker command above is being used for testing the master branch on [https://github.com/kasperskytte/autotax](https://github.com/kasperskytte/autotax), the latest test log of the master branch can be seen [here](https://github.com/KasperSkytte/AutoTax/blob/master/test_result.log)).

# Generating input full-length 16S sequences
*AutoTax* is made to take input sequences obtained from the method described in [Karst et al, 2018](https://www.nature.com/articles/nbt.4045). The sequences need to be processed first using the Perl scripts in the `/fSSU-pipelines` subfolder. 

It is also possible to use full-length 16S sequences obtained from other methods such as that described in [Callahan et al, 2019](https://doi.org/10.1101/392332) based on PacBio circular consensus sequencing as long as the the error-rate is near-zero.

# See also
In the future full ribosomal operon taxonomic databases may be possible using Nanopore Sequencing and Unique Molecular Tagging, see https://www.biorxiv.org/content/10.1101/645903v2.

# Notes
If you want to use vsearch instead of usearch, feel free to adjust the script accordingly and test it. We would love to hear about the results if so, we have only used and tested usearch10 and usearch11. 
