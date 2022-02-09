# AutoTax

*AutoTax* is a workflow that automatically generates de novo taxonomy from full length 16S rRNA amplicon sequence variants (FL-ASVs). This allows generation of ecosystem-specific de novo taxonomic databases based on any environmental sample(s). It does so by combining several different software tools, listed below, into a single BASH script that otherwise only requires a single FASTA file as input. For a more detailed description of *AutoTax*, please refer to the paper [Dueholm et al, 2020](https://doi.org/10.1128/mBio.01557-20).

Table of Contents
=================

* [What the script does](#what-the-script-does)
* [Installation and requirements](#installation-and-requirements)
   * [Required software](#required-software-tools)
   * [Required database files](#required-database-files)
* [Environment variables](#environment-variables)
* [Usage](#usage)
   * [Customization](#customization)
* [Running AutoTax from a container (recommended)](#running-autotax-from-a-container-recommended)
   * [Running getsilvadb.sh through a container](#running-getsilvadbsh-through-a-container)
   * [Important notes when running AutoTax through a container](#important-notes-when-running-autotax-through-a-container)
* [Unit tests](#unit-tests)
* [Generating input full-length 16S sequences](#generating-input-full-length-16s-sequences)
* [See also](#see-also)
* [vsearch to replace usearch](#vsearch-to-replace-usearch)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc)

# What the script does
In brief, the script performs the following steps:
 - Check user input, files and folders, and check for installed R packages, installing missing ones
 
**Generate/identify FL-ASVs**
 
 - Orient the sequences based on the SILVA taxonomic database (usearch)
 - Dereplicate the input sequences (both strands), and determine the coverage of each unique sequence (usearch)
 - Denoise the dereplicated sequences using UNOISE3, with `minsize = 2` by default (usearch)
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
 - FL-ASVs in FASTA format with DADA2 formatted taxonomy in the headers (R)
 - QIIME formatted table (R)
 - CSV files of the individual tables mentioned earlier as well as the combined, complete taxonomy for each FL-ASV (R)

# Installation and requirements
The easiest and recommended way to run AutoTax is through the official Docker container, either by using Docker (privileged, but most convenient) or Singularity (non-privileged), see [Running AutoTax from a container (recommended)](#running-autotax-from-a-container-recommended). This ensures complete reproducibility, and our unit tests have only been designed for running *AutoTax* through the container. 

Alternatively, install *AutoTax* natively by downloading the `autotax.bash` script and then make sure the required tools are installed and available in the [PATH variable](https://opensource.com/article/17/6/set-path-linux):

```
wget https://raw.githubusercontent.com/KasperSkytte/AutoTax/main/autotax.bash
```
or clone the git repository (recursively to include submodules):
```
git clone --recursive https://github.com/KasperSkytte/AutoTax.git
cd AutoTax
```

## Required software tools
 - usearch (11)
 - SINA (1.6 or later)
 - GNU parallel (20161222-1)
 - [filterShortSeqs](https://github.com/KasperSkytte/filtershortseqs/), credit goes to [Nick Green](https://github.com/nickgreensgithub/find_longest_sequences). The initial R implementation was extremely inefficient
 - R (3.5 or later) with the following packages installed (the script will attempt to install if missing):
   - Biostrings (from Bioconductor through `BiocManager::install()`)
   - doParallel
   - stringr (and stringi)
   - data.table
   - tidyr
   - dplyr
 - standard linux tools `awk`, `grep`, and `cat` (already included in most Linux distributions)

## Required database files
*AutoTax* is tailored for the SILVA database, which is also required. SILVA and SILVA typestrains database files in both UDB and ARB format are needed. A zip file with all 4 files for SILVA releases 132+138 can be found on figshare [here](https://doi.org/10.6084/m9.figshare.9994568), but won't be updated in the future. Instead use the `getsilvadb.sh` script  which will download all the required files for a chosen release version directly from https://www.arb-silva.de/, and then automagically reformat, extract typestrains, and generate UDB databases for usearch. This is perhaps also easiest through a container, see [Running getsilvadb.sh through a container](#running-getsilvadbsh-through-a-container).

Now it's important to **make sure the paths to these files are set correctly**. This is done by setting a few [environment variables](#environment-variables).

# Environment variables
Before running *AutoTax* it's important to set a few options and filepaths to the respective database files. Inspect the [autotax.bash](https://github.com/KasperSkytte/AutoTax/blob/main/autotax.bash) for defaults. This is done by setting the following environment variables in the current shell (fx using [`export`](https://www.computerhope.com/unix/bash/export.htm)):

- `silva_db`: Path to the SILVA `.arb` database file
- `silva_udb`: Path to the SILVA SSURef database file in `.udb` format
- `typestrains_udb`: Path to the typestrains database file in `.udb` format
- `denovo_prefix`: A character string which will be the prefix for de novo taxonomy, resulting in e.g. `denovo_s_23` if set to `denovo` (default)
- `denoise_minsize`: The minimum abundance of each unique input sequence. Input sequences with lower abundance than this threshold will be discarded. Passed on directly to [UNOISE3](https://drive5.com/usearch/manual/cmd_unoise3.html) during the denoise step. Set this to `1` to skip denoising, e.g. if input sequences are already pre-processed, or output from a previous autotax run etc, in which case the pipeline will fail due to 0 sequences output from this step.
- `usearch_global_threads`: Any `usearch_global` command will be split into smaller separate jobs using GNU parallel as the multithreading implementation in usearch does not scale linearly. It's much faster to run many smaller jobs. This sets the max number of threads each parallel command will use. Increase this if you lack memory.

# Usage
Make sure the script is executable with `chmod +x autotax.bash`. 

Type `bash autotax.bash -h` to show available options and version:
```
$ bash autotax.bash -h
Pipeline for extracting Full-length 16S rRNA Amplicon Sequence Variants (FL-ASVs) from full length 16S rRNA gene DNA sequences and generating de novo taxonomy
Version: 1.7.4
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

## Customization
The `autotax.bash` script essentially consists of individual functions that can be used independently. This is what makes it possible to run unit tests on a BASH script, but it also makes it possible to source the individual functions manually from another script to create custom workflows or resuming from a previous run. Simply adding `. autotax.bash` to the script won't run autotax, but will load the functions.

# Running AutoTax from a container (recommended)
To run AutoTax through a docker container first install [Docker Engine - Community](https://docs.docker.com/install/linux/docker-ce/ubuntu/) as described there. A prebuilt image `autotax` based on Ubuntu Linux 20.04 with all the required software and dependencies preinstalled (exact versions that are tested and working as intended) is provided with:
```
docker pull ghcr.io/kasperskytte/autotax:latest
```

Alternatively build the image manually from the git repo. You can also use a specific (and locked) image to ensure complete reproducibilty by using semver tags instead of just pulling the latest image every time, for instance: `ghcr.io/kasperskytte/autotax:v1.7.4`.

The image also contains the autotax git repository itself (most recent from main branch) located at `/opt/autotax/`. Now run AutoTax with the current working directory mounted inside the container at `/autotax`:
```
docker run \
  -it \
  -v ${PWD}:/autotax \
  -e denovo_prefix="midas" \
  ghcr.io/kasperskytte/autotax:latest -h
```

Running the AutoTax docker container using [Singularity](https://sylabs.io/) instead is also possible (3.9 or later to be compatible with https://ghcr.io) and is as simple as:
```
singularity run \
  --bind ${PWD}:/autotax \
  --pwd /autotax \
  --env denovo_prefix="midas" \
  docker://ghcr.io/kasperskytte/autotax:latest -h
```

Singularity has the advantage that it doesn't require elevated privileges by default like docker does. You can find a convenience script I've made to install singularity [here](https://github.com/KasperSkytte/bioscripts#install_singularitysh).

Setting the required [Environment variables](#environment-variables) is done by using the `-e` option with docker, fx `-e denovo_prefix="midas" -e denoise_minsize=4`, and similarly `--env` with singularity.

## Running getsilvadb.sh through a container
Downloading the SILVA database files automagically is easiest through a container. With docker you can easily do so by just adding `--entrypoint getsilvadb.sh`. With singularity you have to use `exec` instead of `run`:
```
$ singularity exec --bind ${PWD}:/autotax --pwd /autotax docker://ghcr.io/kasperskytte/autotax:latest getsilvadb.sh -h
INFO:    Using cached SIF image
This script downloads a desired release version of the SILVA database and makes it ready for AutoTax.
Version: 1.0
Options:
  -h    Display this help text and exit.
  -r    (required) The desired SILVA release version, fx "138.1".
  -o    Output folder. (Default: refdatabases/)
  -t    Max number of threads to use. (Default: all available except 2)
  -v    Print version and exit.
```

## Important notes when running AutoTax through a container
As [usearch](http://drive5.com/usearch/) is non-free software it is not included in the image. You must buy it or use the free 32-bit version which is preinstalled (limited to 4GB memory, most likely not sufficient). Place the executable in the same folder that is mounted inside the container and name it `usearch11` or mount it directly at `/autotax/usearch11`. Please respect the [usearch software license](http://drive5.com/usearch/license64comm.html).

When running through a container all paths must relative to the working directory. Absolute paths (i.e. starts with `/`) won't work as the container file system is separate from the host file system. 

Furthermore, the output folders `temp` and `output` will likely be owned by root if you haven't done the [Post-installation steps for Linux](https://docs.docker.com/engine/install/linux-postinstall/), so it's a good idea to either run docker with appropriate user ID mapping (`--user` option), or adjust the ownership afterwards with fx:

```
sudo chown -R $(id -u ${USER}):$(id -g ${USER}) temp/ output/
```

# Unit tests
AutoTax is being unit tested by the [Bash Automated Testing System](https://github.com/bats-core/bats-core). To run the tests, preferably before running with your own data, you can do so with the `autotax.bash -b` argument. The test result is printed to the terminal as well as a log file `test_result.log`. If you want to run through docker, you can run the tests properly with the following command:

```
$ docker run -it --rm --name autotax -v ${PWD}:/autotax ghcr.io/kasperskytte/autotax:latest -b
1..33
ok 1 Variable set: version
ok 2 silva_db database file
ok 3 silva_udb database file
ok 4 typestrains_udb database file
ok 5 Variable set: denovo_prefix
ok 6 Variable set: MAX_THREADS
ok 7 usearch11 in $PATH
ok 8 sina in $PATH
ok 9 R in $PATH
ok 10 Rscript in $PATH
ok 11 Check installed R packages
ok 12 Echo with timestamp
ok 13 Shell is BASH
ok 14 Error if temp/ folder exists
ok 15 Error if output/ folder exists
ok 16 Check input data
ok 17 Step: Orient
ok 18 Step: Dereplication
ok 19 Step: Denoise
ok 20 Step: Find longest and rename
ok 21 Step (optional): Add additional FLASVs to DB
ok 22 Step: Global alignment against SILVA
ok 23 Step: Trim and strip alignment
ok 24 Step: Sort FLASVs by ID (i.e. highest coverage)
ok 25 Step: Obtaining the taxonomy of the best hit in the SILVA database
ok 26 Step: Obtaining the taxonomy of species (>98.7% id) in the SILVA typestrains database
ok 27 Step: Cluster at species level
ok 28 Step: Cluster at genus level
ok 29 Step: Cluster at family level
ok 30 Step: Cluster at order level
ok 31 Step: Cluster at class level
ok 32 Step: Cluster at phylum level
ok 33 Step: Merge and output taxonomy
```

The exact docker command above is being used for testing the main branch on [https://github.com/kasperskytte/autotax](https://github.com/kasperskytte/autotax). The latest test log of the main branch can always be seen [here](https://github.com/KasperSkytte/AutoTax/blob/main/test_result.log)).

Beware that these tests are done by verifying that the output from the individual steps are identical to the output from a previous, manually verified, run, and that this is specific to the particular version of SILVA used (currently release 138.1).

# Generating input full-length 16S sequences
*AutoTax* is originally made to take input sequences obtained from the method described in [Karst et al, 2018](https://www.nature.com/articles/nbt.4045). The sequences may need to be processed first using the Perl scripts in the `/fSSU-pipelines` subfolder.

It is also possible to use full-length 16S sequences obtained from other methods such as that described in [Callahan et al, 2019](https://doi.org/10.1101/392332) based on PacBio circular consensus sequencing as long as the the error-rate is near-zero and primers are stripped.

# See also
In the future full ribosomal operon taxonomic databases may be possible using Nanopore Sequencing and Unique Molecular Tagging, see https://www.biorxiv.org/content/10.1101/645903v2.

# vsearch to replace usearch
A vsearch version of AutoTax is work in (slow) progress and is available at the development branch. But bear in mind that we have only used and tested usearch10 and usearch11 and the results presented in the paper has been done with usearch11, so results will likely be different compared to when using usearch (though still reliable in its own way, but don't mix output from both usearch and vsearch). `vsearch` is also installed in the docker image by default. CONTRIBUTIONS ARE WELCOME. I'm just a PhD student doing many other projects now.
