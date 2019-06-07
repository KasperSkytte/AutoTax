#!/bin/bash

#####################
#### Description ####
#####################

USAGE="$(basename "$0") [-h] [-1 fastq -2 fastq -b textfile -a path -d alignment -c value -t] -- fSSU pipeline: Extraction of readbins, de-novo assembly and trimming.

where:
    -h  Show this help text
    -1  Read1 fastq file
    -2  Read2 fastq file
    -b  Sample barcodes
    -a  Path to CLC Assembly Cell installation
    -d  Path to aligned SSU reference database in fasta format with
        50000 alignment positions.
    -c  Number of CPUs to use
    -t  Flag to run pipeline on test data

Run pipeline on test data:
./fSSU-pipeline-DNA_v1.2.h -t
./fSSU-pipeline-DNA_v1.2.h -1 testdata/mock_R1.fq -2 testdata/mock_R2.fq -b testdata/barcodes.txt -a CLC_assembly_cell -d testdata/mock16S.alignment -c 2

Run pipleline on real data:
./fSSU-pipeline-DNA_v1.2.sh R1.fq -2 R2.fq -b barcodes.txt -a CLC_assembly_cell -d SILVA_132_SSURef_Nr99_tax_silva_full_align_trunc.fasta -c 2

Requirements - if your programs are named differently, then change the reference in the respective bash scripts
mothur (version 1.37.6; https://www.mothur.org/wiki/Download_mothur)
parallel (version 20150422; https://www.gnu.org/software/parallel/)
CLC assembly cell (version 5.0.3; https://www.qiagenbioinformatics.com/products/clc-assembly-cell/)
perl (v5.18.2; https://www.perl.org/get.html)
"

###################
#### Arguments ####
###################
EXTRACT_SCRIPT=./scripts/fSSU_extract_27f1391r_v1.1.pl
ASSEMBLY_SCRIPT=./scripts/fSSU_assembly_v1.2.sh
ADP1=CTGAGCCAKRATCRAACYCT
ADP2=TGYACWCACCGCCCGTC
RUNDIR=`pwd`
NCPU=100

while getopts ':hz1:2:b:a:d:c:t' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    1) RAW_R1=$OPTARG;;
    2) RAW_R2=$OPTARG;;
    b) BARCODES=$OPTARG;;
    a) CLCFOLDER=$OPTARG;;
    d) ALIGNMENT_16S=$OPTARG;;
    c) NCPU=$OPTARG;;
    t) RAW_R1=testdata/mock_R1.fq;RAW_R2=testdata/mock_R2.fq;BARCODES=testdata/barcodes.txt;CLCFOLDER=CLC_assembly_cell;ALIGNMENT_16S=testdata/mock16S.alignment;; 
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

MISSING="is unset but required. Exiting.
"
if [ -z ${RAW_R1+x} ]; then echo "-f $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${RAW_R2+x} ]; then echo "-r $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${BARCODES+x} ]; then echo "-b $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${CLCFOLDER+x} ]; then echo "-a $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${ALIGNMENT_16S+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi;

##################
#### Pipeline ####
##################

echo ""
echo "#########################"
echo "### Extract read-bins ###"
echo "#########################"
perl $EXTRACT_SCRIPT -f $RAW_R1 -r $RAW_R2 -q -b $BARCODES -m 1000 -a 2
mv A1_Read_Stats.txt A2_Read_Stats.txt Anchor_Stats.txt stats.txt ./AnchorReads/

echo ""
echo "####################################"
echo "### De novo assembly of readbins ###"
echo "####################################"
mkdir assembly
SAMPLES=`find ./AnchorReads/ -maxdepth 1 -mindepth 1 -type d`
for EACH_SAMPLE in $SAMPLES
do
  NAME=`basename $EACH_SAMPLE`
  echo "##########################################"
  echo "### Processing sample $NAME ###"
  echo "##########################################"
  sh $ASSEMBLY_SCRIPT $EACH_SAMPLE $CLCFOLDER assembly $ADP1 $ADP2 $NCPU
  sed -i "s/^>/>$NAME\_/" assembly/${NAME}_all_sequences.fa
done

