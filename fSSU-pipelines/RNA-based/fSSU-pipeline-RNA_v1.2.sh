#!/bin/bash

#####################
#### Description ####
#####################

usage="$(basename "$0") [-h] [-1 file -2 file -b file -a path -d alignment -t] -- fSSU pipeline: Extraction of readbins, de-novo assembly and trimming.

where:
    -h  Show this help text
    -z  Flag to run pipeline for DNA amplicon data
    -1  Read1 fastq file
    -2  Read2 fastq file
    -b  Sample barcodes
    -a  Path to CLC Assembly Cell installation
    -d  Path to aligned SSU reference database in fasta format with
        50000 alignment positions.
    -t  Flag to run pipeline on test data

Run pipeline on test data:
./fSSU-pipeline-RNA_v1.2.sh -t
./fSSU-pipeline-RNA_v1.2.sh -1 mock_R1.fq -2 mock_R2.fq -b barcodes.txt -a CLC_assembly_cell -d scripts/mock16S.alignment

Run pipleline on real data:
./fSSU-pipeline-RNA_v1.2 -1 R1.fq -2 R2.fq -b barcodes.txt -a CLC_assembly_cell -d SILVA_132_SSURef_Nr99_tax_silva_full_align_trunc.fasta

Requirements - if your programs are named differently, then change the reference in the respective bash scripts
mothur (version 1.37.6; https://www.mothur.org/wiki/Download_mothur)
parallel (version 20150422; https://www.gnu.org/software/parallel/)
CLC assembly cell (version 5.0.3; https://www.qiagenbioinformatics.com/products/clc-assembly-cell/)
perl (v5.18.2; https://www.perl.org/get.html)
"

###################
#### Arguments ####
###################
EXTRACT_SCRIPT=./scripts/fSSU.extract.N15_v5.1.pl
ASSEMBLY_SCRIPT=./scripts/fSSU.CLC.assembly_v2.0.sh
ADP1=TTTTTTTTTT
ADP2=TTTCTGTTGGTGCTGATATTGCCC
RUNDIR=`pwd`

while getopts ':hz1:2:b:a:d:t' option; do
  case $option in
    h) echo "$usage"; exit 1;;
    z) EXTRACT_SCRIPT=./scripts/fSSU.DNA.extract_v1.0.pl;ADP1=TACGGYTACCTTGTTACGACTT;ADP2=AGAGTTTGATCMTGGCTCAG;;
    1) RAW_R1=$OPTARG;;
    2) RAW_R2=$OPTARG;;
    b) BARCODES=$OPTARG;;
    a) CLCFOLDER=$OPTARG;;
    d) ALIGNMENT_16S=$OPTARG;;
    t) RAW_R1=mock_R1.fq;RAW_R2=mock_R2.fq;BARCODES=barcodes.txt;CLCFOLDER=CLC_assembly_cell;ALIGNMENT_16S=scripts/mock16S.alignment;; 
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

missing="is unset but required. Exiting.
"
if [ -z ${RAW_R1+x} ]; then echo "-f $missing"; echo "$usage"; exit 1; fi; 
if [ -z ${RAW_R2+x} ]; then echo "-r $missing"; echo "$usage"; exit 1; fi; 
if [ -z ${BARCODES+x} ]; then echo "-b $missing"; echo "$usage"; exit 1; fi; 
if [ -z ${CLCFOLDER+x} ]; then echo "-a $missing"; echo "$usage"; exit 1; fi; 
if [ -z ${ALIGNMENT_16S+x} ]; then echo "-d $missing"; echo "$usage"; exit 1; fi;

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
  sh $ASSEMBLY_SCRIPT $EACH_SAMPLE $CLCFOLDER assembly $ADP1 $ADP2
  NAME=`basename $EACH_SAMPLE`
  sed -i "s/^>/>$NAME\_/" assembly/${NAME}_all_sequences.fa
done

echo ""
echo "###############################################"
echo "### Extraction of full length SSU sequences ###"
echo "###############################################"
mkdir postprocessing
cat assembly/*_all_sequences.fa > ./postprocessing/all_samples_raw.fa
$CLCFOLDER/clc_adapter_trim -r ./postprocessing/all_samples_raw.fa -s ./postprocessing/all_samples_trim1200.fa -a $ADP1 -a $ADP2 -c 7 -m 1200
sed -i "s/-/_/g" ./postprocessing/all_samples_trim1200.fa
mothur "#align.seqs(fasta=./postprocessing/all_samples_trim1200.fa, reference=$ALIGNMENT_16S, flip=t, threshold=0.9999, processors=24, outputdir=./postprocessing/); screen.seqs(start=1007, end=43230, processors=24, outputdir=./postprocessing/)"
awk '!/^>/ {$0=substr($0, 1007, 42224)}1' ./postprocessing/all_samples_trim1200.good.align > ./postprocessing/all_samples.cut.align
cat ./postprocessing/all_samples.cut.align | tr -d "-" | tr -d "." > SSU.full.length.fa

#######################
#### Testing tools ####
#######################
#cut -f1 ./postprocessing/all_samples_trim1200.bad.accnos | seqtk subseq ./postprocessing/all_samples_trim1200.fa - > test.fa # Extract discarded sequences.

