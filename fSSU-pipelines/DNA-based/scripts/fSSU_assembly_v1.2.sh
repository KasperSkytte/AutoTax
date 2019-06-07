#!/bin/bash

# fSSU assembly using CLC assembly cell

# Arguments
DATAFOLDER=${1:-data}
CLCFOLDER=${2:-clcfolder}
ASSEMBLYWD=${3:-assembly}
ADP1=$4
ADP2=$5
NCPU=${6:-2}

# Preparation
mkdir -p $ASSEMBLYWD
ASSEMBLYFOLDER=`basename $DATAFOLDER`
mkdir $ASSEMBLYWD/QUAL
QUAL=$ASSEMBLYWD/QUAL
mkdir $ASSEMBLYWD/ADAP
ADAP=$ASSEMBLYWD/ADAP
mkdir $ASSEMBLYWD/$ASSEMBLYFOLDER
OUTDIR=$ASSEMBLYWD/$ASSEMBLYFOLDER

export OUTDIR
export ADAP
export QUAL
export CLCFOLDER
export ADP1
export ADP2

echo ""
echo "### Adapter trimming ###"
echo $ADP1 $ADP2
find  $DATAFOLDER/ -name '*.fastq' |\
parallel --nice 10 --progress -j$NCPU '$CLCFOLDER/clc_adapter_trim -r {} -s $ADAP/{/.}_adap.fastq -a $ADP1 -a $ADP2 -c 7' >/dev/null

echo ""
echo "### Quality trimming ###"
find  $ADAP/ -name '*.fastq' |\
parallel --nice 10 --progress -j$NCPU '$CLCFOLDER/clc_quality_trim -r {} -o $QUAL/{/.}_qual.fastq -c 20 -m 100' >/dev/null


echo ""
echo "### de novo Assembly ###"
find  $QUAL/ -name '*.fastq' |\
parallel --nice 10 --progress -j$NCPU --rpl '{assemblyname} s:.*/::; s/_adap_qual.fastq$//' ' $CLCFOLDER/clc_assembler -o $OUTDIR/{assemblyname}_out.fa -p no -q {} -m 1200 --cpus 1' >/dev/null

echo ""
echo "### Contig Renaming ###"
find  $OUTDIR -name '*_out.fa' |\
parallel --nice 10 --progress -j$NCPU --rpl '{assemblyname} s:.*/::; s/_out.fa$//' 'sed -i "s/>/>{assemblyname}\_/g" {}'

echo ""
echo "### Create sequence list ###"
find  $OUTDIR/ -name '*_out.fa' | parallel --nice 10 --progress -j$NCPU 'cat' > ${OUTDIR}_all_sequences.fa

#echo ""
echo "### Removing temp files ###"

rm -r $QUAL $ADAP
