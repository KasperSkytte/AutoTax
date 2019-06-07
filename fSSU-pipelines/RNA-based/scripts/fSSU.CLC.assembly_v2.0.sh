#!/bin/bash

# fSSU assembly using CLC assembly cell

# Arguments
DATAFOLDER=${1:-data}
CLCFOLDER=${2:-clcfolder}
ASSEMBLYWD=${3:-assembly}
ADP1=${4:-TTTTTTTTTT}
ADP2=${5:-TTTCTGTTGGTGCTGATATTGCCC}

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
echo "### Quality trimming ###"
find  $DATAFOLDER/ -name '*.fastq' |\
parallel --progress -j80% '$CLCFOLDER/clc_quality_trim -r {} -o $QUAL/{/.}_qual.fastq -c 20 -m 100' >/dev/null

echo ""
echo "### Adapter trimming ###"
echo $ADP1 $ADP2
find  $QUAL/ -name '*.fastq' |\
parallel --progress -j80% '$CLCFOLDER/clc_adapter_trim -r {} -s $ADAP/{/.}_adap.fastq -a $ADP1 -a $ADP2 -c 7' >/dev/null

echo ""
echo "### de novo Assembly ###"
find  $ADAP/ -name '*.fastq' |\
parallel --progress -j80% 'assembly_name=`basename {} _qual_adap.fastq`; $CLCFOLDER/clc_assembler -o $OUTDIR/${assembly_name}_out.fa -p no -q {} -m 1200 --cpus 1' >/dev/null

echo ""
echo "### Contig Renaming ###"
find  $OUTDIR -name '*_out.fa' |\
parallel --progress -j80% 'assembly_name=`basename {} _out.fa`; sed -i "s/>/>$assembly_name\_/g" {}'

echo ""
echo "### Create sequence list ###"
find  $OUTDIR/ -name '*_out.fa' | parallel --progress -j80% 'cat {} >> ${OUTDIR}_all_sequences.fa'

#echo ""
echo "### Removing temp files ###"

rm -r $QUAL $ADAP
