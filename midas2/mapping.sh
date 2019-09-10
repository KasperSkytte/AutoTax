#!/bin/bash
midas2="MiDAS2_subset_stripped.fa"
midas3="../../../runs/MiDAS3_beta3_20190904/output/ESVs_w_sintax.fa"
usearch11 -fasta_stripgaps MiDAS2_subset_full-length.fasta -fastaout $midas2
usearch11 -usearch_global $midas3 -threads 36 \
  -db $midas2 -id 0 -maxrejects 0 -maxaccepts 0 \
  -top_hit_only -strand plus -blast6out 3v2_ESV.b6
