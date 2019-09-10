# MiDAS3.3
fSSU's were first run with the script as-is and 9521 ESVs were generated.
Then the addon sequences were used as input with the 9521 ESVs generated first as database (-d), and the denoising step was skipped (just commented out in the script). 

Afterwards the addon sequences were mapped back to the resulting ESVs_w_sintax.fa and the replacements file was filled out manually based on the mapping. The mapping was done with `usearch -usearch_global` with the `-id 0 -strand both -maxaccepts 0 -maxrejects 0 -top_hit_only` settings. Addon sequences that got classified with a SILVA typestrain were skipped, while those that got a denovo name was filled at all levels according to the article (PMID) for the particular sequence. 

Finally a subset of the MiDAS2 database containing only those manually added to MiDAS2 was mapped to the ESVs and entered accordingly in the replacements file using the scripts in the midas2 folder.
