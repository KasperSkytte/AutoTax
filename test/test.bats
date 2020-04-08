#!/usr/bin/env bats
load /autotax/autotax.sh
mkdir -p /autotax/temp
mkdir -p /autotax/output
export verified_run_dir=/autotax/test/verified_run/ #WITH / AT THE END!

@test "Shell is BASH" {
	#expect error if any arguments are passed to function
	run checkBASH test
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run checkBASH
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Variable set: VERSION" {
  [ ${VERSION} ]
}

@test "Variable set: silva_db" {
  [ ${silva_db} ]
}

@test "Variable set: silva_udb" {
  [ ${silva_udb} ]
}

@test "Variable set: typestrains_db" {
  [ ${typestrains_db} ]
}

@test "Variable set: typestrains_udb" {
  [ ${typestrains_udb} ]
}

@test "Variable set: denovo_prefix" {
  [ ${denovo_prefix} ]
}

@test "Variable set: MAX_THREADS" {
  [ ${MAX_THREADS} ]
}

@test "usearch11 in \$PATH" {
	usearch11=$(which usearch11)
	[ ${usearch11} ]
}

@test "sina in \$PATH" {
	sina=$(which sina)
	[ ${sina} ]
}

@test "R in \$PATH" {
	R=$(which R)
	[ ${R} ]
}

@test "Rscript in \$PATH" {
	Rscript=$(which Rscript)
	[ ${Rscript} ]
}

@test "Check installed R packages" {
	skip "WIP"
}

@test "Echo with timestamp" {
	#expect error if no arguments passed to function
	run echoWithHeader
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run echoWithHeader test
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#check if output pattern matches the format "[2020-03-16 10:08:30]: test"
	pattern="^\[[0-9\ -:]*\]: test$"
	[[ ${lines[0]} =~ $pattern ]]
}

@test "Error if temp/ folder exists" {
	#expect error if no arguments passed to function
	run checkFolder
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect error
	mkdir -p temp
	run checkFolder temp
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]
}

@test "Error if output/ folder exists" {
	#expect error if no arguments passed to function
	run checkFolder
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect error
	mkdir -p output
	run checkFolder output
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]
}

@test "Check input data" {
	skip "WIP"
}

@test "Step: Orient" {
	#test input file name
	local in=/autotax/example_data/10k_fSSUs.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/fSSUs_oriented.fa

	#expect error if no arguments passed to function
	run orient
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run orient -i $in -d $silva_udb -o $out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q $out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Dereplication" {
	#test input file name
	local in=temp/fSSUs_oriented.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/uniques_wsize.fa

	#expect error if no arguments passed to function
	run derep
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run derep -i ${verified_run_dir}$in -o $out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q $out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Denoise" {
	#test input file name
	local in=temp/uniques_wsize.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/preESVs.fa

	#expect error if no arguments passed to function
	run denoise
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run denoise -i ${verified_run_dir}$in -o $out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q $out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Find longest and rename" {
	#test input file name
	local in=temp/preESVs.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/ESVs.fa

	#expect error if no arguments passed to function
	run findLongest
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run findLongest -i ${verified_run_dir}$in -o $out -t $MAX_THREADS
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q $out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step (optional): Add additional ESVs to DB" {
	#test input file name
	local in=temp/ESVs.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/ESVs.fa

	#test database file name
	local db=temp/ESVs.fa

	#expect error if no arguments passed to function
	run addESVs
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#add identical, already generated ESVs, and expect no new unique/redundant ESVs
	run addESVs -i $in -d ${verified_run_dir}$db -o $out -t $MAX_THREADS
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect no difference between input and output as $in+$db are the same file
	run diff -q $out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#add new unique fSSUs
	local in=/autotax/example_data/100_addonESVs.fa
	local out=temp/addESVs.fa
	run addESVs -i $in -d ${verified_run_dir}$db -o $out -t $MAX_THREADS
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q $out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Global alignment against SILVA" {
	#test input file name
	local in=temp/ESVs.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/ESVs_SILVA_aln.fa

	#test database file name
	local db=$silva_db

	#test log file
	local log=temp/sinaAlign_log.txt

	#expect error if no arguments passed to function
	run sinaAlign
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run sinaAlign -i $in -d $db -o $out -t $MAX_THREADS -l $log
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	#run diff -q $out ${verified_run_dir}$out
	#echo $output >&2 #redirect to stderr for debugging
	#[ "$status" -eq 0 ]
}

@test "Step: trim" {
	trimStripAlignment -i temp/ESVs_SILVA_aln.fa -o temp/ESVs_SILVA_aln_trimmed.fa
}

@test "Step: sort" {
	sortESVs -i temp/ESVs_SILVA_aln_trimmed.fa -o temp/ESVs_SILVA_aln_trimmed_sorted.fa
}

@test "Step: Obtaining the taxonomy of the best hit in the SILVA database" {
	searchTaxDB -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -d $silva_udb -o temp/tax_SILVA.txt -t $MAX_THREADS
}

@test "Step: Obtaining the taxonomy of species (>98.7% id) in the SILVA typestrains database" {
	searchTaxDB_typestrain -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -d $typestrains_udb -o temp/tax_typestrains.txt -t $MAX_THREADS
}

@test "Step: Cluster at species level (98.7%)" {
	clusterSpecies -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_ESV-S.txt -c temp/SILVA_ESV-S_centroids.fa
}

@test "Step: Cluster at genus level (94.5%)" {
	clusterGenus -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_S-G.txt -c temp/SILVA_S-G_centroids.fa
}

@test "Step: Cluster at family level (86.5%)" {
	clusterFamily -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_G-F.txt -c temp/SILVA_G-F_centroids.fa
}

@test "Step: Cluster at order level (82.0%)" {
	clusterOrder -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_F-O.txt -c temp/SILVA_F-O_centroids.fa
}

@test "Step: Cluster at class level (78.5%)" {
	clusterClass -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_O-C.txt -c temp/SILVA_O-C_centroids.fa
}

@test "Step: Cluster at phylum level (75.0%)" {
	clusterPhylum -i temp/ESVs_SILVA_aln_trimmed_sorted.fa -o temp/SILVA_C-P.txt -c temp/SILVA_C-P_centroids.fa
}

@test "Rstuff" {
  Rstuff
  skip
}