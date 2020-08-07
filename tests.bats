#!/usr/bin/env bats
load autotax

@test "Variable set: VERSION" {
  [ -n ${VERSION} ]
}

@test "silva_db database file" {
  #expect variable is set
  [ -n ${silva_db} ]

  #expect non-empty and readable file
  checkFiles $silva_db
}

@test "silva_udb database file" {
  #expect variable is set
  [ -n ${silva_udb} ]

  #expect non-empty and readable file
  checkFiles $silva_udb
}

@test "typestrains_udb database file" {
  #expect variable is set
  [ -n ${typestrains_udb} ]

  #expect non-empty and readable file
  checkFiles $typestrains_udb
}

@test "Variable set: denovo_prefix" {
  [ -n ${denovo_prefix} ]
}

@test "Variable set: MAX_THREADS" {
  [ -n ${MAX_THREADS} ]
}

@test "usearch11 in \$PATH" {
	usearch11=$(which usearch11)
	[ -n ${usearch11} ]
}

@test "sina in \$PATH" {
	sina=$(which sina)
	[ -n ${sina} ]
}

@test "R in \$PATH" {
	R=$(which R)
	[ -n ${R} ]
}

@test "Rscript in \$PATH" {
	Rscript=$(which Rscript)
	[ -n ${Rscript} ]
}

@test "Check installed R packages" {
	run checkRPkgs
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
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

@test "Error if temp/ folder exists" {
	#expect error if no arguments passed to function
	run checkFolder
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect error
	mkdir -p ${test_run_dir}temp
	run checkFolder ${test_run_dir}temp
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]
}

@test "Error if output/ folder exists" {
	#expect error if no arguments passed to function
	run checkFolder
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect error
	mkdir -p ${test_run_dir}output
	run checkFolder ${test_run_dir}output
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]
}

@test "Check input data" {
	#expect error
	unset $DATA
	run checkInputData
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	local DATA="non-empty string"
	run checkInputData
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Orient" {
	#test input file name
	local in=${test_dir}example_data/5k_fSSUs.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/fSSUs_oriented.fa

	#expect error if no arguments passed to function
	run orient
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run orient -i $in -d $silva_udb -o ${test_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Dereplication" {
	#test input file name
	local in=${verified_run_dir}temp/fSSUs_oriented.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/uniques_wsize.fa

	#expect error if no arguments passed to function
	run derep
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run derep -i $in -o ${test_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Denoise" {
	#test input file name
	local in=${verified_run_dir}temp/uniques_wsize.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/preFLASVs.fa

	#expect error if no arguments passed to function
	run denoise
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run denoise -i $in -o ${test_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Find longest and rename" {
	#test input file name
	local in=${verified_run_dir}temp/preFLASVs.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/FLASVs.fa

	#expect error if no arguments passed to function
	run findLongest
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run findLongest -i $in -o ${test_run_dir}$out -t $MAX_THREADS
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step (optional): Add additional FLASVs to DB" {
	#test input file name
	local in=${verified_run_dir}temp/FLASVs.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/FLASVs.fa

	#test database file name
	local db=${verified_run_dir}temp/FLASVs.fa

	#expect error if no arguments passed to function
	run addFLASVs
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#add identical, already generated FLASVs, and expect no new unique/redundant FLASVs
	run addFLASVs -i $in -d $db -o ${test_run_dir}$out -t $MAX_THREADS
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect no difference between input and output as $in+$db are the same file
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#add new unique fSSUs
	local in=${test_dir}example_data/5k_fSSUs_addons.fa
	local out=temp/FLASVs_waddons.fa
	run addFLASVs -i $in -d $db -o ${test_run_dir}$out -t $MAX_THREADS
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}FLASVs_waddons.fa
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Global alignment against SILVA" {
	#test input file name
	local in=${verified_run_dir}temp/FLASVs.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/FLASVs_SILVA_aln.fa

	#test database file name
	local db=$silva_db

	#test log file
	local log=temp/sinaAlign_log.txt

	#expect error if no arguments passed to function
	run sinaAlign
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run sinaAlign -i $in -d $db -o ${test_run_dir}$out -t $MAX_THREADS -l ${test_run_dir}$log
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	#due to multithreading, the output can be sorted differently between runs
	run diff -q <(sort -V ${test_run_dir}$out) <(sort -V ${verified_run_dir}$out)
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Trim and strip alignment" {
  #test input file name
	local in=${verified_run_dir}temp/FLASVs_SILVA_aln.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/FLASVs_SILVA_aln_trimmed.fa
	
	#expect error if no arguments passed to function
	run trimStripAlignment
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run trimStripAlignment -i $in -o ${test_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	#due to multithreading, the output is sorted differently between runs
	run diff -q <(sort -V ${test_run_dir}${out}) <(sort -V ${verified_run_dir}${out})
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Sort FLASVs by ID (i.e. highest coverage)" {
  #test input file name
	local in=${verified_run_dir}temp/FLASVs_SILVA_aln_trimmed.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/FLASVs_SILVA_aln_trimmed_sorted.fa
	
	#expect error if no arguments passed to function
	run sortFLASVs
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run sortFLASVs -i $in -o ${test_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Obtaining the taxonomy of the best hit in the SILVA database" {
  #test input file name
	local in=${verified_run_dir}temp/FLASVs_SILVA_aln_trimmed_sorted.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/tax_SILVA.txt

	#test database file name
	local db=$silva_udb

	#expect error if no arguments passed to function
	run searchTaxDB
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run searchTaxDB -i $in -d $db -o ${test_run_dir}$out -t $MAX_THREADS
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	#due to multithreading, the output can be sorted differently between runs
	run diff -q <(sort -V ${test_run_dir}$out) <(sort -V ${verified_run_dir}$out)
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Obtaining the taxonomy of species (>98.7% id) in the SILVA typestrains database" {
  #test input file name
	local in=${verified_run_dir}temp/FLASVs_SILVA_aln_trimmed_sorted.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/tax_typestrains.txt

	#test database file name
	local db=$typestrains_udb

	#expect error if no arguments passed to function
	run searchTaxDB_typestrain
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run searchTaxDB_typestrain -i $in -d $db -o ${test_run_dir}$out -t $MAX_THREADS
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	#due to multithreading, the output can be sorted differently between runs
	run diff -q <(sort -V ${test_run_dir}$out) <(sort -V ${verified_run_dir}$out)
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Cluster at species level" {
  #test input file name
	local in=${verified_run_dir}temp/FLASVs_SILVA_aln_trimmed_sorted.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/SILVA_FLASV-S.txt

	#test centroids file name
	local centroids=temp/SILVA_FLASV-S_centroids.fa

	#expect error if no arguments passed to function
	run clusterSpecies
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run clusterSpecies -i $in -o ${test_run_dir}$out -c ${test_run_dir}$centroids
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Cluster at genus level" {
  #test input file name
	local in=${verified_run_dir}temp/FLASVs_SILVA_aln_trimmed_sorted.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/SILVA_S-G.txt

	#test centroids file name
	local centroids=temp/SILVA_S-G_centroids.fa

	#expect error if no arguments passed to function
	run clusterGenus
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run clusterGenus -i $in -o ${test_run_dir}$out -c ${test_run_dir}$centroids
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Cluster at family level" {
  #test input file name
	local in=${verified_run_dir}temp/FLASVs_SILVA_aln_trimmed_sorted.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/SILVA_G-F.txt

	#test centroids file name
	local centroids=temp/SILVA_G-F_centroids.fa

	#expect error if no arguments passed to function
	run clusterFamily
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run clusterFamily -i $in -o ${test_run_dir}$out -c ${test_run_dir}$centroids
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Cluster at order level" {
  #test input file name
	local in=${verified_run_dir}temp/FLASVs_SILVA_aln_trimmed_sorted.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/SILVA_F-O.txt

	#test centroids file name
	local centroids=temp/SILVA_F-O_centroids.fa

	#expect error if no arguments passed to function
	run clusterOrder
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run clusterOrder -i $in -o ${test_run_dir}$out -c ${test_run_dir}$centroids
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Cluster at class level" {
  #test input file name
	local in=${verified_run_dir}temp/FLASVs_SILVA_aln_trimmed_sorted.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/SILVA_O-C.txt

	#test centroids file name
	local centroids=temp/SILVA_O-C_centroids.fa

	#expect error if no arguments passed to function
	run clusterClass
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run clusterClass -i $in -o ${test_run_dir}$out -c ${test_run_dir}$centroids
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Cluster at phylum level" {
  #test input file name
	local in=${verified_run_dir}temp/FLASVs_SILVA_aln_trimmed_sorted.fa

	#test output file name, dont use "output" as it is reserved by BATS
	local out=temp/SILVA_C-P.txt

	#test centroids file name
	local centroids=temp/SILVA_C-P_centroids.fa

	#expect error if no arguments passed to function
	run clusterPhylum
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 1 ]

	#expect no error
	run clusterPhylum -i $in -o ${test_run_dir}$out -c ${test_run_dir}$centroids
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]
}

@test "Step: Merge and output taxonomy" {
	local tempfolder=${verified_run_dir}temp
	local outputfolder=${test_run_dir}output
	local denovoprefix=$denovo_prefix

	#expect error if no arguments passed to function
	#run mergeTaxonomy
	#echo $output >&2 #redirect to stderr for debugging
	#[ "$status" -eq 1 ]

	#expect no error
	run mergeTaxonomy -t $tempfolder -o $outputfolder -p $denovoprefix
	echo $output >&2 #redirect to stderr for debugging
	[ "$status" -eq 0 ]

	#expect identical result compared to a previous, verified run
	#run diff -q ${test_run_dir}$out ${verified_run_dir}$out
	#echo $output >&2 #redirect to stderr for debugging
	#[ "$status" -eq 0 ]
}