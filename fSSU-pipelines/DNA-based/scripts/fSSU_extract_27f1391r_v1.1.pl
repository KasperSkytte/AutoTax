###############################################################################
#
#    fSSU_DNA_extract_v1.1.pl
#
#    
#    Copyright (C) 2018 Mads Albertsen & Søren Karst
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;
use Cwd;

#core Perl modules
use Getopt::Long;

#locally-written modules
BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

my $version = "4.1";

# get input params
my $global_options = checkParams();

my $R1;
my $R2;
my $stats;
my $afq;
my $sfq;
my $fastq;
my $bc;
my $max_seq;
my $print_map;
my $print_raw;
my $only_stats;

$R1 = &overrideDefault("R1.fq",'R1');
$R2 = &overrideDefault("R2.fq",'R2');
$stats = &overrideDefault("stats.txt",'stats');
$afq = &overrideDefault(5,'afq');
$sfq = &overrideDefault(25,'sfq');
$fastq = &overrideDefault(0,'fastq');
$bc = &overrideDefault("none",'barcodes');
$max_seq = &overrideDefault(1000,'maxseq');
$print_map = &overrideDefault(0,'printmap');
$print_raw = &overrideDefault(0,'printraw');
$only_stats = &overrideDefault(0,'onlystats');

# Barcode variables
my %h_a1_barcodes;                                                   # Barcodes of sub-anchor 1
my %h_a2_barcodes;                                                   # Barcodes of sub-anchor 2
my %h_barcodes;                                                      # Combined barcodes

# Global - variables
my $line_nr = 0;
my $sequence = "";
my $read_index = "";
my $out_file = "error";
my $print_map_header = "error";
my $print_raw_header = "error";
my $ext = "fasta";
my $header = "";
my $print_read = 0;
my $raw_print = 0;
my $progress_count = 0;
my $progress_total = 0;

# Global - hash
my %h_a1_index;                                                      # Counts and saves the abundance of each a1 index
my %h_a2_index;                                                      # Counts and saves the abundance of each a2 index
my %h_a1_read;                                                       # Counts and saves the abundance of each a1 read
my %h_a2_read;                                                       # Counts and saves the abundance of each a2 read
my %h_a1_keep;                                                       # Marks good a1 reads for extraction (links: a1 -> index)
my %h_a2_keep;                                                       # Marks good a2 reads for extraction (links: a2 -> index)
my %h_index;                                                         # Counts and saves the complete anchor index
my %h_barcode_index;                                                 # Save the barcode of the anchor (links: barcode -> index)
my %h_barcode_index_error;                                           # Flag anchors with similar index, but different barcodes, for removal
my %h_index_occurance;                                               # Store the occurance of each sub-anchor. Sorted by global index abundance
my %h_index_keep;                                                    # Anchors passing filtering 
my %h_read_out;                                                      # Count the number of passed reads exported from each index
my %h_raw_out;                                                       # Count the number of raw reads exported from each index

## Stats - hash
my %h_stats_barcode_index_total_potential;                           # Total number of, barcode specific, potential total anchor reads
my %h_stats_barcode_index_unique_potential;                          # Total number of, barcode specific, potential unique anchor reads
my %h_stats_barcode_index_unique_keep;                               # Total number of, barcode specific, anchor reads kept

## Stats - variables
my $stats_reads_total = 0;                                           # Total number of reads
my $stats_anchor_total_potential = 0;                                # Total number of potential anchor reads
my $stats_anchor_total_barcode = 0;                                  # Total number of anchor reads with barcode
my $stats_anchor_unique_keep = 0;                                    # Count the number of good indexes
my $stats_anchor_total_keep = 0;                                     # Count the anchor reads associated with good indexes

my $stats_filter_index_unique_abundance = 0;                         # Number of index removed by index abundance
my $stats_filter_index_unique_occurance = 0;                         # Number of index removed by index occurance abundance
my $stats_filter_a_read_unique_abundance = 0;                        # Number of index removed, where both sub-index have less than X reads associated
my $stats_filter_a1_read_unique_abundance = 0;                       # Number of index removed, where a1-index have less than X reads associated
my $stats_filter_a2_read_unique_abundance = 0;                       # Number of index removed, where a2-index have less than X reads associated
my $stats_filter_index_total_abundance = 0;                          # Number of index reads removed by index abundance
my $stats_filter_index_total_occurance = 0;                          # Number of index reads removed by index occurance abundance
my $stats_filter_a_read_total_abundance = 0;                         # Number of index reads removed, where both sub-index have less than X reads associated
my $stats_filter_a1_read_total_abundance = 0;                        # Number of index reads removed, where a1-index have less than X reads associated
my $stats_filter_a2_read_total_abundance = 0;                        # Number of index reads removed, where a2-index have less than X reads associated


my $stats_a1_reads_total_potential = 0;                              # Total number of potential sub-anchor 1 reads
my $stats_a1_reads_total_barcode = 0;                                # Total number of sub-anchor 1 reads with barcode
my $stats_a1_reads_unique_filter = 0;                                # Unique sub-anchor 1 reads after abundance filtering
my $stats_a1_reads_total_filter = 0;                                 # Total number of sub-anchor 1 reads after abundance filtering
my $stats_a1_reads_unique_in_a1_index = 0;                           # Unique sub-anchor 1 reads after abundance filtering found in a1 index (in anchor)
my $stats_a1_reads_total_in_a1_index = 0;                            # Total number of sub-anchor 1 reads after abundance filtering found in a1 index (in anchor)
my $stats_a1_reads_keep_total = 0;                                   # Total number of sub-anchor 1 reads in passed anchors

my $stats_a2_reads_total_potential = 0;                              # Total number of potential sub-anchor 2 reads
my $stats_a2_reads_total_barcode = 0;                                # Total number of sub-anchor 2 reads with barcode
my $stats_a2_reads_unique_filter = 0;                                # Unique sub-anchor 2 reads after abundance filtering
my $stats_a2_reads_total_filter = 0;                                 # Total number of sub-anchor 2 reads after abundance filtering
my $stats_a2_reads_unique_in_a2_index = 0;                           # Unique sub-anchor 2 reads after abundance filtering found in a2 index (in anchor)
my $stats_a2_reads_total_in_a2_index = 0;                            # Total number of sub-anchor 2 reads after abundance filtering found in a2 index (in anchor)
my $stats_a2_reads_keep_total = 0;                                   # Total number of sub-anchor 2 reads in passed anchors

######################################################################
# CODE HERE
######################################################################

### Importing barcodes if supplied ###################################

if($bc ne "none"){
	open(BC_fh, $bc) or die("Cannot read file: $bc\n");             # Load barcodes from file. Assumes that the forward and reverse barcode is joined e.g. AAAAAAAAGGGGGGGG
	while (my $line = <BC_fh>){
		chomp $line;
		$h_barcodes{$line} = 0;
		my $a1 = substr($line, 0, 8);
		my $a2 = substr($line, 8, 8);
		$h_a1_barcodes{$a1} = 0;
		$h_a2_barcodes{$a2} = 0;
	}
	close BC_fh;
} else{
###	$h_barcodes{"TAAGGCGAATAGAGGC"} = 0; $h_a1_barcodes{"TAAGGCGA"} = 0; $h_a2_barcodes{"ATAGAGGC"} = 0; ### Load a set of standard barcodes
	$h_barcodes{"CTAGTACGTGCCTCTT"} = 0; $h_a1_barcodes{"CTAGTACG"} = 0; $h_a2_barcodes{"TGCCTCTT"} = 0; ### Load a set of standard barcodes
}

# Make an output folder for each barcode
my $dir = cwd();
my $sub_folder = "$dir/AnchorReads";
if (-d $sub_folder){
	die("Aborting: The output folder already exists!\n");
	} else{
	mkdir $sub_folder;
	foreach my $barcode (keys %h_barcodes) {
		mkdir "$sub_folder/$barcode";
		}
	}

### Identifying Anchors in Read 1 ###################################

open(INR1, $R1) or die("Cannot read file: $R1\n");       

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $nice_timestamp = sprintf ( "%02d:%02d:%02d", $hour,$min,$sec);
print "\n$nice_timestamp Identifying Anchors\n";       

while (my $line1 = <INR1>){                                                                        # The anchors are in R1
	chomp $line1;
	$line_nr++;
	if ($line_nr == 2){
		$stats_reads_total++; 
###		if ($line1 =~ m/^GAAGA/) {$stats_anchor_total_potential++;}                              # Count the number of sequences from potential anchors - e.g. they start with the correct 5 bases
		if ($line1 =~ m/^TG.AC.CA/) {$stats_anchor_total_potential++;}                              # Count the number of sequences from potential anchors - e.g. they start with the correct 5 bases
	     # Identify sub-anchor 1 index
		my $potential_a1 = substr($line1, 15, 35);		                                       # Extract position -3 and +10 of theoretical A1 barcode + Index (18, 40), this is done due to size variation in the index
		my $a1_barcode = "None";
		my $a1_index = 0;
		foreach my $key (keys %h_a1_barcodes){                  
			if ($potential_a1 =~ m/$key/) {                                                     # The barcode should be in the target area
				my @sl = split($key, $potential_a1, 2);		                                  # Split the sequence to retrieve the index
				if ($sl[1] =~ m/ATCTTC/) {	
					my @sl2 = split(/ATCTTC/, $sl[1], 2);
					my $a1_length = length($sl2[0]);		
					next if($a1_length < 14 or $a1_length > 16);                               # Reject if the index length is outside pre-defined criterias			
					$a1_index = $sl2[0];                                                      # The sub-anchor needs to be reverse complemented for extraction of reads. This is done up-front here.
			          $a1_index =~ tr/ACGTacgt/TGCAtgca/;                        
					$a1_index = reverse($a1_index);
					$a1_barcode = $key;                                                       # Store the barcode of the index - used to check if both barcodes were correct
 				}
			}        
		} 

	     # Identify sub-anchor 2 index
		my $potential_a2 = substr($line1, 73, 43);		     	                             # Extract position -10 and +10 of theoretical A1 barcode + Index (83, 23)
		my $a2_barcode = "None";
		my $a2_index = 0;
		foreach my $key (keys %h_a2_barcodes){
			if ($potential_a2 =~ m/$key/) {
				my @sl = split($key, $potential_a2, 2);				
				if ($sl[0] =~ m/ATCCAT/) {	
					my @sl2 = split(/ATCCAT/, $sl[0], 2);
					my $a2_length = length($sl2[1]);					
					next if($a2_length < 14 or $a2_length > 16);
					$a2_index = $sl2[1];
					$a2_barcode = $key;
 				}
							
			}        
		} 

		my $a_barcode = $a1_barcode.$a2_barcode;	
		next if(!exists($h_barcodes{$a_barcode}));                                               # Next if incorrect combination of barcodes
		$stats_anchor_total_barcode++;                                                           # Count anchors with correct barcode   
	
		# Test for duplication of the index
		next if($a1_index eq $a2_index);                                                         ### NOTE: Could also test if the "normal", e.g. non-reverse complemented anchors are the same? Could that indicate an error?

		# Save the potential good sub anchor index
		$h_a1_index{$a1_index}++;                                                                # Save and count the abundance of sub-anchor 1
		$h_a2_index{$a2_index}++;                                                                # Save and count the abundance of sub-anchor 2

		# Save the potential good anchor index
		my $a_index = $a1_index."-".$a2_index;
		$h_index{$a_index}++;                          	                                       # Save and count the abundance of the Anchor index
		$h_barcode_index{$a_index} = $a_barcode;                                                 # Save the barcode of the anchor

		# Test if multiple barcodes have the same anchor
		if(exists($h_barcode_index{$a_index})){                                                  # Test for barcode conflict - same anchor different barcodes
			if($a_barcode ne $h_barcode_index{$a_index}){                                       # NOTE: Flag these for removal <- Maybe it should just take the most abundant one - should be rare occurance. But, seems not like a problem in the data anyway..?
				$h_barcode_index_error{$a_index} = 1;
			}
		}
		
		# Count index occurance pr. barcode
		$h_stats_barcode_index_total_potential{$a_barcode}++;                                    # Total barcode specific anchors - NOTE: these could be combined to one statistic and calculated later, although not more effective code-wise
		$h_stats_barcode_index_unique_potential{$a_barcode}{$a_index} = 1;                       # Unique barcode specific anchors
	}
     if ($line_nr == 4){$line_nr = 0};
}

close INR1;

### Identifying sub-Anchors in reads - located in read 2 ############

open(INR2, $R2) or die("Cannot read file: $R2\n");                   

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$nice_timestamp = sprintf ( "%02d:%02d:%02d", $hour,$min,$sec);
print "$nice_timestamp Identifying sub-Anchors in reads\n";

$line_nr = 0;

while (my $line2 = <INR2>){                                  
	chomp $line2;    
	$line_nr++;
	if ($line_nr == 2){       		
### Probably not working?
		if ($line2 =~ m/^AGGGG/) {$stats_a1_reads_total_potential++;}                            # Count the total potential a1 associated reads
		if ($line2 =~ m/^CTCCA/) {$stats_a2_reads_total_potential++;}                            # Count the total potential a2 associated reads

          # Identify sub-anchor 1 reads
		foreach my $a1_barcode (keys %h_a1_barcodes){
			$a1_barcode =~ tr/ACGTacgt/TGCAtgca/;  
			$a1_barcode = reverse($a1_barcode); 
			if ($line2 =~ m/$a1_barcode/) {
				my @sl = split($a1_barcode, $line2, 2);				
					my $a1_length = length($sl[0]);					
					next if($a1_length < 14 or $a1_length > 16);
					my $a1_index = substr($sl[0], 0, $a1_length);
				    	$h_a1_read{$a1_index}++;
					$stats_a1_reads_total_barcode++;  
			}
		} 

          # Identify sub-anchor 1 reads
		foreach my $a2_barcode (keys %h_a2_barcodes){
			if ($line2 =~ m/$a2_barcode/) {
				my @sl = split($a2_barcode, $line2, 2);				
					my $a2_length = length($sl[0]);					
					next if($a2_length < 14 or $a2_length > 16);
					my $a2_index = substr($sl[0], 0, $a2_length);
				    	$h_a2_read{$a2_index}++;
					$stats_a2_reads_total_barcode++;	
			}
		} 
	}
     if ($line_nr == 4){$line_nr = 0};
}

close INR2;

### Filter Anchors ###################################################

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$nice_timestamp = sprintf ( "%02d:%02d:%02d", $hour,$min,$sec);
print "$nice_timestamp Filtering Anchors\n";                         

foreach my $index (sort { $h_index{$b} <=> $h_index{$a}} keys %h_index) {                          # Sort the anchor index by abundance
	my @sl =	split("-",$index);
	my $a1_index = $sl[0];                                                                        # Retrive each sub-anchor
	my $a2_index = $sl[1];    
     $h_index_occurance{$a1_index}++;                                                              # Count the occurance of sub-anchor 1
	$h_index_occurance{$a2_index}++;                                                              # Count the occurance of sub-anchor 2                       
	next if ($h_index{$index} < $afq);                                                            # Discard if abundance of the specific index is less than X
	next if (exists($h_barcode_index_error{$index}));                                             # Discard if two different barcodes were assigned to the same anchor
	if(!exists($h_a1_read{$a1_index})){ $h_a1_read{$a1_index} = 0;}
	if(!exists($h_a2_read{$a2_index})){ $h_a2_read{$a2_index} = 0;}
	next if ($h_a1_read{$a1_index} < $sfq);                                                       # Discard if the number of reads associated to sub-anchor 1 index is less than X.
	next if ($h_a2_read{$a2_index} < $sfq);                                                       # Discard if the number of reads associated to sub-anchor 2 index is less than X.
	if($h_index_occurance{$a1_index} == 1 && $h_index_occurance{$a2_index} == 1){                 # Test if any of the sub-anchors have been seen before				
		$h_a1_keep{$a1_index} = $index;                                                          # Mark good a1 reads for extraction (links: a1 -> index)
		$h_a2_keep{$a2_index} = $index;
		$stats_anchor_unique_keep++;                                                             # Count the number of good indexes
		$stats_anchor_total_keep += $h_index{$index};                                            # Count the anchor reads associated with good indexes
		$h_index_keep{$index} = $h_index{$index};                                                # Save the good indexes for faster lookup later - NOTE: is this used somewhere?
		$h_stats_barcode_index_unique_keep{$h_barcode_index{$index}}++;                          # Count the number of occurances of each barcode in good anchors
	}
}

### Anchor Stats #####################################################
open(aSTATS, ">Anchor_Stats.txt") or die("Cannot create file: Anchor_Stats.txt\n");
print aSTATS "Anchor\ta_count\tKeep\ta1_count\ta2_count\ta1_seq\ta2_seq\ta1_or\ta2_or\ta1_oc\ta2_oc\n";

my %h_index_occurance_temp;

foreach my $index (sort { $h_index{$b} <=> $h_index{$a}} keys %h_index) {                          # Sort the anchors by abundance
	my @sl =	split("-",$index);
	my $a1_index = $sl[0];                                                                        # Retrive each sub-anchor
	my $a2_index = $sl[1];      
	$h_index_occurance_temp{$a1_index}++;
	$h_index_occurance_temp{$a2_index}++;
     my $keep = "No";
	if($h_index_occurance_temp{$a1_index} == 1 && $h_index_occurance_temp{$a2_index} == 1){$keep = "Yes";}        # Test if any of the sub-anchors have been seen before
	if(!exists($h_a1_read{$a1_index})){ $h_a1_read{$a1_index} = 0;}
	if(!exists($h_a2_read{$a2_index})){ $h_a2_read{$a2_index} = 0;}
	if ($h_index{$index} < $afq){$keep = "No"};                                                   # Discard if abundance is less than X
	if (exists($h_barcode_index_error{$index})){$keep = "No"};                                    # Discard if two different barcodes w
	if ($h_a1_read{$a1_index} < $sfq){$keep = "No"};                                              # Discard if the number of reads associated to anchor 1 is less than X.
	if ($h_a2_read{$a2_index} < $sfq){$keep = "No"};                                              # Discard if the number of reads associated to anchor 2 is less than X.
	print aSTATS "$index\t$h_index{$index}\t$keep\t$h_a1_index{$a1_index}\t$h_a2_index{$a2_index}\t$h_a1_read{$a1_index}\t$h_a2_read{$a2_index}\t$h_index_occurance_temp{$a1_index}\t$h_index_occurance_temp{$a2_index}\t$h_index_occurance{$a1_index}\t$h_index_occurance{$a2_index}\n";	
}

close aSTATS;


### Anchor filter Stats #############################################

%h_index_occurance_temp = ();

foreach my $index (sort { $h_index{$b} <=> $h_index{$a}} keys %h_index) {                          # Sort the anchors by abundance
	my @sl =	split("-",$index);
	my $a1_index = $sl[0];                                                                        # Retrive each sub-anchor index
	my $a2_index = $sl[1];     
	$h_index_occurance_temp{$a1_index}++;
	$h_index_occurance_temp{$a2_index}++;
	if ($h_index{$index} < $afq){$stats_filter_index_unique_abundance++; $stats_filter_index_total_abundance += $h_index{$index}; next;};                         # Discard if index abundance is less than X
	if($h_index_occurance_temp{$a1_index} != 1 or $h_index_occurance_temp{$a2_index} != 1){$stats_filter_index_unique_occurance++; $stats_filter_index_total_occurance += $h_index{$index}; next;}                  # Discard if not first occurance 
	if(($h_a1_read{$a1_index} < $sfq) and ($h_a2_read{$a2_index} < $sfq)){$stats_filter_a_read_unique_abundance++; $stats_filter_a_read_total_abundance += $h_index{$index}; next;}                                  # Discard if both less than X
	if ($h_a1_read{$a1_index} < $sfq){$stats_filter_a1_read_unique_abundance++; $stats_filter_a1_read_total_abundance += $h_index{$index}};                         # Discard if the number of reads associated to anchor 1 is less than X.
	if ($h_a2_read{$a2_index} < $sfq){$stats_filter_a2_read_unique_abundance++; $stats_filter_a2_read_total_abundance += $h_index{$index}};                         # Discard if the number of reads associated to anchor 2 is less than X.
}

### Read Stats #####################################################
# A1

open(a1STATS, ">A1_Read_Stats.txt") or die("Cannot create file: A1_Read_Stats.txt\n");
print a1STATS "Sub-anchor\tSeqCount\tKept\tInAnchor\n";

foreach my $a1_index (sort { $h_a1_read{$b} <=> $h_a1_read{$a}} keys %h_a1_read) {                # Sort the read sub-anchors by abundance
	next if($h_a1_read{$a1_index} == 0);
	my $keep = 0;
	if (exists($h_a1_keep{$a1_index})){ $keep = 1;}       				                        # Check if the sub-anchor is part of a good anchor
	my $read_in_anchor = 0;
	if (exists($h_a1_index{$a1_index})){ $read_in_anchor = 1;}                                    # Check if the sub-anchor is seen in any anchor
	if ($h_a1_read{$a1_index} > $sfq-1){
		$stats_a1_reads_unique_filter++;
		$stats_a1_reads_total_filter += $h_a1_read{$a1_index}; 
		if (exists($h_a1_index{$a1_index})){                                                     # Check if the sub-anchor is seen in any anchor 
			$stats_a1_reads_unique_in_a1_index++;
			$stats_a1_reads_total_in_a1_index += $h_a1_read{$a1_index};
			if (exists($h_a1_keep{$a1_index})){
				$stats_a1_reads_keep_total += $h_a1_read{$a1_index};				
			}
		}
		       
	}	                                                
	print a1STATS "$a1_index\t$h_a1_read{$a1_index}\t$keep\t$read_in_anchor\n";
}

close a1STATS;

# A2

open(a2STATS, ">A2_Read_Stats.txt") or die("Cannot create file: A2_Read_Stats.txt\n");
print a2STATS "Sub-anchor\tSeqCount\tKept\tInAnchor\n";

foreach my $a2_index (sort { $h_a2_read{$b} <=> $h_a2_read{$a}} keys %h_a2_read) {                 # Sort the read sub-anchors by abundance
	next if($h_a2_read{$a2_index} == 0);
	my $keep = 0;
	if (exists($h_a2_keep{$a2_index})){ $keep = 1;}       				                        # Check if the sub-anchor is part of a good anchor
	my $read_in_anchor = 0;
	if (exists($h_a2_index{$a2_index})){ $read_in_anchor = 1;}                                    # Check if the sub-anchor is seen in any anchor
	if ($h_a2_read{$a2_index} > $sfq-1){
		$stats_a2_reads_unique_filter++;
		$stats_a2_reads_total_filter += $h_a2_read{$a2_index}; 
		if (exists($h_a2_index{$a2_index})){                                                     # Check if the sub-anchor is seen in any anchor 
			$stats_a2_reads_unique_in_a2_index++;
			$stats_a2_reads_total_in_a2_index += $h_a2_read{$a2_index};
			if (exists($h_a2_keep{$a2_index})){
				$stats_a2_reads_keep_total += $h_a2_read{$a2_index};				
			}
		}       
	}	                                                
	print a2STATS "$a2_index\t$h_a2_read{$a2_index}\t$keep\t$read_in_anchor\n";
}

close a2STATS;

### Extract reads from good Anchors ##################################

if ($only_stats == 0){

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$nice_timestamp = sprintf ( "%02d:%02d:%02d", $hour,$min,$sec);
print "$nice_timestamp Extracting reads and writing to files\n";     

if ($fastq == 1){$ext = "fastq";}

# Open the R1 and R2 files
open(INR1, $R1) or die("Cannot read file: $R1\n");                                                 # Index is located in Read 1
open(INR2, $R2) or die("Cannot read file: $R2\n");                                                 # Associated reads are located in Read 2

if ($print_map == 1){open(OUTr, ">map.$ext") or die("Cannot write file: map.$ext\n");}                 # Make a single combined file with reads from passing index's - stamped with barcode and index 
if ($print_raw == 1){open(OUTx, ">map_raw.$ext") or die("Cannot write file: map_raw.$ext\n");}         # Make a single combined file with ALL reads from sub-anchors - stamped with barcode and index 

$line_nr = 0;

while ( defined(my $line1 = <INR1>) and  defined(my $line2 = <INR2>) ) {                           # The Read 1 (index) and Read 2 (sequence) files are read in parallel.
	chomp $line1;
	chomp $line2; 
	$line_nr++;
	if ($line_nr == 1){$header = $line1}                                                          # Store the header of the read in case it needs to be printed

	# Check for reads for index that passes filtering criterias
	if ($line_nr == 2){
		$sequence = $line1;                                                                      # Store the sequence in case it needs to be printed

		# A1
		foreach my $a1_barcode (keys %h_a1_barcodes){
			$a1_barcode =~ tr/ACGTacgt/TGCAtgca/;  
			$a1_barcode = reverse($a1_barcode); 
			if ($line2 =~ m/$a1_barcode/) {                                         # Test if the barcode is found within the read
				$read_index = "A1";
				my @sl = split($a1_barcode, $line2, 2);				         # Split based on flanking positions to allow variability in the index size
					my $a1_length = length($sl[0]);					
					next if($a1_length < 14 or $a1_length > 16);                               # Remove any index that is outside +/- 1 of the expected index size (15) 
					my $a1_index = substr($sl[0], 0, $a1_length);
					if (exists($h_a1_keep{$a1_index})){ 
						$h_read_out{$a1_index}++;                                              # Count the number of reads for each sub-anchor
						if($h_read_out{$a1_index} <= ($max_seq/2)){                            # Only print X reads in total from each anchor
							$print_read = 1; 
							$out_file = "$sub_folder/$h_barcode_index{$h_a1_keep{$a1_index}}/$h_a1_keep{$a1_index}.$ext";
							$print_map_header = "\@".$h_barcode_index{$h_a1_keep{$a1_index}}."-".$h_a1_keep{$a1_index}."-".$h_read_out{$a1_index}."-A1";
						} 
					}
					if ($print_raw == 1 && exists($h_a1_index{$a1_index})){ 
						if($h_a1_index{$a1_index} >= $sfq){ 
							$h_raw_out{$a1_index}++; 
							if($h_raw_out{$a1_index} <= ($max_seq/2)){
								$raw_print = 1; 
								$print_raw_header = "\@".$a1_barcode."-".$a1_index."-".$h_raw_out{$a1_index}."-A1";
							}
						}
					}
			}
		} 
		# A2
		foreach my $a2_barcode (keys %h_a2_barcodes){
			if ($line2 =~ m/$a2_barcode/) {                                         # Test if the barcode is found within the read
				$read_index = "A2";
				my @sl = split($a2_barcode, $line2, 2);				         # Split based on flanking positions to allow variability in the index size
					my $a2_length = length($sl[0]);					
					next if($a2_length < 14 or $a2_length > 16);                               # Remove any index that is outside +/- 1 of the expected index size (15) 
					my $a2_index = substr($sl[0], 0, $a2_length);					
					if (exists($h_a2_keep{$a2_index})){ 
						$h_read_out{$a2_index}++;                                              # Count the number of reads for each sub-anchor
						if($h_read_out{$a2_index} <= ($max_seq/2)){                            # Only print X reads in total from each anchor
							$print_read = 1; 
							$out_file = "$sub_folder/$h_barcode_index{$h_a2_keep{$a2_index}}/$h_a2_keep{$a2_index}.$ext";
							$print_map_header = "\@".$h_barcode_index{$h_a2_keep{$a2_index}}."-".$h_a2_keep{$a2_index}."-".$h_read_out{$a2_index}."-A2";
						} 
					}
					if ($print_raw == 1 && exists($h_a2_index{$a2_index})){ 
						if($h_a2_index{$a2_index} >= $sfq){ 
							$h_raw_out{$a2_index}++; 
							if($h_raw_out{$a2_index} <= ($max_seq/2)){
								$raw_print = 1; 
								$print_raw_header = "\@".$a2_barcode."-".$a2_index."-".$h_raw_out{$a2_index}."-A2";								
							}
						}
					}
			}
		}
	}
	
	# Write the good sequenes to their respective files
     if ($line_nr == 4){
		if ($print_read == 1){                                                                   # Output the sequence if matches a good sub-anchor
			open(OUT, ">>$out_file") or die("Cannot create file: $out_file\n");
			if ($fastq == 1){ print OUT "$header\n$sequence\n+\n$line1\n";} else {$header =~ s/^\@/\>/; print  OUT "$header\n$sequence\n";}
			close OUT;
			# Optional: Print all good anchor reads to a single file for further analysis
			if ($print_map == 1){ if ($fastq == 1){print OUTr "$print_map_header\n$sequence\n+\n$line1\n";} else {$print_map_header =~ s/^\@/\>/;print OUTr "$print_map_header\n$sequence\n";}}
			$out_file = "error";
			$print_map_header = "error";
		}
		# Optional: Print all sub-anchor reads to a single file for further analysis
		if ($raw_print == 1){ if ($fastq == 1){print OUTx "$print_raw_header\n$sequence\n+\n$line1\n";} else {$print_raw_header =~ s/^\@/\>/; print OUTx "$print_raw_header\n$sequence\n";}}
		$line_nr = 0;
		$print_read = 0;
		$raw_print = 0;
	
	# For progress reporting
		$progress_count++;
		$progress_total++;
		if ($progress_count >= ($stats_reads_total / 1000)){
			$progress_count = 0;
			my $pct = sprintf("%.0f", $progress_total/$stats_reads_total*100);
			($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
			$nice_timestamp = sprintf ( "%02d:%02d:%02d", $hour,$min,$sec);
			print "$nice_timestamp $pct %\n";
		}
	}
}

close INR1;
close INR2;
if ($print_map == 1){close OUTr;}
if ($print_raw == 1){close OUTx;}

}

### Calculate and export stats ########################################

open(STATS, ">$stats") or die("Cannot create file: $stats\n");

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$nice_timestamp = sprintf ( "%02d:%02d:%02d", $hour,$min,$sec);
print "$nice_timestamp Calulating stats and writing them to $stats\n";                          # To get nice time-stamps

my $ttime =localtime(time);
print STATS "$ttime\n\n";

print STATS "Script parameters:\n";
print STATS "Version: $version\nRead 1 file: $R1\nRead 2 file: $R2\nMinimum index frequency: $afq\nMinimum sub-index read frequency: $sfq\nMaximum number of sequences to output: $max_seq\n\n";

my $pct_stats_anchor_total_potential = sprintf("%.0f",$stats_anchor_total_potential/$stats_reads_total*100);
my $pct_stats_anchor_total_barcode = sprintf("%.0f",$stats_anchor_total_barcode/$stats_reads_total*100);
my $stats_anchor_unique = keys %h_index;

my $stats_filter_index_unique_abundance_left = $stats_anchor_unique - $stats_filter_index_unique_abundance;
my $stats_filter_index_unique_occurance_left = $stats_filter_index_unique_abundance_left - $stats_filter_index_unique_occurance;
my $stats_filter_a_read_unique_abundance_left = $stats_filter_index_unique_occurance_left - $stats_filter_a_read_unique_abundance;
my $stats_filter_a1_read_unique_abundance_left = $stats_filter_a_read_unique_abundance_left - $stats_filter_a1_read_unique_abundance;
my $stats_filter_a2_read_unique_abundance_left = $stats_filter_a1_read_unique_abundance_left - $stats_filter_a2_read_unique_abundance;

my $stats_filter_index_total_abundance_left = $stats_anchor_total_barcode - $stats_filter_index_total_abundance;
my $stats_filter_index_total_occurance_left = $stats_filter_index_total_abundance_left - $stats_filter_index_total_occurance;
my $stats_filter_a_read_total_abundance_left = $stats_filter_index_total_occurance_left - $stats_filter_a_read_total_abundance;
my $stats_filter_a1_read_total_abundance_left = $stats_filter_a_read_total_abundance_left - $stats_filter_a1_read_total_abundance;
my $stats_filter_a2_read_total_abundance_left = $stats_filter_a1_read_total_abundance_left - $stats_filter_a2_read_total_abundance;

my $pct_stats_filter_index_unique_occurance_left = sprintf("%.0f",$stats_filter_index_unique_occurance_left/$stats_filter_index_unique_abundance_left*100);
my $pct_stats_filter_a_read_unique_abundance_left = sprintf("%.0f",$stats_filter_a_read_unique_abundance_left/$stats_filter_index_unique_abundance_left*100);
my $pct_stats_filter_a1_read_unique_abundance_left = sprintf("%.0f",$stats_filter_a1_read_unique_abundance_left/$stats_filter_index_unique_abundance_left*100);
my $pct_stats_filter_a2_read_unique_abundance_left = sprintf("%.0f",$stats_filter_a2_read_unique_abundance_left/$stats_filter_index_unique_abundance_left*100);

my $pct_stats_filter_index_total_abundance_left = sprintf("%.0f",$stats_filter_index_total_abundance_left/$stats_reads_total*100);
my $pct_stats_filter_index_total_occurance_left = sprintf("%.0f",$stats_filter_index_total_occurance_left/$stats_reads_total*100);
my $pct_stats_filter_a_read_total_abundance_left = sprintf("%.0f",$stats_filter_a_read_total_abundance_left/$stats_reads_total*100);
my $pct_stats_filter_a1_read_total_abundance_left = sprintf("%.0f",$stats_filter_a1_read_total_abundance_left/$stats_reads_total*100);
my $pct_stats_filter_a2_read_total_abundance_left = sprintf("%.0f",$stats_filter_a2_read_total_abundance_left/$stats_reads_total*100);

my $pct_stats_a1_reads_total_potential = sprintf("%.0f",$stats_a1_reads_total_potential/$stats_reads_total*100);
my $pct_stats_a1_reads_total_barcode = sprintf("%.0f",$stats_a1_reads_total_barcode/$stats_reads_total*100);
my $stats_a1_reads_unique = keys %h_a1_index;
my $stats_a1_reads_unique_filter_removed = $stats_a1_reads_unique - $stats_a1_reads_unique_filter;
my $pct_stats_a1_reads_unique_filter = sprintf("%.0f",$stats_a1_reads_unique_filter/$stats_a1_reads_unique*100);
my $pct_stats_a1_reads_total_filter = sprintf("%.0f",$stats_a1_reads_total_filter/$stats_reads_total*100);
my $stats_a1_reads_unique_in_a1_index_removed = $stats_a1_reads_unique_filter - $stats_a1_reads_unique_in_a1_index;
my $pct_stats_a1_reads_unique_in_a1_index = sprintf("%.0f",$stats_a1_reads_unique_in_a1_index/$stats_a1_reads_unique_filter*100);
my $pct_stats_a1_reads_total_in_a1_index = sprintf("%.0f",$stats_a1_reads_total_in_a1_index/$stats_reads_total*100);
my $stats_a1_reads_keep_unique = keys %h_a1_keep;
my $pct_stats_a1_reads_keep_unique = sprintf("%.0f",$stats_a1_reads_keep_unique/$stats_a1_reads_unique_filter*100);
my $pct_stats_a1_reads_keep_total = sprintf("%.0f",$stats_a1_reads_keep_total/$stats_reads_total*100);

my $pct_stats_a2_reads_total_potential = sprintf("%.0f",$stats_a2_reads_total_potential/$stats_reads_total*100);
my $pct_stats_a2_reads_total_barcode = sprintf("%.0f",$stats_a2_reads_total_barcode/$stats_reads_total*100);
my $stats_a2_reads_unique = keys %h_a2_index;
my $stats_a2_reads_unique_filter_removed = $stats_a2_reads_unique - $stats_a2_reads_unique_filter;
my $pct_stats_a2_reads_unique_filter = sprintf("%.0f",$stats_a2_reads_unique_filter/$stats_a2_reads_unique_filter*100);
my $pct_stats_a2_reads_total_filter = sprintf("%.0f",$stats_a2_reads_total_filter/$stats_reads_total*100);
my $stats_a2_reads_unique_in_a2_index_removed = $stats_a2_reads_unique_filter - $stats_a2_reads_unique_in_a2_index;
my $pct_stats_a2_reads_unique_in_a2_index = sprintf("%.0f",$stats_a2_reads_unique_in_a2_index/$stats_a2_reads_unique_filter*100);
my $pct_stats_a2_reads_total_in_a2_index = sprintf("%.0f",$stats_a2_reads_total_in_a2_index/$stats_reads_total*100);
my $stats_a2_reads_keep_unique = keys %h_a2_keep;
my $pct_stats_a2_reads_keep_unique = sprintf("%.0f",$stats_a2_reads_keep_unique/$stats_a2_reads_unique_filter*100);
my $pct_stats_a2_reads_keep_total = sprintf("%.0f",$stats_a2_reads_keep_total/$stats_reads_total*100);

print STATS "--- Global statistics ---:\n";
print STATS "\t\t$stats_reads_total\t100%\tTotal Reads\n";

print STATS "\n--- Anchor statistics ---:\n";
print STATS "Unique\t\tTotal\t\n";
print STATS "\t\t$stats_anchor_total_potential\t$pct_stats_anchor_total_potential%\tAnchor reads\n";
print STATS "$stats_anchor_unique\t\t$stats_anchor_total_barcode\t$pct_stats_anchor_total_barcode%\tAnchor reads w. barcode\n";
print STATS "$stats_filter_index_unique_abundance_left\t100%\t$stats_filter_index_total_abundance_left\t$pct_stats_filter_index_total_abundance_left%\tAfter anchor abundance filtering ($stats_filter_index_unique_abundance index seen less than $afq times)\n";
print STATS "$stats_filter_index_unique_occurance_left\t$pct_stats_filter_index_unique_occurance_left%\t$stats_filter_index_total_occurance_left\t$pct_stats_filter_index_total_occurance_left%\tAfter anchor occurance filtering ($stats_filter_index_unique_occurance index not first sub-index observation)\n";
print STATS "$stats_filter_a_read_unique_abundance_left\t$pct_stats_filter_a_read_unique_abundance_left%\t$stats_filter_a_read_total_abundance_left\t$pct_stats_filter_a_read_total_abundance_left%\tAfter A1 + A2 read abundance filtering ($stats_filter_a_read_unique_abundance index where both reads were seen less than $sfq times)\n";
print STATS "$stats_filter_a1_read_unique_abundance_left\t$pct_stats_filter_a1_read_unique_abundance_left%\t$stats_filter_a1_read_total_abundance_left\t$pct_stats_filter_a1_read_total_abundance_left%\tAfter A1 read abundance filtering ($stats_filter_a1_read_unique_abundance index where A1 reads were seen less than $sfq times)\n";
print STATS "$stats_filter_a2_read_unique_abundance_left\t$pct_stats_filter_a2_read_unique_abundance_left%\t$stats_filter_a2_read_total_abundance_left\t$pct_stats_filter_a2_read_total_abundance_left%\tAfter A2 read abundance filtering ($stats_filter_a2_read_unique_abundance index where A2 reads were seen less than $sfq times)\n";

print STATS "\n--- A1 read statistics ---:\n";
print STATS "Unique\t\tTotal\t\n";
print STATS "\t\t$stats_a1_reads_total_potential\t$pct_stats_a1_reads_total_potential%\tA1 reads\n";
print STATS "$stats_a1_reads_unique\t\t$stats_a1_reads_total_barcode\t$pct_stats_a1_reads_total_barcode%\tA1 reads w. barcode\n";
print STATS "$stats_a1_reads_unique_filter\t100%\t$stats_a1_reads_total_filter\t$pct_stats_a1_reads_total_filter%\tAfter abundance filtering ($stats_a1_reads_unique_filter_removed A1 index were seen less than $sfq times)\n";
print STATS "$stats_a1_reads_unique_in_a1_index\t$pct_stats_a1_reads_unique_in_a1_index%\t$stats_a1_reads_total_in_a1_index\t$pct_stats_a1_reads_total_in_a1_index%\tAssociated with an anchor index ($stats_a1_reads_unique_in_a1_index_removed was not seen in any index)\n";
print STATS "$stats_a1_reads_keep_unique\t$pct_stats_a1_reads_keep_unique%\t$stats_a1_reads_keep_total\t$pct_stats_a1_reads_keep_total%\tIn passed anchors\n";

print STATS "\n--- A2 read statistics ---:\n";
print STATS "Unique\t\tTotal\t\n";
print STATS "\t\t$stats_a2_reads_total_potential\t$pct_stats_a2_reads_total_potential%\tA2 reads\n";
print STATS "$stats_a2_reads_unique\t\t$stats_a2_reads_total_barcode\t$pct_stats_a2_reads_total_barcode%\tA2 reads w. barcode\n";
print STATS "$stats_a2_reads_unique_filter\t100%\t$stats_a2_reads_total_filter\t$pct_stats_a2_reads_total_filter%\tAfter abundance filtering ($stats_a2_reads_unique_filter_removed A2 index were seen less than $sfq times)\n";
print STATS "$stats_a2_reads_unique_in_a2_index\t$pct_stats_a2_reads_unique_in_a2_index%\t$stats_a2_reads_total_in_a2_index\t$pct_stats_a2_reads_total_in_a2_index%\tAssociated with an anchor index ($stats_a2_reads_unique_in_a2_index_removed was not seen in any index)\n";
print STATS "$stats_a2_reads_keep_unique\t$pct_stats_a2_reads_keep_unique%\t$stats_a2_reads_keep_total\t$pct_stats_a2_reads_keep_total%\tIn passed anchors\n";

close STATS;

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$nice_timestamp = sprintf ( "%02d:%02d:%02d", $hour,$min,$sec);
print "$nice_timestamp Done, enjoy!\n\n";                          # To get nice time-stamps

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "R1|f:s", "R2|r:s", "stats|t:s", "afq|a:s", "sfq|s:s", "fastq|q:+", "barcodes|b:s", "maxseq|m:s", "printmap|p:+", "printraw|x:+", "onlystats|j:+");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );
    
	#if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    #if(!exists $options{'infile'} ) { print "**ERROR: $0 : \n"; exec("pod2usage $0"); }

    return \%options;
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

__DATA__

=head1 NAME

    F16S.extract_v5.0.pl

=head1 COPYRIGHT

   copyright (C) 2016 Mads Albertsen

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION



=head1 SYNOPSIS

fSSU.extract_v5.0.pl  -f -r [-h -a -s -t -q -b -m -p -x]

 [-help -h]           Displays this basic usage information.
 [-R1 -f]             Forward read in fastq format (default: R1.fq).
 [-R2 -r]             Reverse read in fastq format (default: R2.fq).
 [-barcodes -b]       File specifying barcodes to be used (default: 1 standard barcode).
 [-afq -a]            Minimum anchor frequency cutoff (default: 5).
 [-sfq -s]            Minimum sequence freaquency per sub-Anchor (default: 25).
 [-stats -t]          Stats file. 
 [-fastq -q]          Output fastq files instead of fasta (flag).
 [-maxseq -m]         Maximum number of reads to output (default: 1000).
 [-printmap -p]       Print all good anchor reads to a single file (flag).
 [-printraw -x]       Print all sub-anchor reads to a single file (flag).
 [-onlystats -j]      Only print stats (flag).

=cut
