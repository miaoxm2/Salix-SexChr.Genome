#!/usr/bin/perl

#USAGE: Find the alleles that only exsits in females but not in males.
#INPUT file should be females first and male first. after vcftools filter, gatk generate dp 
#DATE: 2024-06-11
#Xiaomeng

use strict;
use warnings;

# Open the VCF file for reading
open  F1, "<","$ARGV[0]" or die "Cannot open file: $!";
my $pos=4;
my $num=3;
# Process each line of the VCF file
while (my $line = <F1>) {
    chomp $line;
    # Skip header lines
    next if $line =~ /^CHROM/;

    # Split the line into fields
    my @fields = split /\t/, $line;

    my $fdp=0;
    my $mdp=0;
    my %female;
    for (my $i = $pos; $i < $pos+$num*2; $i=$i+2) {
        my ($dp1,$dp2) = (split/,/, $fields[$i])[0,1];
        my $dp=$dp1+$dp2;
        if ($dp>0) {
        	$fdp++;
        }
    }
    if ($fdp==$num) {
    	for (my $i = $pos+$num*2; $i < $pos+$num*4; $i=$i+2) {
        	my ($dp1,$dp2) = (split/,/, $fields[$i])[0,1];
        	my $dp=$dp1+$dp2;
        	if ($dp==0) {
            	$mdp++;
       		}
       	}
    }
    if ($mdp==$num) {
    	print "$line\n";
    }
}

# Close the VCF file
close F1;


