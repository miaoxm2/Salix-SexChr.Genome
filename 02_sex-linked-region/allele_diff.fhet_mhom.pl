#!/usr/bin/perl

#TODO: Find the alleles that heter in females and homo males. and the alleles of male belong to two alles of females
#USAGE:
#perl $script [input_file] [num_females] [num_males] > [output]
#INPUT file should be females first and male first. after vcftools filter, gatk generate gt 
#previous step: gatk VariantsToTable -F CHROM -F POS -F TYPE -GF GT -GF AD -V $filter2.bi.recode.vcf -O $filter2.bi.gt.tbl


#DATE: 2024-06-17
#UPDATED: 2024-11-27
#Xiaomeng

use strict;
use warnings;

my $num=10;

if (@ARGV < 3) {
    die "Usage: $0 input_file num_females num_males\n";
}


my $file = $ARGV[0]; #GT file to process
# Open the file for reading
open F1, '<', $file or die "Could not open file '$file': $!\n";

my $fnum = $ARGV[1]; #the number of female sampels, first!
my $mnum = $ARGV[2]; #the number of male samples
my $pos=3; #starting column of the first sample's GT. Usually the forth.


# Process each line of the GT file
while (my $line = <F1>) {
    chomp $line;
    # Skip header lines
    next if $line =~ /^CHROM/;

    # Split the line into fields
    my @fields = split /\t/, $line;

    my $is_unique = 0;
    my $female=0;
    my %allele;
    for (my $i = $pos; $i < $pos+$fnum*2; $i=$i+2) {

        my ($gt1,$gt2) = (split/[\/|]/, $fields[$i])[0,1];
        if ($gt1 ne $gt2) {
            $female++;
            $allele{$gt1}++;
            $allele{$gt2}++;

        }
    }
    if ($female==$fnum) {
    for (my $i = $pos+$mnum*2; $i < $pos+$mnum*4; $i=$i+2) {
        my ($gt1,$gt2) = (split/[\/|]/, $fields[$i])[0,1];
        if ($gt1 eq $gt2 && $gt1 ne '.' && $gt2 ne '.' ) {
            if (exists $allele{$gt1}) { 
                $is_unique++;
            }
        }

    }
}
    if ($is_unique==$mnum) {
        print "$line\n";
           
    }
}

# Close the VCF file
close F1;

