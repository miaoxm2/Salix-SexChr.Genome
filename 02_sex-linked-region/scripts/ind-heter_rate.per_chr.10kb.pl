#!/usr/bin/perl
#USAGE: ind-heter_rate.per_chr.10kb.pl
#count the number of heterozygous per 10kb in each individual
#INPUT file should be filtered vcf splitted into each chr. create individual line manually.
#DATE: 2024-08-26
#Xiaomeng

use List::Util qw(max);

# Define input VCF file
open  F1, "<","$ARGV[0]" or die "Cannot open file: $!";

# Initialize variables
my %heterozygosity_counts;
my $current_chr = '';
my $current_window_start = 1;
my $current_window_end = $window_size;

#create headline
@individuals=("She1-F05", "She1-F06", "She1-F08", "She1-F09", "She1-F11", "She1-F16", "She1-F22", "She1-F32", "She1-F36", "She1-F40", "She1-M05", "She1-M06", "She1-M08", "She1-M09", "She1-M11", "She1-M16", "She1-M22", "She1-M32", "She1-M36", "She1-M40");
print "chr\tcoord\t";
foreach $individual(@individuals){
    print "$individual\t";
}
print "\n";

#read vcf file
while (my $line = <F1>) {
    chomp($line);
    next if $line =~ /^#/;
    my @fields = split(/\t/, $line);
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genotypes) = @fields;
        my $kb=int($pos/10000)*10000; #create window as 10kb.
        push @kb,$kb;
        # Count heterozygosity in each individual
        for my $i (0 .. $#genotypes) {
            my ($gt_call) = split(':', $genotypes[$i]);
            if ($gt_call =~ /0\/1|1\/0/) {
                my $ind_pos = join('#',$individuals[$i],$kb);
                $heterozygosity_counts{$ind_pos}++;
            }
    }
}

#use coord as ref, print each coord and ind.
for (my $var = 0; $var < (max @kb)+1; $var=$var+10000) {
    print "$var\t";
    foreach $individual(@individuals){
        my $ind_pos = join('#',$individual,$var);
        print "$heterozygosity_counts{$ind_pos}\t";
    }
    print "\n";
}
