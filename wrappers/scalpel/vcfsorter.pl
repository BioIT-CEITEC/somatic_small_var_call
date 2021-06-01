#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

####################################################################################
# vcfsorter.pl
#
# Copyright (C) 2011 German Gaston Leparc
#
# Modified from vcfsorter to handle some goofy things ~dfj3
#
# sorts VCF by reference genome
#
# Usage:    vcfsorter.pl -d genome.dict -v input.vcf > newfile.vcf
#
# Options: -d [FILE] genome.dict file                                     [REQUIRED]
#          -v [FILE] Input variant file in vcf format                     [REQUIRED]
#          -o [FILE] Output sorted vcf file                               [OPTIONAL]
#          -l Lenient: only warn if duplicates found (Default is to fail) [OPTIONAL]
#          -r Remove duplicate entries and only report 1st (force -l)     [OPTIONAL]
#
####################################################################################

my %opts;
getopts('d:v:o:lr', \%opts);

sub usage{
    die(qq/
Usage:    vcfsorter.pl -d genome.dict -v input.vcf > newfile.vcf
Options: -d [FILE] genome.dict file                                     [REQUIRED]
         -v [FILE] Input variant file in vcf format                     [REQUIRED]
         -o [FILE] Output sorted vcf file                               [OPTIONAL]
         -l Lenient: only warn if duplicates found (Default is to fail) [OPTIONAL]
         -r Remove duplicate entries and only report 1st (force -l)     [OPTIONAL]
\n/);
}

my ($dict_file, $vcf_file, $outfile);
my $stringency = 'STRICT';

if (exists $opts{d} && exists $opts{v}){
    $dict_file = $opts{d};
    $vcf_file = $opts{v};
}
else{
    print STDERR "\nSubmit required arguments!\n";
    usage();
}

if (exists $opts{o}){
    $outfile = $opts{o};
}

if (exists $opts{l} || exists $opts{r}){
    $stringency = 'LENIENT';
    if (exists $opts{r}){
	$stringency = 'RMDUP';
    }
}

#---------------------------------------- LOAD IN FASTA DICT INTO MEMORY
open(DICT,$dict_file) or die "Can't open $dict_file!\n";
my @contig_order;
my $c=0;
while(<DICT>){
    if($_=~ /\@SQ/){
        my ($contig) = $_ =~ /SN:(\S+)/;
        $contig_order[$c]=$contig;
        ++$c;
        #print $contig,"\n";
    }
}
close(DICT);

#---------------------------------------- PARSE VCF FILE & OUTPUT SORTED VCF

open(VCF,$vcf_file) or die "Can't open $vcf_file!\n";

my %vcf_hash;
my $header = '';

while(<VCF>){
    if($_=~/^\#/){
        $header .= $_; # store header and comment fields
        next;
    }

    chomp($_);

    my @data = split(/\t/,$_);
    my $contig = $data[0];
    my $start = $data[1];
    my $variant = $data[3]."to".$data[4];
    my $line = $_;

    #print $contig,":",$start,"\n";
    if(exists $vcf_hash{$contig}{$start}{$variant}){
        if($stringency eq 'STRICT'){
            print STDERR "ERROR: Duplicate entry found!\n",$line,"\n",$vcf_hash{$contig}{$start}{$variant},"\n";
            die;
        }
        elsif($stringency eq 'LENIENT'){
            print STDERR "WARNING: Duplicate entry found! Records will be reported in the order they were in the file!\n";
            $vcf_hash{$contig}{$start}{$variant} .= "\n" . $line;
        }
        elsif($stringency eq 'RMDUP'){
            print STDERR "WARNING: Duplicate entry found! Only the first entry will be reported!\n";
        }
        else{
            print STDERR "ERROR: Sorry, malformed stringency\n";
            die;
        }
    }
    else{
        $vcf_hash{$contig}{$start}{$variant}=$line;
    }
}
close(VCF);

#------------------ print out the VCF in the order of the reference genome

my $out;
if (defined $outfile){
    open($out, ">$outfile");
}

#print standard VCF header
if (defined $outfile){
    print $out $header;
}
else{
    print $header;
}

foreach my $contig (@contig_order){ # sort by contig order
    #print $contig,"\n";
    foreach my $start (sort {$a <=> $b} keys %{$vcf_hash{$contig}}){ # sort numerically by coordinates
        #print $start,"\n";
        foreach my $variant (keys %{$vcf_hash{$contig}{$start}}){ # if overlapping mutation, print each variant
            if (defined $outfile){
                print $out $vcf_hash{$contig}{$start}{$variant},"\n";
            }
            else{
                print $vcf_hash{$contig}{$start}{$variant},"\n";
            }
        }
    }
}

