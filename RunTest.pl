#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV != 1) {
    die "\n  ./RunTest.pl [gene_dir]\n\n";
}


my $hmmer_dir_name = './trans-hmmer/';

my $hmmbuild   = $hmmer_dir_name.'src/hmmbuild';
my $hmmsearcht = $hmmer_dir_name.'src/hmmsearcht';


if (!(-d $hmmer_dir_name)) { 
    system("./Build.pl");
    if (!(-e $hmmbuild)) {
	die "\n  ERROR:  hmmer build appears to have failed\n\n";
    }
}


my $gene_dir_name = $ARGV[0];

$gene_dir_name =~ s/\/$//;
$gene_dir_name = $gene_dir_name.'/';

if (!(-d $gene_dir_name)) {

    # In case only the gene name was provided, we'll
    # see if there's a human entry for this gene.
    my $test_gene_dir_name = 'inputs/homo_sapiens/'.$gene_dir_name;

    if (!(-d $test_gene_dir_name)) {
	die "\n  ERROR:  Failed to locate gene directory '$gene_dir_name'\n\n";
    }

    $gene_dir_name = $test_gene_dir_name;
    
}


$gene_dir_name =~ /([^\/]+)\/$/;
my $gene = $1;


my $hmm_file_name = $gene_dir_name.$gene.'.hmm';
if (!(-e $hmm_file_name)) {
    
    if (!(-e $hmmbuild)) {
	die "\n  ERROR:  Failed to locate hmmbuild executable ($hmmbuild)\n\n";
    }
    
    my $msa_file_name = $gene_dir_name.'isos.afa';
    
    if (!(-e $msa_file_name)) {
	die "\n  ERROR:  Failed to locate alignment file '$msa_file_name'\n\n";
    }
    
    my $build_cmd = "$hmmbuild $hmm_file_name $msa_file_name";
    if (system("$build_cmd")) {
	die "\n  ERROR:  hmmbuild command '$build_cmd' failed\n\n";
    }
    
}

if (!(-e $hmm_file_name)) {
    die "\n  ERROR:  Failed to locate HMM file '$hmm_file_name'\n\n";
}



my $dna_file_name = $gene_dir_name.'dna.fa';
    
if (!(-e $dna_file_name)) {
    die "\n  ERROR:  Failed to locate DNA file '$dna_file_name'\n\n";
}



if (!(-e $hmmsearcht)) {
    die "\n  ERROR:  Failed to locate hmmsearcht executable ($hmmsearcht)\n\n";
}



my $hmmsearcht_cmd = "$hmmsearcht $hmm_file_name $dna_file_name";
if (system($hmmsearcht_cmd)) {
    die "\n  ERROR:  hmmsearcht command '$hmmsearcht_cmd' failed\n\n";
}



1;

