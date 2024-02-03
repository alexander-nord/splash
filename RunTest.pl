#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV != 1) {
	die "\n  ./RunTest.pl [gene_dir]\n\n";
}


if (!(-d 'hmmer')) { 
	system("./Build.pl"); 
}


my $gene_dir_name = $ARGV[0];

$gene_dir_name =~ s/\/$//;
$gene_dir_name = $gene_dir_name.'/';

if (!(-d $gene_dir_name)) {
	die "\n  ERROR:  Failed to locate gene directory '$gene_dir_name'\n\n";
}


$gene_dir_name =~ /([^\/]+)\/$/;
my $gene = $1;


my $hmm_file_name = $gene_dir_name.$gene.'.hmm';
if (!(-e $hmm_file_name)) {

	my $hmmbuild = './hmmer/src/hmmbuild';
	if (!(-e $hmmbuild)) {
		die "\n  ERROR:  Failed to locate hmmbuild executable\n\n";
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


my $hmmsearch = './hmmer/src/hmmsearch';
if (!(-e $hmmsearch)) {
	die "\n  ERROR:  Failed to locate hmmsearch executable\n\n";
}


my $hmmsearch_cmd = "$hmmsearch $hmm_file_name $dna_file_name";
if (system($hmmsearch_cmd)) {
	die "\n  ERROR:  hmmsearch command '$hmmsearch_cmd' failed\n\n";
}



1;

