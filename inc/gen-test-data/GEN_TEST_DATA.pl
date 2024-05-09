#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

    

# Check that we're confident in the locations of the helper tools
# before anything else
my $Mirage_to_Splash  = 'MirageOutToSplashIn.pl';
die "\n  ERROR:  Failed to locate Mirage2 result parsing script ($Mirage_to_Splash)\n\n"
    if (!(-e $Mirage_to_Splash));
my $genome_downloader = 'DownloadHMRGenomes.sh';
die "\n  ERROR:  Failed to locate genome downloading script ($genome_downloader)\n\n"
    if (!(-e $genome_downloader));




if (@ARGV != 2) { die "\n  USAGE:  ./GEN_TEST_DATA.pl [Mirage2-Results] [out-dir-name]\n\n"; }




my $out_dir_name = $ARGV[1];
if (-e $out_dir_name) {
    die "\n  ERROR:  Directory '$out_dir_name' already exists\n\n";
} elsif (system("mkdir \"$out_dir_name\"")) {
    die "\n  ERROR:  Failed to create directory '$out_dir_name'\n\n";
}


my $mirage_results_dir_name = $ARGV[0];
$mirage_results_dir_name = $mirage_results_dir_name.'/'
    if ($mirage_results_dir_name !~ /\/$/);
if (!(-d $mirage_results_dir_name.'Species-MSAs/')) {
    die "\n  ERROR:  Directory '$mirage_results_dir_name' does not have a 'Species-MSAs' directory...?\n\n";
}


if (system("perl $Mirage_to_Splash \"$mirage_results_dir_name\"")) {
    system("rm -rf \"$out_dir_name\"");
    die "\n  ERROR:  Failed to parse Mirage2 results from directory '$mirage_results_dir_name'\n\n";
}


if (system("./$genome_downloader")) {
    system("rm -rf \"$out_dir_name\"");
    die "\n  ERROR:  Failed to download (or unpack) genomic sequence data\n\n";
}


system("mv \"protein-data\"      \"$out_dir_name\""         );
system("mv \"genomes\"           \"$out_dir_name\""         );
system("mv \"species-to-genome\" \"$out_dir_name\/genomes\"");



print "\n  Test data generated!\n  Direct splash testing towards '$out_dir_name'\n\n";


1;
