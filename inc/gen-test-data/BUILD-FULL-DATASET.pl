#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub GetScriptDir { return './' if ($0 !~ /\//); $0 =~ /^(.+\/)[^\/]+$/; return $1; }
my $SCRIPT_DIR = GetScriptDir();

    

# Check that we're confident in the locations of the helper tools
# before anything else
my $Mirage_to_Splash  = $SCRIPT_DIR.'MirageOutToSplashIn.pl';
die "\n  ERROR:  Failed to locate Mirage2 result parsing script ($Mirage_to_Splash)\n\n"
    if (!(-e $Mirage_to_Splash));

my $genome_downloader = $SCRIPT_DIR.'DownloadHMRGenomes.sh';
die "\n  ERROR:  Failed to locate genome downloading script ($genome_downloader)\n\n"
    if (!(-e $genome_downloader));

my $PANTHER_formatter = $SCRIPT_DIR.'FormatPANTHER.pl';
die "\n  ERROR:  Failed to locate PANTHER formatting script ($PANTHER_formatter)\n\n"
    if (!(-e $PANTHER_formatter));

my $PANTHER_downloader = $SCRIPT_DIR.'DownloadPANTHER.sh';
die "\n  ERROR:  Failed to locate PANTHER downloading script ($PANTHER_downloader)\n\n"
    if (!(-e $PANTHER_downloader));




if (@ARGV < 2 || @ARGV > 3) { die "\n  USAGE:  ./BUILD-FULL-DATASET.pl {OPT:--panther} [Mirage2-Results] [out-dir-name]\n\n"; }



my $panther = 0;
my $out_dir_name;
my $mirage_results_dir_name;
if (scalar(@ARGV) == 2) {
    $mirage_results_dir_name = $ARGV[0];
    $out_dir_name            = $ARGV[1];
} elsif (lc($ARGV[0]) =~ /^\-?\-?panther$/) {
    $panther = 1;
    $mirage_results_dir_name = $ARGV[1];
    $out_dir_name            = $ARGV[2];
} else {
    die "\n  ERROR:  Failed to make sense of commandline arguments...?\n\n";
}


$out_dir_name =~ s/\/$//;


if (-e $out_dir_name) {
    die "\n  ERROR:  Directory '$out_dir_name' already exists\n\n";
} elsif (system("mkdir \"$out_dir_name\"")) {
    die "\n  ERROR:  Failed to create directory '$out_dir_name'\n\n";
}


$mirage_results_dir_name = $mirage_results_dir_name.'/'
    if ($mirage_results_dir_name !~ /\/$/);
if (!(-d $mirage_results_dir_name.'Species-MSAs/')) {
    die "\n  ERROR:  Directory '$mirage_results_dir_name' does not have a 'Species-MSAs' directory...?\n\n";
}


if (system("perl $Mirage_to_Splash \"$mirage_results_dir_name\"")) {
    system("rm -rf \"$out_dir_name\"");
    die "\n  ERROR:  Failed to parse Mirage2 results from directory '$mirage_results_dir_name'\n\n";
}


if (system("$genome_downloader")) {
    system("rm -rf \"$out_dir_name\"");
    die "\n  ERROR:  Failed to download (or unpack) genomic sequence data\n\n";
}


system("mv \"mirage-based-data\" \"$out_dir_name\"");
system("mv \"genomes\"           \"$out_dir_name\"");




if ($panther) {

    die "\n  ERROR:  Failed to download (or unpack) the PANTHER18.0 HMM database\n\n"
        if (system("$PANTHER_downloader"));

    die "\n  ERROR:  Failed to locate PANTHER dataset (likely unexpected naming, looking for 'PANTHER18.0')\n\n"
        if (!(-d 'PANTHER18.0'));    

    die "\n  ERROR:  Failed to process PANTHER into easily Splash-able format\n\n"
        if (system("perl $PANTHER_formatter \"PANTHER18.0\" \"$out_dir_name/panther-data\""));

    system("rm -rf \"PANTHER18.0\"");

}




print "\n  Test data generated!\n  Direct splash testing towards '$out_dir_name'\n\n";


1;
