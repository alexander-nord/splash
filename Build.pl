#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV == 1 && lc($ARGV[0]) =~ /^\-?\-?miniprot/)
{
    system("git clone https://github.com/lh3/miniprot");
    chdir("miniprot");
    system("make");
    exit(0);
}



my $splash_bathsearch = 'bathsearch.c';
if (!(-e $splash_bathsearch)) {
    die "\n  ERROR: I can't find bathsearch.c!\n\n";
}


my $splash_pipeline = 'inc/p7_pipeline.c';
if (!(-e $splash_pipeline)) {
    die "\n  ERROR: There should be a 'p7_pipeline.c' file in 'inc/'? (mainly force-disables frameshift-aware pipeline for splash)\n\n";
}


my $splash_hmmer_h = 'inc/hmmer.h';
if (!(-e $splash_hmmer_h)) {
    die "\n  ERROR: There should be a 'hmmer.h' file in 'inc/'\n\n";
}


my $bath_dir_name = 'BATH/';


if (!(-d $bath_dir_name)) {

    my $bath_git_url = 'https://github.com/traviswheelerlab/BATH';
    if (system("git clone $bath_git_url")) {
	die "\n  ERROR: Failed to clone into $bath_git_url\n\n";
    }
    
    chdir($bath_dir_name);

    system("git clone https://github.com/TravisWheelerLab/easel");
    chdir("easel/");
    system("git checkout BATH");
    chdir("../");

    system("autoconf");
    system("./configure");

    system("cp ../$splash_bathsearch src/");
    system("cp ../$splash_pipeline   src/");
    system("cp ../$splash_hmmer_h    src/");
    system("make");
    
} else {

    my $bath_src_name = $bath_dir_name.'src/';
    system("cp $splash_bathsearch $bath_src_name");
    system("cp $splash_pipeline   $bath_src_name");
    system("cp $splash_hmmer_h    $bath_src_name");
    chdir($bath_src_name);
    system("make -B");

}

1;
