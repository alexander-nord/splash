#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


my $hmmsearcht_template = 'hmmsearcht.c';
if (!(-e $hmmsearcht_template)) {
    die "\n  ERROR: I can't find hmmsearcht.c!\n\n";
}


my $hmmer_dir_name = 'trans-hmmer/';


if (!(-d $hmmer_dir_name)) {

    my $trans_hmmer_git_url = 'https://github.com/traviswheelerlab/trans-hmmer';
    if (system("git clone $trans_hmmer_git_url")) {
	die "\n  ERROR: Failed to clone into $trans_hmmer_git_url\n\n";
    }
    
    chdir($hmmer_dir_name);
    system("git submodule init");
    system("git submodule update --remote");
    system("cp ../alt-configure.ac configure.ac");
    system("autoconf");
    system("./configure");
    system("cp ../$hmmsearcht_template src/");
    system("make");
    
} else {

    my $hmmer_src_name = $hmmer_dir_name.'src/';
    system("cp $hmmsearcht_template $hmmer_src_name");
    chdir($hmmer_src_name);
    system("make -B");

}

1;
