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
    
    my $hmmer_tgz = $hmmer_dir_name;
    $hmmer_tgz =~ s/\/$/\.tgz/;

    if (!(-e $hmmer_tgz)) {
	die "\n  ERROR: Failed to locate trans-hmmer tarball '$hmmer_tgz'\n\n";
    }
    
    if (system("tar -xzf $hmmer_tgz")) {
	die "\n  ERROR: Failed to unpack trans-hmmer tarball '$hmmer_tgz'\n\n";
    }
    
    chdir($hmmer_dir_name);
    system("autoconf");
    system("./configure");
    system("make");
    
    chdir('..');
    
}

my $hmmer_src_name = $hmmer_dir_name.'src/';
system("cp $hmmsearcht_template $hmmer_src_name");
chdir($hmmer_src_name);
system("make -B");

1;
