#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


my $hmmsearch_template = 'hmmsearch.c';
if (!(-e $hmmsearch_template)) {
	die "\n  F-ING H!: I can't find hmmsearch.c!\n\n";
}


my $hmmer_dir_name = 'hmmer/';
my $hmmer_src_name = $hmmer_dir_name.'src/';


if (!(-d $hmmer_dir_name)) {

	my $hmmer_tar_name = 'hmmer.tar';

	if (!(-e $hmmer_tar_name)) {
		die "\n  F-ING H!: I can't find HMMER tarball ($hmmer_tar_name)\n\n";
	}

	if (system("tar -xf $hmmer_tar_name")) {
		die "\n  ERROR: Failed to un-tar $hmmer_tar_name\n\n";
	}

	chdir($hmmer_dir_name);
	system("autoconf");
	system("./configure");
	system("make");

	chdir('..');

}

system("cp $hmmsearch_template $hmmer_src_name");
chdir($hmmer_src_name);
system("make -B");

1;