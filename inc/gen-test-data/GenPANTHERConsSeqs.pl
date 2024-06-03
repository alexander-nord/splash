#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


my $hmmer_dir_name = $0;
$hmmer_dir_name =~ s/[^\/]+$/hmmer\//;

if (!(-d $hmmer_dir_name))
{
	use Cwd qw();
	my $pwd = Cwd::cwd();

	system("git clone https://github.com/EddyRivasLab/hmmer \"$hmmer_dir_name\"");

	chdir($hmmer_dir_name);
	system("git clone https://github.com/EddyRivasLab/easel");
	system("autoconf");
	system("./configure");
	system("make");

	chdir($pwd);
}


my $HMMEMIT = $hmmer_dir_name.'src/hmmemit';
die "\n  ERROR:  Failed to locate hmmemit (looking for '$HMMEMIT')\n\n"
	if (!(-e $HMMEMIT));



if (@ARGV != 1) { die "\n  USAGE:  ./GenPANTHERConsSeqs.pl [PANTHER-data/]\n\n"; }



my $panther_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate PANTHER directory '$panther_dir_name'\n\n"
	if (!(-d $panther_dir_name));
$panther_dir_name = $panther_dir_name.'/' if (!(-d $panther_dir_name));


opendir(my $PantherDir,$panther_dir_name)
	|| die "\n  ERROR:  Failed to open PANTHER directory '$panther_dir_name'\n\n";

while (my $family = readdir($PantherDir))
{
	next if ($family =~ /^\./);
	$family =~ s/\/$//;

	my $fam_dir_name = $panther_dir_name.$family.'/';

	my $fam_hmm_file_name = $fam_dir_name.$family.'.hmm';
	next if (!(-e $fam_hmm_file_name));

	my $cons_seq_file_name = $fam_dir_name.$family.'.cons.fa';
	
	my $hmmemit_cmd = $HMMEMIT." -c -o \"$cons_seq_file_name\" \"$fam_hmm_file_name\"";

	print "  WARNING:  hmmemit command '$hmmemit_cmd' failed!\n"
		if (system($hmmemit_cmd));

}

closedir($PantherDir);


1;
