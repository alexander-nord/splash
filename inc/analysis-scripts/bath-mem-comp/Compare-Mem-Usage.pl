#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub GetMaxRSS;



if (@ARGV != 2) { die "\n  USAGE:  ./Compare-Mem-Usage.pl [splashBATH-Output] [BATH-Output]\n\n"; }



my $splash_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate Splash results directory '$splash_dir_name'\n\n"
	if (!(-d $splash_dir_name));
$splash_dir_name = $splash_dir_name.'/' if ($splash_dir_name !~ /\/$/);

my $bath_dir_name = $ARGV[1];
die "\n  ERROR:  Failed to locate BATH results directory '$bath_dir_name'\n\n"
	if (!(-d $bath_dir_name));
$bath_dir_name = $bath_dir_name.'/' if ($bath_dir_name !~ /\/$/);



my $splash_rbg_name = $splash_dir_name.'Results-by-Gene/';
die "\n  ERROR:  Failed to locate Splash Results-by-Gene directory '$splash_rbg_name'\n\n" 
	if (!(-d $splash_rbg_name));

my $bath_rbg_name = $bath_dir_name.'Results-by-Gene/';
die "\n  ERROR:  Failed to locate BATH Results-by-Gene directory '$bath_rbg_name'\n\n" 
	if (!(-d $bath_rbg_name));


opendir(my $SplashRBG,$splash_rbg_name)
	|| die "\n  ERROR:  Failed to open Splash Results-by-Gene directory '$splash_rbg_name'\n\n";

my @FamList;
while (my $family = readdir($SplashRBG)) 
{
	next if ($family =~ /^\./);
	$family =~ s/\/$//;

	next if (!(-d $splash_rbg_name.$family));
	next if (!(-d   $bath_rbg_name.$family));

	push(@FamList,$family);
}

closedir($SplashRBG);


die "\n  ERROR:  No overlapping families found between '$bath_rbg_name' and '$splash_rbg_name'\n\n"
	if (scalar(@FamList) == 0);


@FamList = sort @FamList;

foreach my $species (('human','chicken','stickleback','fruit_fly')) 
{

	my $out_file_name = $species.'-mem-usage.csv';
	open(my $OutFile,'>',$out_file_name)
		|| die "\n  ERROR:  Failed to open output file '$out_file_name'\n\n";

	print $OutFile "Family,BATH(MB),splash(MB),diff\n";


	my @SpeciesMemDiffs;

	foreach my $family (@FamList)
	{
		my $splash_err_name = $splash_rbg_name.$family.'/'.$family.'.'.$species.'.err';
		my $splash_max_rss  = GetMaxRSS($splash_err_name);

		my $bath_err_name = $bath_rbg_name.$family.'/'.$family.'.'.$species.'.err';
		my $bath_max_rss  = GetMaxRSS($bath_err_name);

		my $mb_diff = int(10.0 * ($splash_max_rss - $bath_max_rss))/10.0;
		print $OutFile "$family,$bath_max_rss,$splash_max_rss,$mb_diff\n";

		push(@SpeciesMemDiffs,$mb_diff);
	
	}

	close($OutFile);


	@SpeciesMemDiffs = sort {$a <=> $b} @SpeciesMemDiffs;

	my $min = $SpeciesMemDiffs[0];
	my $max = $SpeciesMemDiffs[scalar(@SpeciesMemDiffs)-1];
	my $np  = $SpeciesMemDiffs[int(90 * scalar(@SpeciesMemDiffs) / 100)];
	my $nnp = $SpeciesMemDiffs[int(99 * scalar(@SpeciesMemDiffs) / 100)];
	my $med = $SpeciesMemDiffs[int(scalar(@SpeciesMemDiffs)/2)];

	$min = ' '.$min while (length("$min") < 7);
	$med = ' '.$med while (length("$med") < 7);
	$np  = ' '.$np  while (length("$np" ) < 7);
	$nnp = ' '.$nnp while (length("$nnp") < 7);
	$max = ' '.$max while (length("$max") < 7);

	print "\n";
	print "$species\n";
	print "  Min RSS Diff    : $min MB\n";
	print "  Med RSS Diff    : $med MB\n";
	print "  90th Percentile : $np MB\n";
	print "  99th Percentile : $nnp MB\n";
	print "  Max RSS Diff    : $max MB\n";
	print "\n";

}





1;




#########################################################################################
#
#  Subroutine: GetMaxRSS
#
sub GetMaxRSS
{
	my $err_file_name = shift;

	open(my $ErrFile,'<',$err_file_name)
		|| die "\n  ERROR:  Failed to open file '$err_file_name'\n\n";

	my $max_rss_mb = '-';
	while (my $line = <$ErrFile>)
	{
		if ($line =~ /Maximum resident set size \(kbytes\): (\d+)/)
		{
			$max_rss_mb = int(10.0 * $1 / 1024) / 10.0;
			last;
		}
	}

	close($ErrFile);

	return $max_rss_mb;

}

