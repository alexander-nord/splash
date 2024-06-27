#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV != 1) { die "\n  USAGE:  ./Get-Time-and-Mem-Data.pl [MP-vs-Splash-Output]\n\n"; }


my $in_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate directory '$in_dir_name'\n\n"
	if (!(-d $in_dir_name));
$in_dir_name = $in_dir_name.'/' if ($in_dir_name !~ /\/$/);


opendir(my $InDir, $in_dir_name);
my @ComparisonFiles;
while (my $file_name = readdir($InDir))
{
	if ($file_name =~ /-comparison\.csv$/)
	{
		push(@ComparisonFiles,$file_name);
	}
}
closedir($InDir);


die "\n  ERROR:  No comparison files located\n\n"
	if (scalar(@ComparisonFiles) == 0);


foreach my $comparison_file_name (sort @ComparisonFiles)
{

	$comparison_file_name =~ /^(\S+)-comparison\.csv$/;
	my $species = $1;

	$comparison_file_name = $in_dir_name.$comparison_file_name;

	open(my $CompFile,'<',$comparison_file_name)
		|| die "\n  ERROR:  Failed to open file '$comparison_file_name'\n\n";
	
	<$CompFile>; # Eat the header
	
	my $total_spl_time = 0;
	my $total_mp_time  = 0;

	my @SplTimes;
	my @SplMem;
	my @MPTimes;
	my @MPMem;
	
	while (my $line = <$CompFile>)
	{
		$line =~ s/\n|\r//g;
		next if (!$line);

		my @Data = split(/\,/,$line);

		my $spl_time = $Data[ 2];
		my $spl_mem  = $Data[ 3];
		my  $mp_time = $Data[10];
		my  $mp_mem  = $Data[11];

		if ($spl_time ne '-')
		{
			$total_spl_time += $spl_time;
			push(@SplTimes,$spl_time);
			push(@SplMem  ,$spl_mem );
		}

		if ($mp_time ne '-')
		{
			$total_mp_time += $mp_time;
			push(@MPTimes,$mp_time);
			push(@MPMem  ,$mp_mem );
		}

	}

	close($CompFile);


	@SplTimes = sort {$a <=> $b} @SplTimes;
	@SplMem   = sort {$a <=> $b} @SplMem;
	@MPTimes  = sort {$a <=> $b} @MPTimes;
	@MPMem    = sort {$a <=> $b} @MPMem;


	my $spl_min_time = $SplTimes[0];
	my $spl_min_mem  = $SplMem[0];
	my  $mp_min_time = $MPTimes[0];
	my  $mp_min_mem  = $MPMem[0];


	my $spl_25p_time = $SplTimes[int(scalar(@SplTimes)/4)];
	my $spl_25p_mem  = $SplMem[int(scalar(@SplTimes)/4)];
	my  $mp_25p_time = $MPTimes[int(scalar(@MPTimes)/4)];
	my  $mp_25p_mem  = $MPMem[int(scalar(@MPTimes)/4)];


	my $spl_med_time = $SplTimes[int(scalar(@SplTimes)/2)];
	my $spl_med_mem  = $SplMem[int(scalar(@SplTimes)/2)];
	my  $mp_med_time = $MPTimes[int(scalar(@MPTimes)/2)];
	my  $mp_med_mem  = $MPMem[int(scalar(@MPTimes)/2)];


	my $spl_75p_time = $SplTimes[int(3 * scalar(@SplTimes)/4)];
	my $spl_75p_mem  = $SplMem[int(3 * scalar(@SplTimes)/4)];
	my  $mp_75p_time = $MPTimes[int(3 * scalar(@MPTimes)/4)];
	my  $mp_75p_mem  = $MPMem[int(3 * scalar(@MPTimes)/4)];


	my $spl_max_time = $SplTimes[scalar(@SplTimes)-1];
	my $spl_max_mem  = $SplMem[scalar(@SplTimes)-1];
	my  $mp_max_time = $MPTimes[scalar(@MPTimes)-1];
	my  $mp_max_mem  = $MPMem[scalar(@MPTimes)-1];


	print "./Box-and-Whisker.py $species splash   time $spl_min_time $spl_25p_time $spl_med_time $spl_75p_time $spl_max_time\n";
	print "./Box-and-Whisker.py $species splash   mem  $spl_min_mem  $spl_25p_mem  $spl_med_mem  $spl_75p_mem  $spl_max_mem\n";
	print "./Box-and-Whisker.py $species miniprot time  $mp_min_time  $mp_25p_time  $mp_med_time  $mp_75p_time  $mp_max_time\n";
	print "./Box-and-Whisker.py $species miniprot mem   $mp_min_mem   $mp_25p_mem   $mp_med_mem   $mp_75p_mem   $mp_max_mem\n";


	my $total_spl_mins = int(10.0 * $total_spl_time / 60.0) / 10.0;
	my $total_mp_mins = int(10.0 * $total_mp_time / 60.0) / 10.0;


	print "\n";
	print "$species\n";
	print "  splashBATH :  Total Time  : $total_spl_time seconds\n";
	print "                            = $total_spl_mins minutes\n";
	print "                Median Time : $spl_med_time seconds\n";
	print "                Median Mem  : $spl_med_mem MB\n";
	print "  miniprot   :  Total Time  : $total_mp_time seconds\n";
	print "                            = $total_mp_mins minutes\n";
	print "                Median Time : $mp_med_time seconds\n";
	print "                Median Mem  : $mp_med_mem MB\n";
	print "\n";

}


1;


