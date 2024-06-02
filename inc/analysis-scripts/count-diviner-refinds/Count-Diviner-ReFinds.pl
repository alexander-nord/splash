#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub ForwardOverlap;
sub ReverseOverlap;



if (@ARGV != 2) { die "\n  USAGE: ./Count-Diviner-ReFinds.pl [Splash-Diviner-Output/] [diviner-based-input/]\n\n"; }



my $spl_results_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate Splash results directory '$spl_results_dir_name'\n\n"
	if (!(-d $spl_results_dir_name));

$spl_results_dir_name = $spl_results_dir_name.'/' if ($spl_results_dir_name !~ /\/$/);


my $rbg_dir_name = $spl_results_dir_name.'Results-by-Gene/';
die "\n  ERROR:  Failed to locate 'Results-by-Gene' directory '$rbg_dir_name'\n\n"
	if (!(-d $rbg_dir_name));


my $input_dir_name = $ARGV[1];
die "\n  ERROR:  Failed to locate Diviner-based input directory '$input_dir_name'\n\n"
	if (!(-d $input_dir_name));

$input_dir_name = $input_dir_name.'/' if ($input_dir_name !~ /\/$/);



opendir(my $InputDir,$input_dir_name)
	|| die "\n  ERROR:  Failed to open input directory '$input_dir_name'\n\n";

my $num_covered = 0;
my $num_missed  = 0;

while (my $gene = readdir($InputDir))
{

	next if ($gene =~ /^\./);
	$gene =~ s/\/$//;


	my $gene_targets_file_name = $input_dir_name.$gene.'/novel-exon-ranges.out';
	next if (!(-e $gene_targets_file_name));


	my $spl_hit_file_name = $rbg_dir_name.$gene.'/'.$gene.'.homo_sapiens.out';
	next if (!(-e $spl_hit_file_name));
	

	open(my $GeneTargetsFile,'<',$gene_targets_file_name)
		|| die "\n  ERROR:  Failed to open gene target list file '$gene_targets_file_name'\n\n";

	
	my $target_chr = <$GeneTargetsFile>;
	$target_chr =~ /CHR\s+:\s+(\S+)/;
	$target_chr = $1;

	my $revcomp = 0;
	if ($target_chr =~ /\[revcomp\]/)
	{
		$revcomp = 1;
		$target_chr =~ s/\[revcomp\]//;
	}


	my @DivRanges;
	my @Covered;
	while (my $line = <$GeneTargetsFile>)
	{
		if ($line =~ /RANGE : (\S+)/)
		{
			push(@DivRanges,$1);
			push(@Covered,0);
		}
	}


	close($GeneTargetsFile);



	open(my $SplHits,'<',$spl_hit_file_name)
		|| die "\n  ERROR:  Failed to open Splash hit file '$spl_hit_file_name'\n\n";

	while (my $line = <$SplHits>)
	{

		next if ($line !~ /Exon Set \d+ \((\d+)/);
		my $num_exons = $1;


		<$SplHits>; # Model Positions


		$line = <$SplHits>;
		$line =~ /Target Seq Name (\S+)/;

		my $chr = $1;
		next if ($chr ne $target_chr);


		my $line = <$SplHits>;
		$line =~ /Coords (\d+)\.\.(\d+)/;
		my $full_range_start = $1;
		my $full_range_end   = $2;

		if ($revcomp)
		{
			next if ($full_range_start < $full_range_end);
		}
		else
		{
			next if ($full_range_start > $full_range_end);
		}



		for (my $exon_id=1; $exon_id<=$num_exons; $exon_id++)
		{

			$line = <$SplHits>;
			$line =~ /\/ (\d+\.\.\d+)/;
			my $spl_exon_range = $1;

			for (my $target_id=0; $target_id<scalar(@DivRanges); $target_id++)
			{

				next if ($Covered[$target_id]);
			
				if ( $revcomp && ReverseOverlap($DivRanges[$target_id],$spl_exon_range))
				{
					$Covered[$target_id] = 1;
					last;
				}
				if (!$revcomp && ForwardOverlap($DivRanges[$target_id],$spl_exon_range))
				{
					$Covered[$target_id] = 1;
					last;
				}

			}

		}

	}

	close($SplHits);


	for (my $target_id=0; $target_id<scalar(@DivRanges); $target_id++)
	{
		if ($Covered[$target_id]) 
		{
			$num_covered++;
		}
		else
		{
			$num_missed++;
		}
	}


}

closedir($InputDir);



print "  $num_covered Covered\n";
print "  $num_missed Missed\n";



1;








sub ForwardOverlap
{
	my $range1 = shift;
	my $range2 = shift;

	$range1 =~ /^(\d+)\.\.(\d+)$/;
	my $start1 = $1;
	my $end1   = $2;

	$range2 =~ /^(\d+)\.\.(\d+)$/;
	my $start2 = $1;
	my $end2   = $2;


	return 1 if ($start1 <= $start2 && $end1 >= $start2);
	return 1 if ($start2 <= $start1 && $end2 >= $start1);
	return 1 if ($start1 >= $start2 && $end1 <= $end2);
	return 1 if ($start2 >= $start1 && $end2 <= $end1);

	return 0;
}

sub ReverseOverlap
{
	my $range1 = shift;
	my $range2 = shift;

	$range1 =~ /^(\d+)\.\.(\d+)$/;
	my $start1 = $1;
	my $end1   = $2;

	$range2 =~ /^(\d+)\.\.(\d+)$/;
	my $start2 = $1;
	my $end2   = $2;


	return 1 if ($start1 >= $start2 && $end1 <= $start2);
	return 1 if ($start2 >= $start1 && $end2 <= $start1);
	return 1 if ($start1 <= $start2 && $end1 >= $end2);
	return 1 if ($start2 <= $start1 && $end2 >= $end1);

	return 0;
}




