#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub CheckRangeOverlap;
sub ExtractSplashSpliceCoords;



if (@ARGV != 2) { die "\n  USAGE:  ./Compute-Coord-Diffs.pl [Splash-Validation-Output/] [mirage-based-input/]\n\n"; }


my $splash_validation_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate Splash output directory '$splash_validation_dir_name'\n\n"
	if (!(-d $splash_validation_dir_name));

$splash_validation_dir_name = $splash_validation_dir_name.'/' 
	if ($splash_validation_dir_name !~ /\/$/);


my $mirage_based_dir_name = $ARGV[1];
die "\n  ERROR:  Failed to locate Mirage2-based input directory '$mirage_based_dir_name'\n\n"
	if (!(-d $mirage_based_dir_name));

$mirage_based_dir_name = $mirage_based_dir_name.'/'
	if ($mirage_based_dir_name !~ /\/$/);


# Generate a CSV file with the coordinate ranges of the exons,
# as described Splash.
#
my $splash_coords_file_name = ExtractSplashSpliceCoords($splash_validation_dir_name);


# Get the number of times we see each difference in splice junction
# position (number of nucleotides), as well as the number of times
# we have Splash exons that don't overlap with Mirage2 exons.
#
my ($diff_counts_ref,$num_non_overlapping) 
	= GetOverlapData($splash_coords_file_name,$mirage_based_dir_name);

my %DiffCounts = %{$diff_counts_ref};


my $diff_counts_file_name = 'junction-nucl-diff-counts.csv';

open(my $DiffCountsFile,'>',$diff_counts_file_name)
	|| die "\n  ERROR:  Failed to open output file '$diff_counts_file_name'\n\n";


print $DiffCountsFile "Nucleotide-Difference,Exon-Count\n";

foreach my $diff_amount (sort {$a <=> $b} keys %DiffCounts)
{
	my $diff_count = $DiffCounts{$diff_amount};
	print $DiffCountsFile "$diff_amount,$diff_count\n";
}

close($DiffCountsFile);


print "\n $num_non_overlapping Splash exon";
print "s" if ($num_non_overlapping != 1);
print " had no overlap with a Mirage exon\n\n";



1;












##############################################################################
#
#  Subroutine: GetOverlapData
#
sub GetOverlapData
{

	$splash_coords_file_name = shift;
	$mirage_based_dir_name   = shift;


	opendir(my $MirageBasedDir,$mirage_based_dir_name)
		|| die "\n  ERROR:  Failed to open Mirage2-based input directory '$mirage_based_dir_name'\n\n";


	my %DiffCounts;
	my $num_non_overlapping = 0;


	while (my $gene = readdir($MirageBasedDir))
	{

		next if ($gene =~ /^\./);
		$gene =~ s/\/$//;

		my $gene_dir_name = $mirage_based_dir_name.$gene.'/';
		next if (!(-d $gene_dir_name));


		opendir(my $GeneDir,$gene_dir_name)
			|| die "\n  ERROR:  Failed to open gene directory '$gene_dir_name'\n\n";

		while (my $filename = readdir($GeneDir))
		{
			

			next if ($filename !~ /^(\S+)\.map\.out$/);
			my $seq_id = $1;


			# This is the lazy way! Yay!
			my $grep_cmd = "grep '$gene,$seq_id,' $splash_coords_file_name |";
			open(my $Grep,$grep_cmd)
				|| die "\n  ERROR:  Failed to grep! (command: $grep_cmd)\n\n";

			my @SplashCoordRanges;
			while (my $line = <$Grep>) 
			{
				if ($line =~ /,(\d+),(\d+)\s*$/)
				{
					push(@SplashCoordRanges,$1.'..'.$2);
				}
			}

			close($Grep);


			next if (scalar(@SplashCoordRanges) == 0);


			# Now pull in the actual (well, Mirage) coordinates
			$filename = $gene_dir_name.$filename;
			open(my $MapFile,'<',$filename)
				|| die "\n  ERROR:  Failed to open mapping coordinate file '$filename'\n\n";

			<$MapFile>; # Original seq name
			<$MapFile>; # Chromosome
			<$MapFile>; # Number of exons (we'll get this implicitly)

			my @MirageCoordRanges;
			while (my $line = <$MapFile>)
			{
				if ($line =~ /^(\d+\.\.\d+)$/)
				{
					push(@MirageCoordRanges,$1);
				}
			}

			close($MapFile);


			# Time to compare!
			for (my $spl_exon_id=0; $spl_exon_id<scalar(@SplashCoordRanges); $spl_exon_id++)
			{

				my $splash_coord_range = $SplashCoordRanges[$spl_exon_id];
				my $overlap_found = 0;
				my $start_diff;
				my $end_diff;


				foreach my $mirage_coord_range (@MirageCoordRanges)
				{
					($overlap_found,$start_diff,$end_diff) 
						= CheckRangeOverlap($splash_coord_range,$mirage_coord_range);
					last if ($overlap_found == 1);
				}


				if ($overlap_found)
				{

					# We don't include the start of the first exon
					if ($spl_exon_id > 0)
					{
						if ($DiffCounts{$start_diff}) { $DiffCounts{$start_diff}++; }
						else                          { $DiffCounts{$start_diff}=1; }
					}

					# We don't include the end of the final exon
					if ($spl_exon_id < scalar(@SplashCoordRanges)-1)
					{
						if ($DiffCounts{$end_diff}) { $DiffCounts{$end_diff}++; }
						else                        { $DiffCounts{$end_diff}=1; }
					}

				}
				else
				{
					$num_non_overlapping++;
				}

			}


		}

		closedir($GeneDir);

	}

	closedir($MirageBasedDir);


	return (\%DiffCounts,$num_non_overlapping);

}









##############################################################################
#
#  Subroutine: CheckRangeOverlap
#
sub CheckRangeOverlap
{
	
	my $range1 = shift;
	my $range2 = shift;


	$range1 =~ /^(\d+)\.\.(\d+)$/;
	my $start1 = $1;
	my $end1   = $2;

	$range2 =~ /^(\d+)\.\.(\d+)$/;
	my $start2 = $1;
	my $end2   = $2;


	# pre-compute the differences
	my $diff1 = abs($start1 - $start2);
	my $diff2 = abs($end1   - $end2  );


	# Just to keep the logic human-readable, let's
	# split things up based on strand direction.
	if ($start1 < $end1)
	{
		if ($start1 <= $start2 && $end1 > $start2) 
		{
			return (1,$diff1,$diff2);
		}
		elsif ($start1 < $end2 && $end1 >= $end2)
		{
			return (1,$diff1,$diff2);
		}
	}
	else
	{
		if ($start1 >= $start2 && $end1 < $start2) 
		{
			return (1,$diff1,$diff2);
		}
		elsif ($start1 > $end2 && $end1 <= $end2)
		{
			return (1,$diff1,$diff2);
		}
	}

	return (0,0,0);

}







##############################################################################
#
#  Subroutine: ExtractSplashSpliceCoords
#
sub ExtractSplashSpliceCoords
{

	my $splash_validation_dir_name = shift;

	my $rbg_dir_name = $splash_validation_dir_name.'Results-by-Gene/';
	opendir(my $RBGDir,$rbg_dir_name)
		|| die "\n  ERROR:  Failed to open Results-by-Gene directory '$rbg_dir_name'\n\n";

	# First pass is so we can sort...
	my %GenesToDirNames;
	while (my $gene = readdir($RBGDir)) 
	{

		next if ($gene =~ /^\./);
		$gene =~ s/\/$//;

		my $gene_dir_name = $rbg_dir_name.$gene.'/';

		next if (!(-e $gene_dir_name.'summary.out'));

		$GenesToDirNames{$gene} = $gene_dir_name;

	}
	closedir($RBGDir);


	my $splash_coords_file_name = 'splash-splice-coords.csv';
	open(my $OutFile,'>',$splash_coords_file_name)
		|| die "\n  ERROR:  Failed to open Splash coordinates output file '$splash_coords_file_name'\n\n";


	# Header
	print $OutFile "Gene,Input-ID,Exon-Identifiers,Chromosome,Start,End\n";


	# Get iteratin'!
	foreach my $gene (sort keys %GenesToDirNames)
	{

		my $gene_dir_name = $GenesToDirNames{$gene};

		my $summary_file_name = $gene_dir_name.'summary.out';
		open(my $SummaryFile,'<',$summary_file_name) 
			|| die "\n  ERROR:  Failed to open gene search summary file '$summary_file_name'\n\n";

		while (my $line = <$SummaryFile>)
		{

			next if ($line !~ /Input ID\s+:\s+(\S+)/);
			my $input_id = $1;


			<$SummaryFile>; # BATH   runtime
			<$SummaryFile>; # splash runtime
			<$SummaryFile>; # model length


			# We'll want to know how many exon sets there are
			# for this query
			$line = <$SummaryFile>;
			$line =~ /Num Exon Sets\s+:\s+(\d+)/;
			my $num_exon_sets = $1;


			for (my $exon_set_id = 1; $exon_set_id <= $num_exon_sets; $exon_set_id++)
			{

				<$SummaryFile>; # Exon Set ID


				# How many exons in this set?
				$line = <$SummaryFile>;
				$line =~ /Num Exons\s+:\s+(\d+)/;
				my $num_exons = $1;


				<$SummaryFile>; # Model range


				# What's the chromosome (and range within that chromosome)
				$line = <$SummaryFile>;
				$line =~ /Chromosome\s+:\s+([^\/]+)\/(\d+)-(\d+)/;
				my $chr       = $1;
				my $chr_start = $2;
				my $chr_end   = $3;

				my $chr_revcomp = 0;
				$chr_revcomp = 1 if ($chr_start > $chr_end);


				<$SummaryFile>; # Range of exon set within the chromosomal sub-region
				<$SummaryFile>; # Missed exons


				for (my $exon_id = 1; $exon_id <= $num_exons; $exon_id++)
				{

					<$SummaryFile>; # Exon ID
					<$SummaryFile>; # Model Range


					# Nucleotide range
					$line = <$SummaryFile>;
					$line =~ /Nucl\. Range\s+:\s+(\d+)\.\.(\d+)/;
					my $exon_start = $1;
					my $exon_end   = $2;


					# Correct it!
					if ($chr_revcomp)
					{
						$exon_start = $chr_start - $exon_start + 1;
						$exon_end   = $chr_start - $exon_end   + 1;
					}
					else
					{
						$exon_start = $chr_start + $exon_start - 1;
						$exon_end   = $chr_start + $exon_end   - 1;
					}


					# Become yelling!
					print $OutFile "$gene,$input_id,$exon_set_id:$exon_id/$num_exons,$chr,$exon_start,$exon_end\n";


					# Eat the splice signal
					<$SummaryFile>;

				}

			}

		}

		close($SummaryFile);

	}

	close($OutFile);


	return $splash_coords_file_name;

}

