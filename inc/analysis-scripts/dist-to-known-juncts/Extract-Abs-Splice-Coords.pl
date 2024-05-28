#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV != 1) { die "\n  USAGE:  ./Extract-Abs-Splice-Coords.pl [Splash-Validation-Output/]\n\n"; }


my $in_dir_name = $ARGV[0];
$in_dir_name = $in_dir_name.'/' if ($in_dir_name !~ /\/$/);
die "\n  ERROR:  Failed to locate input directory '$in_dir_name'\n\n"
	if (!(-d $in_dir_name));

my $rbg_dir_name = $in_dir_name.'Results-by-Gene/';
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


my $out_file_name = 'splice-coords.csv';
open(my $OutFile,'>',$out_file_name)
	|| die "\n  ERROR:  Failed to open output file '$out_file_name'\n\n";


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


1;
