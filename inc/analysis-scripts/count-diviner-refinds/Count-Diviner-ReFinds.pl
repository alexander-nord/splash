#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub ForwardOverlap;
sub ReverseOverlap;



if (@ARGV != 2) { die "\n  USAGE: ./Count-Diviner-ReFinds.pl [Splash-Diviner-Output/] [Diviner-Results/]\n\n"; }



my $spl_results_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate Splash results directory '$spl_results_dir_name'\n\n"
	if (!(-d $spl_results_dir_name));

$spl_results_dir_name = $spl_results_dir_name.'/' if ($spl_results_dir_name !~ /\/$/);


my $spl_rbg_dir_name = $spl_results_dir_name.'Results-by-Gene/';
die "\n  ERROR:  Failed to locate Splash 'Results-by-Gene' directory '$spl_rbg_dir_name'\n\n"
	if (!(-d $spl_rbg_dir_name));



my $div_dir_name = $ARGV[1];
die "\n  ERROR:  Failed to locate Diviner results directory '$div_dir_name'\n\n"
	if (!(-d $div_dir_name));
$div_dir_name = $div_dir_name.'/' if ($div_dir_name !~ /\/$/);

my $div_rbg_dir_name = $div_dir_name.'Results-by-Gene/';
die "\n  ERROR:  Failed to locate Diviner 'Results-by-Gene' directory '$div_rbg_dir_name'\n\n"
	if (!(-d $div_rbg_dir_name));



opendir(my $SplRBG,$spl_rbg_dir_name)
	|| die "\n  ERROR:  Failed to open input directory '$spl_rbg_dir_name'\n\n";

my $num_hmm_matches  = 0;
my $num_nucl_matches = 0;
my $num_misses       = 0;
my $num_range_checks = 0;

while (my $gene = readdir($SplRBG))
{

	$gene =~ s/\/$//;
	next if ($gene =~ /^\./);


	my $gene_spl_dir_name = $spl_rbg_dir_name.$gene.'/';
	die "\n  ERROR:  Missing gene directory '$gene_spl_dir_name'\n\n"
		if (!(-d $gene_spl_dir_name));

	my $gene_div_dir_name = $div_rbg_dir_name.$gene.'/';
	die "\n  ERROR:  Missing gene directory '$gene_div_dir_name'\n\n"
		if (!(-d $gene_div_dir_name));
	

	# Even though this isn't the most efficient approach,
	# we're going to a hit-by-hit independent check of the
	# novel Diviner hits to Splash alignments.

	my $div_hs_file_name = $gene_div_dir_name.'alignments/homo_sapiens.'.$gene.'.out';
	open(my $DivHSFile,'<',$div_hs_file_name)
		|| die "\n  ERROR:  Failed to open Diviner results file '$div_hs_file_name'\n\n";

	while (my $line = <$DivHSFile>)
	{

		next if ($line !~ /Target\s+:\s+homo_sapiens\s+(\S+)/);
		my $div_chr_range = $1;

		$line = <$DivHSFile>;
		next if ($line !~ /Novel exon/);

		<$DivHSFile>; # Exon identifier

		$line = <$DivHSFile>;
		
		my %SpeciesToAminoRange;
		while ($line =~ /:\s+(\S+)\s+\/\s+aminos\s+(\S+)/) 
		{
			$SpeciesToAminoRange{$1} = $2;
			$line = <$DivHSFile>;
		}


		my $div_hmm_range = GetHMMPosRange($gene_div_dir_name.'species.afa',\%SpeciesToAminoRange);


		# Check if we have an overlap in the HMM range.
		# AND... does it overlap with the Diviner-labelled genomic window?
		my ($hmm_range_overlap,$nucl_range_overlap)
			= CheckHitOverlaps($gene_spl_dir_name.$gene.'.homo_sapiens.out',$div_chr_range,$div_hmm_range);

		$num_range_checks++;

			
		print "$gene ($div_hmm_range/$div_chr_range): ";

		if ($hmm_range_overlap)
		{
			$num_hmm_matches++;
			print "HMM Match";

			if ($nucl_range_overlap)
			{
				$num_nucl_matches++;
				print " and Genome Overlap ($nucl_range_overlap)";
			}

			print "\n";

		}
		else
		{
			print "missed\n";
			$num_misses++;
		}

	}

	close($DivHSFile);

}

closedir($SplRBG);



print "\n";
print "  Summary Stats\n";
print "  - Total Num Checks : $num_range_checks\n";
print "  - Num HMM  Matches : $num_hmm_matches\n";
print "  - Num Nucl Matches : $num_nucl_matches\n";
print "  - Total Misses     : $num_misses\n";
print "\n";



1;






sub GetHMMPosRange
{
	my $ali_file_name = shift;
	my $star_ref = shift;

	my %SpeciesToAminoRange = %{$star_ref};

	my @InitMSA;
	my @MSASpecies;
	my $init_msa_len;
	my $num_species = -1;

	open(my $AliFile,'<',$ali_file_name)
		|| die "\n  ERROR:  Failed to open MSA file '$ali_file_name'\n\n";

	my $valid_species = 0;
	while (my $line = <$AliFile>)
	{
		$line =~ s/\n|\r//g;
		next if (!$line);

		if ($line =~ /^\>(\S+)/)
		{
			my $next_species = $1;

			if ($SpeciesToAminoRange{$next_species})
			{
				$valid_species = 1;
				$init_msa_len  = 0;
				$num_species++;
				$MSASpecies[$num_species] = $next_species;
			}
			else
			{
				$valid_species = 0;
			}
		}
		elsif ($valid_species)
		{
			foreach my $char (split(//,uc($line)))
			{
				$InitMSA[$num_species][$init_msa_len] = $char;
				$init_msa_len++;
			}
		}
	}

	close($AliFile);


	my @MSA;
	my $msa_len = 0;
	for (my $init_pos = 0; $init_pos < $init_msa_len; $init_pos++)
	{
		my $aminoful = 0;
		for (my $species_id = 0; $species_id < $num_species; $species_id++)
		{
			my $char = $InitMSA[$species_id][$init_pos];
			$MSA[$species_id][$msa_len] = $char;
			$aminoful = 1 if ($char =~ /[A-Z]/);
		}
		$msa_len += $aminoful;
	}


	
	my $hmm_range_start;
	my $hmm_range_end;

	for (my $species_id = 0; $species_id < $num_species; $species_id++)
	{

		my $species = $MSASpecies[$species_id];

		$SpeciesToAminoRange{$species} =~ /^(\d+)\.\.(\d+)$/;
		my $start_amino = int($1);
		my $end_amino   = int($2);

		my $start_pos;
		my $end_pos=0;

		my $msa_pos   = 0;
		my $amino_pos = 0;
		while ($msa_pos < $msa_len && !$end_pos)
		{
			$amino_pos++ if ($MSA[$species_id][$msa_pos] =~ /[A-Z]/);

			$start_pos = $msa_pos if ($amino_pos == $start_amino);
			$end_pos   = $msa_pos if ($amino_pos == $end_amino  );

			$msa_pos++;
		}

		if (!$species_id || $hmm_range_start > $start_pos)
		{
			$hmm_range_start = $start_pos;
		}
		if (!$species_id || $hmm_range_end < $end_pos)
		{
			$hmm_range_end = $end_pos;
		}

	}


	return $hmm_range_start.'..'.$hmm_range_end;

}






sub CheckHitOverlaps
{
	my $spl_file_name = shift;
	my $div_chr_range = shift;
	my $div_hmm_range = shift;


	$div_chr_range =~ /^([^:]+):(\S+)$/;
	my $div_chr   = $1;
	my $div_range = $2;

	my $revcomp = 0;
	if ($div_chr =~ /\[revcomp\]/)
	{
		$div_chr =~ s/\[revcomp\]//;
		$revcomp = 1;
	}


	open(my $SplFile,'<',$spl_file_name)
		|| die "\n  ERROR:  Failed to open Splash output file '$spl_file_name'\n\n";

	my $hmm_range_overlap = 0;
	while (my $line = <$SplFile>)
	{

		next if ($line !~ /\| = Exon Set \d+ \((\d+)/);
		my $num_exons = $1;

		$line = <$SplFile>;
		$line =~ /\| = Model Positions (\S+)/;
		my $spl_hmm_range = $1;


		$line = <$SplFile>;
		$line =~ /\| = Target Seq Name (\S+)/;
		my $spl_chr = $1;

		next if ($spl_chr ne $div_chr);


		$line = <$SplFile>;
		$line =~ /\| = Nucleotide Coords (\d+)\.\.(\d+)/;

		my $spl_nucl_start = $1;
		my $spl_nucl_end   = $2;

		if ($revcomp)
		{
			next if ($spl_nucl_start < $spl_nucl_end);
		}
		else
		{
			next if ($spl_nucl_start > $spl_nucl_end);
		}



		# We're on the right chromosome facing the right direction!
		# Checks are valid!

		
		$hmm_range_overlap = ForwardOverlap($spl_hmm_range,$div_hmm_range) 
			if (!$hmm_range_overlap);

		
		for (my $exon_id = 1; $exon_id <= $num_exons; $exon_id++)
		{
			$line = <$SplFile>;
			$line =~ /\| = Exon $exon_id: \S+ \/ (\S+)/;
			my $exon_nucl_range = $1;

			if ($revcomp)
			{
				if (ReverseOverlap($div_range,$exon_nucl_range))
				{
					close($SplFile);
					return (1,$spl_chr.':'.$exon_nucl_range);
				}
			}
			elsif (ForwardOverlap($div_range,$exon_nucl_range))
			{
				close($SplFile);
				return (1,$spl_chr.':'.$exon_nucl_range);
			}
		}

	}

	close($SplFile);


	return ($hmm_range_overlap,0);

}






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




