#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub GetCoveredSeqs;
sub ContainsFullSingles;


if (@ARGV != 1) { die "\n  USAGE:  ./Count-Full-Model-Singles.pl [Splash-PANTHER-Output/]\n\n"; }


my $result_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate input directory '$result_dir_name'\n\n"
	if (!(-d $result_dir_name));
$result_dir_name = $result_dir_name.'/' if ($result_dir_name !~ /\/$/);


my $covered_seqs_ref = GetCoveredSeqs($result_dir_name);
my %CoveredSeqs = %{$covered_seqs_ref}; # these are in the form of {species|fam}


#my @CoveredSeqKeys = keys %CoveredSeqs;
#for (my $x=0; $x<10; $x++)
#{
#	print "$CoveredSeqKeys[$x]\n";
#}
#my $sanity_check = scalar(@CoveredSeqKeys);
#print "\n  7796 * 4 =!=  $sanity_check\n\n";


my $rbg_dir_name = $result_dir_name.'Results-by-Gene/';
die "\n  ERROR:  Failed to locate 'Results-by-Gene' directory '$rbg_dir_name'\n\n"
	if (!(-d $rbg_dir_name));


opendir(my $RBG, $rbg_dir_name)
	|| die "\n  ERROR:  Failed to open 'Results-by-Gene' directory '$rbg_dir_name'\n\n";

my %SpeciesToFullSingles;
my %SpeciesToRecoverySingles;

while (my $fam = readdir($RBG))
{
	next if ($fam =~ /^\./);
	$fam =~ s/\/$//;

	my $fam_dir_name = $rbg_dir_name.$fam.'/';
	next if (!(-d $fam_dir_name));

	opendir(my $FamDir, $fam_dir_name)
		|| die "\n  ERROR:  Failed to open family directory '$fam_dir_name'\n\n";

	while (my $file_name = readdir($FamDir))
	{
		next if ($file_name !~ /^$fam\.([^\.]+)\.out$/);
		my $species = lc($1);

		$file_name = $fam_dir_name.$file_name;

		if (ContainsFullSingles($file_name))
		{
			if ($SpeciesToFullSingles{$species}) { $SpeciesToFullSingles{$species}++; }
			else                                 { $SpeciesToFullSingles{$species}=1; }

			if (!$CoveredSeqs{$species.'|'.$fam})
			{
				if ($SpeciesToRecoverySingles{$species}) { $SpeciesToRecoverySingles{$species}++; }
				else                                     { $SpeciesToRecoverySingles{$species}=1; }
			}
			#else
			#{
			#	print "$species\|$fam\n";
			#}

		}

	}

	closedir($FamDir);

}

closedir($RBG);


my $longest_species_name_len = 0;
foreach my $species (sort keys %SpeciesToFullSingles)
{
	if (length($species) > $longest_species_name_len)
	{
		$longest_species_name_len = length($species);
	}
}


foreach my $species (sort keys %SpeciesToFullSingles)
{
	my $num_full_singles = $SpeciesToFullSingles{$species};
	my $recovery_singles = 0;
	$recovery_singles = $SpeciesToRecoverySingles{$species}
		if ($SpeciesToRecoverySingles{$species});

	while (length($species) <= $longest_species_name_len)
	{
		$species = $species.' ';
	}

	print "  $species : $num_full_singles ($recovery_singles had no spliced coverage)\n";

}



1;







##########################################################################
#
#  GetCoveredSeqs
#
sub GetCoveredSeqs
{
	my $result_dir_name = shift;

	opendir(my $ResultDir, $result_dir_name)
		|| die "\n  ERROR:  Failed to open result directory '$result_dir_name'\n\n";

	my %CoveredSeqs;

	while (my $file_name = readdir($ResultDir))
	{
		if ($file_name =~ /^(\S+)-comparison\.csv$/)
		{
			my $species = $1;

			open(my $CompFile,'<',$result_dir_name.$file_name)
				|| die "\n  ERROR:  Failed to open comparison file '$result_dir_name$file_name'\n\n";

			<$CompFile>; # Header

			while (my $line = <$CompFile>)
			{
				$line =~ s/\n|\r//g;
				next if (!$line);

				my @Data = split(/,/,$line);
				if ($Data[4] ne '-')
				{
					$CoveredSeqs{$species.'|'.$Data[0]} = 1;
				}

			}

			close($CompFile);

		}
		elsif ($file_name =~ /(\S+)-Coverage\.csv$/)
		{
			my $species = $1;

			open(my $CoverFile,'<',$result_dir_name.$file_name)
				|| die "\n  ERROR:  Failed to open coverage file '$result_dir_name$file_name'\n\n";

			<$CoverFile>; # Header
			
			while (my $line = <$CoverFile>)
			{
				$line =~ s/\n|\r//g;
				next if (!$line);

				$line =~ /^([^,]+),/;
				$CoveredSeqs{$species.'|'.$1} = 1;
			}

			close($CoverFile);

		}
	}

	closedir($ResultDir);

	return \%CoveredSeqs;

}






##########################################################################
#
#  ContainsFullSingles
#
sub ContainsFullSingles
{
	my $file_name = shift;


	open(my $File,'<',$file_name)
		|| die "\n  ERROR:  Failed to open output file '$file_name'\n\n";


	# Fast-forward past any splashery
	my $model_length = 0;
	while (my $line = <$File>)
	{
		if ($line =~ /Query:\s+\S+\s+\[M=(\d+)\]/)
		{
			$model_length = $1;
		}
		last if ($line =~ /Scores for complete hits:/);
	}


	# UH OH
	if ($model_length == 0)
	{
		close($File);
		return 0;
	}


	while (!eof($File))
	{
		my $line = <$File>;
		if ($line =~ /^\s+score\s+bias\s+Evalue\s+hmm-from\s+hmm-to/)
		{
			<$File>; # underlines

			my $overview = <$File>;
			if ($overview =~ /^\s*\S\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)/)
			{
				my $hmm_from = $1;
				my $hmm_to   = $2;

				if ($hmm_from - 6 <= 0 && $hmm_to + 6 > $model_length) 
				{
					close($File);
					return 1;
				}
			}
			else
			{
				die "\n$file_name\n$overview\n";
			}

		}
	}

	close($File);

	return 0;

}


