#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub ContainsFullSingles;


if (@ARGV != 1) { die "\n  USAGE:  ./Count-Full-Model-Singles.pl [Splash-PANTHER-Output/]\n\n"; }


my $result_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate input directory '$result_dir_name'\n\n"
	if (!(-d $result_dir_name));
$result_dir_name = $result_dir_name.'/' if ($result_dir_name !~ /\/$/);


my $rbg_dir_name = $result_dir_name.'Results-by-Gene/';
die "\n  ERROR:  Failed to locate 'Results-by-Gene' directory '$rbg_dir_name'\n\n"
	if (!(-d $rbg_dir_name));


opendir(my $RBG, $rbg_dir_name)
	|| die "\n  ERROR:  Failed to open 'Results-by-Gene' directory '$rbg_dir_name'\n\n";

my %SpeciesToFullSingles;

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
	while (length($species) <= $longest_species_name_len)
	{
		$species = $species.' ';
	}
	print "  $species : $num_full_singles\n";
}



1;





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


