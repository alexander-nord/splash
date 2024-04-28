#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub HelpAndDie;
sub ParseCommandArguments;
sub ConfirmGenomeDir;
sub ConfirmGenome;
sub ConfirmRequiredTools;




# We'll make sure we have the necessary tools before proceeding
my $SFETCH;
my $SEQSTAT;
my $HMMSEARCHT;
ConfirmRequiredTools();


# If the user didn't provide any arguments, help 'em out
HelpAndDie() if (@ARGV == 0);


# What do they want us to even do?
my %OPTIONS;
ParseCommandArguments();


# Alright, let's figure this stuff out!




1;









###########################################################################
#
#  Subroutine: HelpAndDie
#
sub HelpAndDie
{
	print "\n";
	print "  Use Case 1:  Search all sequences for a gene family\n";
	print "            :\n";
	print "            :  ./Splash.pl {OPT.S} [gene]\n";
	print "            '---------------------------------------------\n";
	print "\n";
	print "  Use Case 2:  Search using a specific file as input\n";
	print "            :\n";
	print "            :  ./Splash.pl {OPT.S} [file.fa]\n";
	print "            '---------------------------------------------\n";
	print "\n";
	print "\n";
	print "  OPT.S: --full-genome : Run test on full genome, rather than\n";
	print "                           extracting a relatively narrow region\n";
	print "                           of sequence associated with a gene's\n";
	print "                           coding region.\n";
	print "\n";
	die   "\n";
}









###########################################################################
#
#  Subroutine: ParseCommandArguments
#
sub ParseCommandArguments
{


	# Unless specified, we're looking for data in this directory
	$OPTIONS{'inputs-dir'} = '../inputs-to-splash/';
	$OPTIONS{'genome-dir'} = '../inputs-to-splash/genomes/';


	my $num_args = scalar(@ARGV);
	for (my $arg_id=0; $arg_id<$num_args; $arg_id++) {

		my $Arg = $ARGV[$arg_id];

		if (lc($Arg) =~ '^-?-?full-genome$') 
		{
			$OPTIONS{'full-genome'} = 1;
		}
		elsif ($arg_id == $num_args-1 && !$OPTIONS{'protein-input'})
		{
			if (-e $Arg) 
			{
				
				# We're being directed to use a specific file
				$OPTIONS{'protein-input'} = $Arg;
				$OPTIONS{'single-seq'} = 1;
				$OPTIONS{'meta-dir'}   = 0;
				$OPTIONS{'family'}     = 0;

				# Is there a genome guide? If not, assume full genome
				$Arg =~ /^(.*\/?[^\/]+)\.[^\.|\/]+$/;	
				my $seq_no_ext = $1;

				$Arg =~ /^(.*\/?[^\.|\/]+)[^\/]+$/;
				my $seq_bare_base = $1;
				
				if (-e $seq_no_ext.'.full-range.out') 
				{
					$OPTIONS{'single-seq-range'} = $seq_no_ext.'.full-range.out';
				} 
				elsif (-e $seq_bare_base) 
				{
					$OPTIONS{'single-seq-range'} = $seq_bare_base.'.full-range.out';
				}
				else
				{
					$OPTIONS{'full-genome'} = 1;
				}

			}
			elsif (-d $Arg)
			{
				# For now, we're going to assume this is a directory full of
				# sub-directories specific to protein families (or maybe some
				# other configuration), but for now we'll just call it a meta-dir
				# and leave the work of understanding what that means for later
				$OPTIONS{'protein-input'} = $Arg;
				$OPTIONS{'meta-dir'}   = 1;
				$OPTIONS{'family'}     = 0;
				$OPTIONS{'single-seq'} = 0;
			}
			elsif (-d $OPTIONS{'inputs-dir'}.'protein-data/'.lc($Arg)) 
			{
				# We're working with a family in our default protein dataset
				$OPTIONS{'protein-input'} = $OPTIONS{'inputs-dir'}.'protein-data/'.lc($Arg).'/';
				$OPTIONS{'family'}     = 1;
				$OPTIONS{'meta-dir'}   = 0;
				$OPTIONS{'single-seq'} = 0;
			}
			else
			{
				die "\n  ERROR:  Failed to make sense of (inferred) query instruction '$Arg'\n\n";
			}
		}
		else
		{
			die "\n  ERROR:  Unrecognized (assumed) protein input argument: '$Arg'\n\n";
		}

	}


	# If we don't have mapping instructions for associating species
	# with genomes, generate them!
	if (!$OPTIONS{'genome-map'}) 
	{
		$OPTIONS{'genome-map'} = ConfirmGenomeDir($OPTIONS{'genome-dir'});
	}



}









###########################################################################
#
#  Subroutine: ConfirmGenomeDir
#
sub ConfirmGenomeDir
{
	
	my $genome_dir_name   = shift;
	my $species_to_genome = $genome_dir_name.'species-to-genome';

	if (!(-e $species_to_genome)) {
		die "\n  ERROR:  Without specified genomes, species-to-genome file is required in genome directory (failed to find '$species_to_genome')\n\n";
	}
	open(my $SpeciesToGenome,'<',$species_to_genome)
		|| die "\n  ERROR:  Failed to open species-to-genome file '$species_to_genome'\n\n";
	
	my @SpeciesGenomeMaps;
	while (my $line = <$SpeciesToGenome>) {
		if ($line =~ /^\s*(\S+)\s*(\S+)\s*$/) {
			my $species = $1;
			my $genome  = ConfirmGenome($genome_dir_name.$2);
			push(@SpeciesGenomeMaps,$species.'|'.$genome);
		}
	}
	close($SpeciesToGenome);


	my $species_genome_map_str = '';
	for (my $i=0; $i<scalar(@SpeciesGenomeMaps); $i++) {
		$species_genome_map_str = $species_genome_map_str.'&' if ($i);
		$species_genome_map_str = $species_genome_map_str.$SpeciesGenomeMaps[$i];
	}

	return $species_genome_map_str;

}









###########################################################################
#
#  Subroutine: ConfirmGenome
#
sub ConfirmGenome
{
	my $genome = shift;

	if (!(-e $genome)) {
		die "\n  ERROR:  Failed to locate genomic file '$genome'\n\n";
	}

	my $ssi_file_name = $genome.'.ssi';
	if (!(-e $ssi_file_name) && system("$SFETCH --index \"$genome\"")) {
		die "\n  ERROR:  Failed to produce an ssi index for genomic file '$genome'\n\n";
	}

	my $seqstat_file_name = $genome.'.seqstat.out';
	if (!(-e $seqstat_file_name) && system("$SEQSTAT -a \"$genome\" > \"$seqstat_file_name\"")) {
		die "\n  ERROR:  Failed to aggregate sequence statistics for genomic file '$genome'\n\n";
	}

	return $genome; # Looks like a genome to me!

}








###########################################################################
#
#  Subroutine: ConfirmRequiredTools
#
sub ConfirmRequiredTools
{

	my $th_base_dir = 'trans-hmmer/';
	die "\n  ERROR:  Failed to locate base translated HMMER3 directory '$th_base_dir'\n          Try running './Build.pl'\n\n"
		if (!(-d $th_base_dir));


	my $easel_dir = $th_base_dir.'easel/';
	die "\n  ERROR:  Failed to locate easel...?  ('$easel_dir' should exist)\n\n"
		if (!(-d $easel_dir));

	my $miniapps_dir = $easel_dir.'miniapps/';
	die "\n  ERROR:  Failed to locate easel/miniapps...? ('$miniapps_dir' should exist)\n\n"
		if (!(-d $miniapps_dir));

	my $SFETCH = $miniapps_dir.'esl-sfetch';
	die "\n  ERROR:  Failed to locate required Easel tool sfetch (looking for '$SFETCH')\n\n"
		if (!(-x $SFETCH));

	my $SEQSTAT = $miniapps_dir.'esl-seqstat';
	die "\n  ERROR:  Failed to locate required Easel tool seqstat (looking for '$SEQSTAT')\n\n"
		if (!(-x $SEQSTAT));


	my $th_src_dir = $th_base_dir.'src/';
	die "\n  ERROR:  Failed to locate translated HMMER3 source directory '$th_src_dir'\n\n"
		if (!(-d $th_src_dir));

	$HMMSEARCHT = $th_src_dir.'hmmsearcht';
	die "\n  ERROR:  Failed to locate hmmsearcht in HMMER3 source directory (looking for '$HMMSEARCHT')\n\n"
		if (!(-x $HMMSEARCHT));


}


