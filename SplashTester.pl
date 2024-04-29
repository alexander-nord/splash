#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub FamilySplash;
sub DetermineTargetSeq;
sub GenomeRangeFileToTargetFile;
sub HelpAndDie;
sub ParseCommandArguments;
sub ConfirmGenomeDir;
sub ConfirmGenome;
sub ConfirmRequiredTools;




# We'll make sure we have the necessary tools before proceeding
my $SFETCH;
my $SEQSTAT;
my $HMMBUILD;
my $HMMSEARCHT;
ConfirmRequiredTools();


# If the user didn't provide any arguments, help 'em out
HelpAndDie() if (@ARGV == 0);


# What do they want us to even do?
my %OPTIONS;
my %SPECIES_TO_GENOME;
my %GENOME_LIST; # So that we don't delete anything important
ParseCommandArguments();


# We survived commandline parsing, so let's build an
# output directory!
die "\n  ERROR:  Failed to create output directory '$OPTIONS{'output-dir'}'\n\n"
	if (system("mkdir \"$OPTIONS{'output-dir'}\""));

# Just in case...
my $ERROR_FILE = $OPTIONS{'output-dir'}.'Splash.err';


# Alright, let's figure this stuff out!
if    ($OPTIONS{'meta-dir'}  ) { die "\n  Not yet supported\n\n";         }
elsif ($OPTIONS{'family'}    ) { FamilySplash($OPTIONS{'protein-input'}); }
elsif ($OPTIONS{'single-seq'}) { die "\n  Not yet supported\n\n";         }
else {   die "\n  ERROR:  Failed to recognize protein input\n\n";         }



1;










###########################################################################
#
#  Subroutine: FamilySplash
#
sub FamilySplash
{

	my $family_dir_name = shift;


	# This is, perhaps, presumptuous...
	$family_dir_name =~ /\/([^\/]+)\/$/;
	my $gene = $1;


	# Read through the family's directory, figuring
	# out what the situation is with our query data
	opendir(my $FamilyDir,$family_dir_name) 
		|| die "\n  ERROR:  Failed to open directory '$family_dir_name' (FamilySplash)\n\n";
	my @InputHMMs;
	my @FilesToHMMBUILD;
	while (my $file_name = readdir($FamilyDir)) 
	{

		next if ($file_name =~ /^\./);

		# Have we stumbled upon a real-life HMM?!
		if (lc($file_name) =~ /\.hmm$/) 
		{
			push(@InputHMMs,$family_dir_name.$file_name);
		}
		elsif (lc($file_name) =~ /^(\S+)\.a?fa[sta]?$/)
		{
			my $file_base_name = $1;
			push(@FilesToHMMBUILD,$family_dir_name.$file_name) 
				if (!(-e $file_base_name.'.hmm'));
		}

	}
	closedir($FamilyDir);



	# Build HMMs for any sequence files that didn't have one
	foreach my $file_to_hmmbuild (@FilesToHMMBUILD)
	{
		
		$file_to_hmmbuild =~ /^(\S+)\.[^\.]+$/;
		my $hmm_file_name = $1.'.hmm';
		my $hmmbuild_cmd  = "$HMMBUILD \"$hmm_file_name\" \"$file_to_hmmbuild\"";

		if (system($hmmbuild_cmd)) {
			die "\n  ERROR:  Failed to build HMM on '$file_to_hmmbuild' (command:'$hmmbuild_cmd')\n\n";
		}

		push(@InputHMMs,$hmmbuild_cmd);

	}



	# If we didn't find any HMMs (or HMM-ables) shout, scream, and cry
	if (scalar(@InputHMMs) == 0) 
	{
		my $err_message = "Failed to locate any HMM or sequence files in directory '$family_dir_name'\n";
		system("echo \"$err_message\" >> $ERROR_FILE");
		print "  WARNING: $err_message\n";
		return;
	}



	# Make an output directory for this family
	my $fam_out_dir_name = $OPTIONS{'output-dir'}.$gene.'/';
	if (system("mkdir \"$fam_out_dir_name\"")) {
		die "\n  ERROR:  Failed to create output directory for family '$gene' (attempted name: '$fam_out_dir_name')\n\n";
	}


	# Iterate over each HMM, determining the appropriate
	# target sequence and running hmmsearcht!
	my $num_fam_errors    = 0;
	my $num_fam_victories = 0;
	foreach my $hmm_file_name (@InputHMMs) 
	{

		my $target_file_name = DetermineTargetSeq($hmm_file_name);
		

		$hmm_file_name =~ /\/([^\/]+)\.hmm$/;
		my $out_file_name = $fam_out_dir_name.$1.'.out';
		my $err_file_name = $fam_out_dir_name.$1.'.err';

		my $hmmsearcht_cmd = "$HMMSEARCHT -o $out_file_name $hmm_file_name $target_file_name 2>$err_file_name";
		if (system($hmmsearcht_cmd)) 
		{
			$num_fam_errors++;
			system("echo \"\% $hmmsearcht_cmd\" >> $ERROR_FILE");
			system("cat $err_file_name >> $ERROR_FILE");
		} 
		else 
		{
			$num_fam_victories++;
		}


		# If we created a target sequence file, kill it!
		system("rm \"$target_file_name\"") unless ($GENOME_LIST{$target_file_name});

	}


}










###########################################################################
#
#  Subroutine: DetermineTargetSeq
#
sub DetermineTargetSeq
{

	my $hmm_file_name = shift;

	# Any chance there's a file related to this one with
	# genome range data to pull in?
	$hmm_file_name =~ /^(.*\/?)([^\/]+)\.hmm$/;
	my $hmm_file_dir_name  = $1;
	my $hmm_file_base_name = $2;


	# If they used 'MirageOutToSplashIn' then we can find out
	# the species easily.
	$hmm_file_base_name =~ /^([^\.]+)\./;
	my $presumptive_species = lc($1);
	if ($OPTIONS{'full-genome'}) 
	{
		if ($SPECIES_TO_GENOME{$presumptive_species}) { 
			return $SPECIES_TO_GENOME{$presumptive_species}; 
		} else {
			die "\n  ERROR:  Failed to locate genome for (presumptive) species '$presumptive_species' (for file '$hmm_file_name')\n\n";
		}
	}


	# Is there a 'genome-ranges.out' file associated with this query?
	my $range_file_name = 0;
	if (-e $hmm_file_dir_name.$hmm_file_base_name.'.genome-range.out')
	{
		$range_file_name = $hmm_file_dir_name.$hmm_file_base_name.'.genome-range.out';
	} 
	elsif (-e $hmm_file_dir_name.$presumptive_species.'.genome_range.out')
	{
		$range_file_name = $hmm_file_dir_name.$presumptive_species.'.genome_range.out';
	}


	# Are we operating on a range file?
	# If not, you're using the full genome, baby!
	if ($range_file_name)
	{
		my $target_file_name = $hmm_file_name;
		$target_file_name =~ s/\.hmm$/\.target\.fa/;

		return GenomeRangeFileToTargetFile($range_file_name,$target_file_name);

	}
	elsif ($SPECIES_TO_GENOME{$presumptive_species})
	{
		return $SPECIES_TO_GENOME{$presumptive_species};
	}

	die "\n  ERROR:  Failed to determine target sequence for query '$hmm_file_name'\n\n";

}










###########################################################################
#
#  Subroutine: GenomeRangeFileToTargetFile
#
sub GenomeRangeFileToTargetFile
{

	my $range_file_name  = shift;
	my $target_file_name = shift;

	open(my $RangeFile,'<',$range_file_name);
	my %RangeFileData;
	while (my $line = <$RangeFile>) {
		if ($line =~ /^\s*(\S+)\s*:\s*(\S+)\s*$/) {
			$RangeFileData{lc($1)} = $2;
		}
	}
	close($RangeFile);



	if (!$RangeFileData{'species'}) {
		die "\n  ERROR:  Range file '$range_file_name' does not have recognizable species entry\n\n";
	}
	my $species = $RangeFileData{'species'};


	if (!$SPECIES_TO_GENOME{$species}) {
		die "\n  ERROR:  Species '$species' has no listed genome (looking at '$range_file_name')\n\n";
	}
	my $genome = $SPECIES_TO_GENOME{$species};


	if (!$RangeFileData{'chr'}) {
		die "\n  ERROR:  Range file '$range_file_name' does not have recognizable chromosome entry\n\n";
	}
	my $chr = $RangeFileData{'chr'};


	# We'll allow a single sequence, I *suppose*
	#if (!$RangeFileData{'range'}) {
	#	die "\n  ERROR:  Range file '$range_file_name' does not have recognizable range entry\n\n";
	#}
	#my $range = $RangeFileData{'range'};


	# Put together the sfetch command
	my $sfetch_cmd = "$SFETCH -o $target_file_name";
	$sfetch_cmd = $sfetch_cmd." -c $RangeFileData{'range'}" if ($RangeFileData{'range'});
	$sfetch_cmd = $sfetch_cmd." $genome $chr";


	# Clear the way!
	system("rm \"$target_file_name\"") if (-e $target_file_name);

	if (system($sfetch_cmd)) {
		die "\n  ERROR:  Subsequence range extraction failed (command:'$sfetch_cmd')\n\n";
	}
	return $target_file_name;

}









###########################################################################
#
#  Subroutine: HelpAndDie
#
sub HelpAndDie
{
	print "\n";
	print "\n";
	print "  Use Case 1:  Search all sequences for a gene family\n";
	print "            :\n";
	print "            :  ./SplashTester.pl {OPT.S} [gene]\n";
	print "            '---------------------------------------------\n";
	print "\n";
	print "  Use Case 2:  Search using a specific file as input\n";
	print "            :\n";
	print "            :  ./SplashTester.pl {OPT.S} [file.fa]\n";
	print "            '---------------------------------------------\n";
	print "\n";
	print "\n";
	print "  OPT.S: --full-genome : Force use of full genome as target sequence.\n";
	print "\n";
	die   "\n";
}









###########################################################################
#
#  Subroutine: ParseCommandArguments
#
sub ParseCommandArguments
{

	# Default output directory name
	my $default_out_dir_name = 'Splash-Results';
	my $out_dir_name = $default_out_dir_name;
	my $attempt = 1;
	while (-d $out_dir_name) {
		$attempt++;
		$out_dir_name = $default_out_dir_name.'-'.$attempt;
	}
	$OPTIONS{'output-dir'} = $out_dir_name.'/';


	# Unless specified, we're looking for data in this directory
	$OPTIONS{'inputs-dir'} = '../inputs-to-splash/';
	$OPTIONS{'genome-dir'} = '../inputs-to-splash/genomes/';


	my $num_args = scalar(@ARGV);
	for (my $arg_id=0; $arg_id<$num_args; $arg_id++) {

		my $Arg = $ARGV[$arg_id];

		if (lc($Arg) =~ '^-?-?full-genome$') 
		{
			$OPTIONS{'full-genome'} = 1; # Overrides existence of '.genome-range.out' files
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
	BuildSpeciesToGenomeMap();



}







###########################################################################
#
#  Subroutine: BuildSpeciesToGenomeMap
#
sub BuildSpeciesToGenomeMap
{

	# Convert the species-to-genome mapping string to a hash
	foreach my $species_genome_pair (split(/\&/,$OPTIONS{'genome-map'})) 
	{

		$species_genome_pair =~ /^([^\|]+)\|(\S+)$/;
		my $species = $1;
		my $genome  = $2;

		if (!(-e $genome)) 
		{
			die "\n  ERROR:  Failed to locate genome file '$genome' (associated with species '$species')\n\n";
		}
		elsif (!(-e $genome.'.ssi') && system("$SFETCH --index \"$genome\""))
		{
			die "\n  ERROR:  Failed to build .ssi index on genome '$genome'\n\n";
		}

		$SPECIES_TO_GENOME{$species} = $genome;
		$GENOME_LIST{$genome} = $species; # Just a little fun!
	
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

	$HMMBUILD = $th_src_dir.'hmmbuild';
	die "\n  ERROR:  Failed to locate hmmbuild in HMMER3 source directory (looking for '$HMMBUILD')\n\n"
		if (!(-x $HMMBUILD));

	$HMMSEARCHT = $th_src_dir.'hmmsearcht';
	die "\n  ERROR:  Failed to locate hmmsearcht in HMMER3 source directory (looking for '$HMMSEARCHT')\n\n"
		if (!(-x $HMMSEARCHT));


}


