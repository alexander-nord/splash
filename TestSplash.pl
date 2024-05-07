#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub Max { my $A = shift; my $B = shift; return $A if ($A>$B); return $B; }
sub Min { my $A = shift; my $B = shift; return $A if ($A<$B); return $B; }


sub FamilySplash;
sub BigBadSplash;
sub DetermineTargetSeq;
sub GenomeRangeFileToTargetFile;
sub HelpAndDie;
sub ReadChromosomeLengths;
sub ParseCommandArguments;
sub ConfirmInputsToSplashDir;
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
my %GENOME_LIST;        # So that we don't delete anything important
ParseCommandArguments();


# How long is a given chromosome on a species' genome?
my %SPECIES_CHR_TO_LEN; 
ReadChromosomeLengths();


# We survived commandline parsing, so let's build an
# output directory!
die "\n  ERROR:  Failed to create output directory '$OPTIONS{'output-dir'}'\n\n"
	if (system("mkdir \"$OPTIONS{'output-dir'}\""));

# Just in case...
my $ERROR_FILE = $OPTIONS{'output-dir'}.'Splash.err';


# Alright, let's figure this stuff out!
if    ($OPTIONS{'meta-dir'}  ) { BigBadSplash($OPTIONS{'protein-input'}); }
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

			# In case runtime was killed while there was a 'target' file,
			# we don't want to treat that as a query
			if ($file_base_name !~ /\.target$/ && !(-e $family_dir_name.$file_base_name.'.hmm')) 
			{
				push(@FilesToHMMBUILD,$family_dir_name.$file_name) 
			}
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

		push(@InputHMMs,$hmm_file_name);

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
			# system("cat $err_file_name >> $ERROR_FILE");

			if ($OPTIONS{'err-kills'})
			{
				unless ($GENOME_LIST{$target_file_name}) {
					system("rm \"$target_file_name\"");
					system("rm \"$target_file_name.ssi\"");
				}
				die "\n  ERROR: Failed during execution of command: $hmmsearcht_cmd\n\n";
			}

		} 
		else 
		{
			$num_fam_victories++;
		}


		# If we created a target sequence file, kill it!
		unless ($GENOME_LIST{$target_file_name}) {
			system("rm \"$target_file_name\"");
			system("rm \"$target_file_name.ssi\"");
		}

	}


	# If this is being run as part of a larger test,
	# we'll want to collect the number of error-full and
	# error-free runs.
	return ($num_fam_victories,$num_fam_errors);


}










###########################################################################
#
#  Subroutine: BigBadSplash
#
sub BigBadSplash
{
	
	my $protein_meta_dir_name = shift;

	# We're going to want to be able to parallelize this,
	# so we'll accumulate a list of directories and then
	# proceed from there
	my @FamilyDirs;

	opendir(my $ProteinMetaDir,$protein_meta_dir_name)
		|| die "\n  ERROR:  Failed to open protein meta-directory '$protein_meta_dir_name'\n\n";
	while (my $family = readdir($ProteinMetaDir)) 
	{

		$family =~ s/\/$//;
		my $family_dir_name = $protein_meta_dir_name.$family.'/';
		next if (!(-d $family_dir_name));

		push(@FamilyDirs,$family_dir_name);

	}
	closedir($ProteinMetaDir);


	my $num_cpus = $OPTIONS{'num-cpus'};
	my $num_fams = scalar(@FamilyDirs);
	if ($num_cpus > $num_fams) {
		$num_cpus = $num_fams;
	}


	my $thread_id = 0;
	my $active_threads = 1;
	my $pid = 0;
	while ($active_threads < $num_cpus) 
	{
		if ($pid = fork)
		{
			die "\n  ERROR: Fork failed\n\n" if (not defined $pid);
			$active_threads++;
		}
		else
		{
			$thread_id = $active_threads;
			last;
		}
	}


	my $start_fam_id =  $thread_id    * int($num_fams/$num_cpus);
	my   $end_fam_id = ($thread_id+1) * int($num_fams/$num_cpus);
	if ($end_fam_id > $num_fams) {
		$end_fam_id = $num_fams;
	}


	my $total_victories = 0; # "Non-errors" probably works better...
	my $total_errors    = 0; # How many times did hmmsearcht exit with non-0?
	for (my $fam_id = $start_fam_id; $fam_id < $end_fam_id; $fam_id++) 
	{

		my $family_dir_name = $FamilyDirs[$fam_id];
		$family_dir_name =~ /\/([^\/]+)\/$/;

		my $family = $1;

		my ($num_fam_victories,$num_fam_errors) = FamilySplash($FamilyDirs[$fam_id]);

		$total_victories += $num_fam_victories;
		$total_errors    += $num_fam_errors;

	}


	# End of the real work!
	exit(0) if ($thread_id);
	while (wait() != -1) {}


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
	$hmm_file_name =~ /^(.*\/)([^\/]+)\.hmm$/;
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
	elsif (-e $hmm_file_dir_name.$presumptive_species.'.genome-range.out')
	{
		$range_file_name = $hmm_file_dir_name.$presumptive_species.'.genome-range.out';
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


	# What's the range once we pull in some extra nucleotides?
	my $range = 0;
	if ($RangeFileData{'range'}) {

		$RangeFileData{'range'} =~ /^(\d+)\.\.(\d+)$/;
		my $range_start = $1;
		my $range_end   = $2;

		if ($OPTIONS{'num-extra-nucls'}) {
			$range_start = Max($range_start-$OPTIONS{'num-extra-nucls'},1                                     );
			$range_end   = Min(  $range_end+$OPTIONS{'num-extra-nucls'},$SPECIES_CHR_TO_LEN{$species.'|'.$chr});
		}

		$range = $range_start.'..'.$range_end;

	}


	# Put together the sfetch command
	my $sfetch_cmd = "$SFETCH -o $target_file_name";
	$sfetch_cmd = $sfetch_cmd." -c $range" if ($range);
	$sfetch_cmd = $sfetch_cmd." $genome $chr 1>/dev/null";


	# Clear the way!
	system("rm \"$target_file_name\"") if (-e $target_file_name);

	if (system($sfetch_cmd)) {
		die "\n  ERROR:  Subsequence range extraction failed (command:'$sfetch_cmd')\n\n";
	}


	# Build a lil' ssi on the target subsequence
	my $sfetch_index_cmd = "$SFETCH --index $target_file_name 1>/dev/null";
	if (system($sfetch_index_cmd)) {
		die "\n  ERROR:  Failed to build .ssi index for target sequence '$target_file_name'\n\n";
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
	print "            :  ./TestSplash.pl {OPT.S} [gene]\n";
	print "            '-------------------------------------------------------\n";
	print "\n";
	print "\n";
	print "  Use Case 2:  Search ALL GENES in an 'inputs-to-splash' directory\n";
	print "            :\n";
	print "            :  ./TestSplash.pl {OPT.S} [path/to/inputs-to-splash]\n";
	print "            '-------------------------------------------------------\n";
	print "\n";
	print "\n";
	print "  OPT.S: --full-genome : Force use of full genome as target sequence.\n";
	print "         --err-kills   : If an hmmsearcht run fails, kill the script\n";
	print "                         (by default we log the error and continue)\n";
	print "\n";
	die   "\n";
}









###########################################################################
#
#  Subroutine: ReadChromosomeLengths
#
sub ReadChromosomeLengths
{
	foreach my $species (sort keys %SPECIES_TO_GENOME) {

		my $genome = $SPECIES_TO_GENOME{$species};

		my $seqstat_file_name = $genome.'.seqstat.out';
		if (!(-e $seqstat_file_name) && system("$SEQSTAT -a \"$genome\" > \"$seqstat_file_name\"")) {
			die "\n  ERROR:  Failed to locate/create expected seqstat file ($seqstat_file_name)\n\n";
		}

		open(my $SeqstatFile,'<',$seqstat_file_name)
			|| die "\n  ERROR:  Failed to open seqstat file '$seqstat_file_name'\n\n";
		while (my $line = <$SeqstatFile>) {
			$line =~ s/\n|\r//g;
			if ($line =~ /^=\s+(\S+)\s+(\d+)\s*$/) {
				my $chr = $1;
				my $len = $2;
				$SPECIES_CHR_TO_LEN{$species.'|'.$chr} = $len;
			}
		}
		close($SeqstatFile);

	}
}









###########################################################################
#
#  Subroutine: ParseCommandArguments
#
sub ParseCommandArguments
{

	# Default output directory name
	my $default_out_dir_name = 'Splash-Test-Results';
	my $out_dir_name = $default_out_dir_name;
	my $attempt = 1;
	while (-d $out_dir_name) {
		$attempt++;
		$out_dir_name = $default_out_dir_name.'-'.$attempt;
	}
	$OPTIONS{'output-dir'} = $out_dir_name.'/';


	# How many CPUs do we want?
	$OPTIONS{'num-cpus'} = 1;


	# Unless specified, we're looking for data in this directory
	$OPTIONS{'inputs-dir'} = '../inputs-to-splash/';
	$OPTIONS{'genome-dir'} = '../inputs-to-splash/genomes/';


	# How much sequence do we want to pull in around any tightly
	# defined nucleotide ranges?
	$OPTIONS{'num-extra-nucls'} = 5000;
	$OPTIONS{'err-kills'}       =    0;


	my $num_args = scalar(@ARGV);
	for (my $arg_id=0; $arg_id<$num_args; $arg_id++) {

		my $Arg = $ARGV[$arg_id];

		if (lc($Arg) =~ '^-?-?full-genome$') 
		{
			$OPTIONS{'full-genome'} = 1; # Overrides existence of '.genome-range.out' files
		}
		elsif (lc($Arg =~ /^-?-?err-kills$/))
		{
			$OPTIONS{'err-kills'} = 1; # First failed run terminates script
		}
		elsif (lc($Arg) =~ /^-?-?cpus$/ || lc($Arg) =~ /^-?-n$/)
		{
			$arg_id++;
			my $num_cpus = $ARGV[$arg_id];
			if ($num_cpus !~ /\D/) 
			{
				$OPTIONS{'num-cpus'} = int($num_cpus);
			} 
			else
			{
				die "\n  ERROR:  Apparent request for multiple CPUs was not followed by an integer? ($ARGV[$arg_id-1],$ARGV[$arg_id])\n\n";
			}
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
				# This is (assumed to be) a path to an 'inputs-to-splash'
				# directory.  Verify!
				# If so, 'protein-input' now names the 'protein-data' path
				$OPTIONS{'protein-input'} = ConfirmInputsToSplashDir($Arg);
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
#  Subroutine: ConfirmInputsToSplashDir
#
sub ConfirmInputsToSplashDir
{
	my $dir_name = shift;
	$dir_name = $dir_name.'/' if ($dir_name !~ /\/$/);

	if (!(-d $dir_name)) {
		die "\n  ERROR:  Presumed 'inputs-to-splash' directory was not as expected? ($dir_name)\n\n";
	}

	my $genome_dir_name = $dir_name.'genomes/';
	if (!(-d $genome_dir_name)) {
		die "\n  ERROR:  Expected genome directory '$genome_dir_name' (ConfirmInputsToSplashDir)\n\n";
	}

	my $protein_meta_dir_name = $dir_name.'protein-data/';
	if (!(-d $protein_meta_dir_name)) {
		die "\n  ERROR:  Expected protein family meta-directory '$protein_meta_dir_name' (ConfirmInputsToSplashDir)\n\n";
	}

	# Looks good enough to consider moving forward!
	$OPTIONS{'genome-dir'} = $genome_dir_name;
	return $protein_meta_dir_name;

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

	$SFETCH = $miniapps_dir.'esl-sfetch';
	die "\n  ERROR:  Failed to locate required Easel tool sfetch (looking for '$SFETCH')\n\n"
		if (!(-x $SFETCH));

	$SEQSTAT = $miniapps_dir.'esl-seqstat';
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


