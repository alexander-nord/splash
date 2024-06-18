#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub Max { my $A = shift; my $B = shift; return $A if ($A>$B); return $B; }
sub Min { my $A = shift; my $B = shift; return $A if ($A<$B); return $B; }


sub AggregateAllResults;
sub RunPostSlurmAggregation;
sub CompileBasicResults;
sub GetSumPctCoverage;
sub FamilySplash;
sub BigBadSplash;
sub DetermineTargetSeq;
sub GenomeRangeFileToTargetFile;
sub HelpAndDie;
sub ReadChromosomeLengths;
sub ParseCommandArguments;
sub GetTildeDirName;
sub ValidateSpeciesGuide;
sub ConfirmGenome;
sub ConfirmRequiredTools;
sub AnticipateMiniprot;
sub ConfirmMiniprot;
sub GetPct;




# We'll make sure we have the necessary tools before proceeding
my $SFETCH;
my $SEQSTAT;
# my $HMMBUILD;   # We been knowing how to BATH!
# my $HMMSEARCHT; # BATH-splashing puppy!
my $BATHBUILD;
my $BATHCONVERT;
my $BATHSEARCH;
ConfirmRequiredTools();


# If the user didn't provide any arguments, help 'em out
HelpAndDie() if (@ARGV < 2);


# We'll pre-load where we *expect* miniprot to be, in case we want it
my $MINIPROT = AnticipateMiniprot();


# What do they want us to even do?
my %OPTIONS;
my %SPECIES_TO_GENOME;
my %GENOME_LIST;        # So that we don't delete anything important
ParseCommandArguments();

# How many threads per job?  
# Currently, we don't support changing this through the commandline.
$OPTIONS{'threads-per-job'} = 8;


# How long is a given chromosome on a species' genome?
my %SPECIES_CHR_TO_LEN; 
ReadChromosomeLengths();


# We survived commandline parsing, so let's build an output directory!
# NOTE that if we're using slurm, we assume that the output directory
#  has already been built for us.
if (!$OPTIONS{'slurm'} && system("mkdir \"$OPTIONS{'output-dir'}\""))
{
	die "\n  ERROR:  Failed to create output directory '$OPTIONS{'output-dir'}'\n\n";
}


# Just in case...
my $ERROR_FILE = $OPTIONS{'output-dir'}.'Splash.err';


# Alright, let's figure this stuff out!
if    ($OPTIONS{'single-gene'}  ) { FamilySplash($OPTIONS{'protein-input'}); }
elsif ($OPTIONS{'protein-input'}) { BigBadSplash($OPTIONS{'protein-input'}); }
else  {     die "\n  ERROR:  Failed to recognize protein input\n\n";         }


1;














###########################################################################
#
#  Subroutine: AggregateAllResults
#
sub AggregateAllResults
{
	my $all_fams_dir_name = shift;
	my $panther_species   = shift;


	# The big 'Summary' file will be written at the end of this
	# subroutine, but we'll build up this CSV file as we go
	# (containing info related to how successful we are at splicing
	# hits to provide full-model coverage)
	my $coverage_csv_file_name;
	if ($panther_species)
	{
		$coverage_csv_file_name = $OPTIONS{'output-dir'}.$panther_species.'-Coverage.csv';
	}
	else
	{
		$coverage_csv_file_name = $OPTIONS{'output-dir'}.'Coverage.csv';
	}

	open(my $CoverageFile,'>',$coverage_csv_file_name)
		|| die "\n  ERROR:  Failed to open coverage-describing output file '$coverage_csv_file_name'\n\n";


	# The header for our coverage file
	print $CoverageFile "Gene,Sequence-ID,Model-Length,Num-Exon-Sets,";
	print $CoverageFile "Model-Coverage-of-Best-Exon-Set,Best-ES-Pct-Coverage,";
	print $CoverageFile "Best-ES-Start-Pos,Best-ES-End-Pos,";
	print $CoverageFile "Sum-ES-Pct-Coverage,BATH-Runtime,Time-in-Splash\n";


	# In case we had weird / discordant output for any inputs,
	# be sure to report that!
	my $discord_file_name = $OPTIONS{'output-dir'}.'Discord.err';


	# We'll want to capture the splicing dinucleotides
	my %ThreePrime;
	my %ThreePrimeOne;
	my %ThreePrimeTwo;
	my %FivePrime;
	my %FivePrimeOne;
	my %FivePrimeTwo;


	# For large-scale debugging, do we see splice signals that don't
	# match standard nucleotide codes?
	my %ValidSpliceNucls;
	$ValidSpliceNucls{'a'} = 1;
	$ValidSpliceNucls{'c'} = 1;
	$ValidSpliceNucls{'g'} = 1;
	$ValidSpliceNucls{'t'} = 1;
	$ValidSpliceNucls{'n'} = 1;
	$ValidSpliceNucls{'-'} = 1;


	# How many sequences had some splicing action?
	# How many covered the full length of the model?
	# How many (of each of the above) incorporated 'missed' exons?
	my $total_input_phmms  = 0;
	my $total_num_spliced  = 0;
	my $num_full_model     = 0;
	my $inputs_with_missed = 0;
	my $full_with_missed   = 0;


	my @AllFamsList;
	opendir(my $AllFamsDir,$all_fams_dir_name)
		|| die "\n  ERROR:  Failed to open all-family meta-directory '$all_fams_dir_name'\n\n";
	while (my $family = readdir($AllFamsDir)) 
	{		
		
		next if ($family =~ /^\./);
		$family =~ s/\/$//;

		if ($panther_species)
		{
			next if (!(-e $all_fams_dir_name.$family.'/'.$panther_species.'-summary.out'));
		}
		else
		{
			next if (!(-e $all_fams_dir_name.$family.'/summary.out'));
		}
	
		push(@AllFamsList,$family);
	
	}
	closedir($AllFamsDir);



	foreach my $family (sort @AllFamsList) 
	{

		my $family_dir_name = $all_fams_dir_name.$family.'/';
		my $summary_file_name;
		if ($panther_species)
		{
			$summary_file_name = $family_dir_name.$panther_species.'-summary.out';
		}
		else 
		{
			$summary_file_name = $family_dir_name.'summary.out';
		}


		# Prepare to learn!
		open(my $SummaryFile,'<',$summary_file_name)
			|| die "\n  ERROR:  Failed to open summary file '$summary_file_name' for family '$family'\n\n";


		# First off, grab the metadata that will let us navigate
		# this summary file smoothtly
		my $num_fam_inputs;
		while (my $line = <$SummaryFile>)
		{
			if ($line =~ /NUM SPLASHED : (\d+)/)
			{
				$num_fam_inputs = $1;
				last;
			}
		}


		$total_input_phmms += $num_fam_inputs;


		# Let's interrogate some stinky ol' exon sets!
		for (my $fam_input_id=0; $fam_input_id<$num_fam_inputs; $fam_input_id++)
		{

			my $line = <$SummaryFile>;
			$line = <$SummaryFile> while ($line !~ /Input ID/);


			$line =~ /Input ID\s+: (\S+)/;
			my $fam_input_name = $1;


			# There was *something* reported -- that's a splicing!
			$total_num_spliced++;


			# How long did we spend running BATH, altogether?
			$line = <$SummaryFile>;
			$line =~ /BATH Runtime\s+: (\S+)/;
			my $fam_input_bath_runtime = $1;


			# How long did we spend inside the *splash* zone?
			$line = <$SummaryFile>;
			$line =~ /Splash Runtime\s+: (\S+)/;
			my $fam_input_splash_runtime = $1;


			# Model length
			$line = <$SummaryFile>;
			$line =~ /Model Length\s+: (\d+)/;
			my $model_length = $1;


			# How many exon sets did we arrive at?
			$line = <$SummaryFile>;
			$line =~ /Num Exon Sets\s+: (\d+)/;
			my $num_exon_sets = $1;


			# Did this hit get us from 0 to gm->M?
			my $full_coverage = 0;
			if ($line =~ /Full model coverage/) {
				$full_coverage = 1;
				$num_full_model++;
			}


			# Even if we know we had full coverage, we'll compile
			# these bits of information independently
			my $best_coverage_len = 0; # Model positions covered by an exon set
			my $best_coverage_start;
			my $best_coverage_end;
			my $best_coverage_set_id; # Probably useless, but *whatever*


			# We'll also want to know the sum coverage achieved by
			# all of the spliced exons
			my @ExonSetModelRanges;


			# In case this hit looked good according to the splice
			# graph but we couldn't make our data play nicely with
			# the model, we'll need to be prepared to adjust our
			# statistics to represent it as a failure (sadly)
			my $num_discord_exon_sets = 0;


			# Let's start running through our dang exon sets
			my $had_missed_exons = 0;
			for (my $exon_set_id = 1; $exon_set_id <= $num_exon_sets; $exon_set_id++) 
			{

				$line = <$SummaryFile>; # Exon set ID
				
				
				$line = <$SummaryFile>; # Num. Exons
				$line =~ /Num Exons\s+: (\d+)/;
				my $num_exons = $1;


				# Model positions covered by this exon set
				$line = <$SummaryFile>;
				$line =~ /Model Range\s+: (\d+)\.\.(\d+)/;
				my $exon_set_model_start = $1;
				my $exon_set_model_end   = $2;
				my $exon_set_coverage    = 1 + ($exon_set_model_end - $exon_set_model_start);


				# NOTE: Because HMMER sometimes likes to chop off the
				#       first/last residues, we'll see if this is
				#       close enough to count as an honorary full model
				#       alignment
				if ($full_coverage == 0 && $exon_set_model_start < 6 && $model_length - $exon_set_model_end < 6) {
					$full_coverage = 1;
					$num_full_model++;
				}



				# Chromosome name and nucleotide range summary for 
				# this exon set
				$line = <$SummaryFile>;
				$line = <$SummaryFile>;


				$line = <$SummaryFile>; # Missed exons?
				my $exon_set_missed_exons = 1 if ($line =~ /\? Yes/);


				# The next line indicates if the secondary search with
				# the full pipeline doesn't match what we fed in.
				# If this is the case, there aren't any exon-by-exon
				# outputs.
				$line = <$SummaryFile>;
				if ($line =~ /EXON SET DISCORD/) {
					$num_discord_exon_sets++;
					next;
				}


				# Since this exon set isn't "discordant" we can include
				# it in our understanding of the total coverage of the
				# model.
				push(@ExonSetModelRanges,$exon_set_model_start.'..'.$exon_set_model_end);


				# Now that we know we didn't have discord with the
				# model behavior, we can report this hit's stats
				# (if applicable)
				if ($exon_set_coverage > $best_coverage_len)
				{
					$best_coverage_set_id = $exon_set_id;
					$best_coverage_start  = $exon_set_model_start;
					$best_coverage_end    = $exon_set_model_end;
					$best_coverage_len    = $exon_set_coverage;
				}
				$had_missed_exons = 1 if ($exon_set_missed_exons);


				# Grab info about splice site signal
				for (my $exon_id = 1; $exon_id <= $num_exons; $exon_id++)
				{

					$line = <$SummaryFile> unless ($exon_id == 1); # Exon ID
					$line = <$SummaryFile>; # Model range
					$line = <$SummaryFile>; # Nucl. range


					$line = <$SummaryFile>; # Splice signal
					if ($line !~ /\S\S\|\S\S/) {
						$line =~ s/\n|\r//g;
					}

					$line =~ /: (\S)(\S)\|(\S)(\S)/;

					my $three_one = lc($1);
					my $three_two = lc($2);
					my $five_one  = lc($3);
					my $five_two  = lc($4);


					# DEBUGGING
					if (!($ValidSpliceNucls{$three_one} && $ValidSpliceNucls{$three_two} && $ValidSpliceNucls{$five_one} && $ValidSpliceNucls{$five_two})) 
					{
						print "Unexpected splice signature observed in '$fam_input_name'\n";
						next;
					}


					my $three = $three_one.$three_two;
					my $five  = $five_one.$five_two;

					if ($three_one ne '-') 
					{
						if ($ThreePrimeOne{$three_one}) { $ThreePrimeOne{$three_one}++; }
						else                            { $ThreePrimeOne{$three_one}=1; }
					}
					if ($three_two ne '-') 
					{
						if ($ThreePrimeTwo{$three_two}) { $ThreePrimeTwo{$three_two}++; }
						else                            { $ThreePrimeTwo{$three_two}=1; }
					}
					if ($three ne '--')
					{
						if ($ThreePrime{$three}) { $ThreePrime{$three}++; }
						else                     { $ThreePrime{$three}=1; }
					}


					if ($five_one ne '-') 
					{
						if ($FivePrimeOne{$five_one}) { $FivePrimeOne{$five_one}++; }
						else                          { $FivePrimeOne{$five_one}=1; }
					}
					if ($five_two ne '-') 
					{
						if ($FivePrimeTwo{$five_two}) { $FivePrimeTwo{$five_two}++; }
						else                          { $FivePrimeTwo{$five_two}=1; }
					}
					if ($five ne '--')
					{
						if ($FivePrime{$five}) { $FivePrime{$five}++; }
						else                   { $FivePrime{$five}=1; }
					}

				}

			}


			# We need to have had at least one non-discordant hit
			# in order to log this as a success.  Otherwise we'll
			# backtrack.
			if ($num_discord_exon_sets != $num_exon_sets)
			{

				# Compute the percent of the model covered by our best
				# exon set
				my $best_coverage_pct = GetPct($best_coverage_len,$model_length);


				# Determine the total coverage of the model by all spliced
				# exons
				my $sum_coverage_pct = GetSumPctCoverage(\@ExonSetModelRanges,$model_length);


				# Hell yeah! We sure were given an input, did something with
				# it, and then produced an output!  Way to go, us!
				print $CoverageFile "$family,$fam_input_name,$model_length,";
				print $CoverageFile "$num_exon_sets,$best_coverage_len,$best_coverage_pct\%,";
				print $CoverageFile "$best_coverage_start,$best_coverage_end,$sum_coverage_pct\%,";
				print $CoverageFile "$fam_input_bath_runtime,$fam_input_splash_runtime\n";


				# Wrap up this input with whether it made use of missed exons
				if ($had_missed_exons)
				{
					$inputs_with_missed++;
					$full_with_missed++ if ($full_coverage);
				}

			}
			else
			{
				# BUMMER!
				$num_full_model-- if ($full_coverage);
				$total_num_spliced--;
				system("echo \"$family/$fam_input_name\" >> $discord_file_name");
			}

		}


		close($SummaryFile);


	}


	# Coverage file is all sewn up!
	close($CoverageFile);



	# For percentages, we'll want these totals
	my $total_3prime_sites = 0;
	foreach my $dinucl (keys %ThreePrime) 
	{
		$total_3prime_sites += $ThreePrime{$dinucl};
	}

	my $total_5prime_sites = 0;
	foreach my $dinucl (keys %FivePrime)
	{
		$total_5prime_sites += $FivePrime{$dinucl};
	}



	# Now that we've gathered all of these summary statistics,
	# it's time to make the good news heard!
	my $final_results_file_name;
	if ($panther_species)
	{
		$final_results_file_name = $OPTIONS{'output-dir'}.$panther_species.'-Splash-Summary.out';
	}
	else
	{
		$final_results_file_name = $OPTIONS{'output-dir'}.'Splash-Summary.out';
	}

	open(my $FinalResults,'>',$final_results_file_name)
		|| die "\n  ERROR:  Failed to open final output file '$final_results_file_name'\n\n";


	my $pct_spliced       = GetPct($total_num_spliced,$total_input_phmms);
	my $pct_full_model    = GetPct($num_full_model,$total_input_phmms);
	my $pct_with_missed   = GetPct($inputs_with_missed,$total_input_phmms);
	my $pct_full_w_missed = GetPct($full_with_missed,$total_input_phmms);


	# First-blush statistics
	print $FinalResults "Total Number of Input pHMMs       : $total_input_phmms\n";
	print $FinalResults "Inputs w/ Successful Splicing     : $total_num_spliced ($pct_spliced\%)\n";
	print $FinalResults "Inputs w/ Full-Model Splicing     : $num_full_model ($pct_full_model\%)\n";
	print $FinalResults "Inputs w/ Hits Using Missed Exons : $inputs_with_missed ($pct_with_missed\%)\n";
	print $FinalResults "Full-Model Hits w/ Missed Exons   : $full_with_missed ($pct_full_w_missed\%)\n";


	# The BAD ending
	if ($total_num_spliced == 0)
	{
		close($FinalResults);
		return;
	}


	# What did splice signals look like?
	print $FinalResults "\n";
	print $FinalResults "3' Splice Signals\n";
	print $FinalResults "- Dinucleotides\n";
	foreach my $dinucl (sort keys %ThreePrime)
	{
		my $count = $ThreePrime{$dinucl};
		my $pct   = GetPct($count,$total_3prime_sites);
		print $FinalResults "  $dinucl : $count ($pct\%)\n";
	}
	print $FinalResults "- [-2] position\n";
	foreach my $nucl (sort keys %ThreePrimeOne)
	{
		my $count = $ThreePrimeOne{$nucl};
		my $pct   = GetPct($count,$total_3prime_sites);
		print $FinalResults "  $nucl  : $count ($pct\%)\n";
	}
	print $FinalResults "- [-1] position\n";
	foreach my $nucl (sort keys %ThreePrimeTwo)
	{
		my $count = $ThreePrimeTwo{$nucl};
		my $pct   = GetPct($count,$total_3prime_sites);
		print $FinalResults "   $nucl : $count ($pct\%)\n";
	}

	print $FinalResults "\n";
	print $FinalResults "5' Splice Signals\n";
	print $FinalResults "- Dinucleotides\n";
	foreach my $dinucl (sort keys %FivePrime)
	{
		my $count = $FivePrime{$dinucl};
		my $pct   = GetPct($count,$total_5prime_sites);
		print $FinalResults "  $dinucl : $count ($pct\%)\n";
	}
	print $FinalResults "- [+1] position\n";
	foreach my $nucl (sort keys %FivePrimeOne)
	{
		my $count = $FivePrimeOne{$nucl};
		my $pct   = GetPct($count,$total_5prime_sites);
		print $FinalResults "  $nucl  : $count ($pct\%)\n";
	}
	print $FinalResults "- [+2] position\n";
	foreach my $nucl (sort keys %FivePrimeTwo)
	{
		my $count = $FivePrimeTwo{$nucl};
		my $pct   = GetPct($count,$total_5prime_sites);
		print $FinalResults "   $nucl : $count ($pct\%)\n";
	}
	print $FinalResults "\n";


	close($FinalResults);


}










###########################################################################
#
#  Subroutine: RunPostSlurmAggregation
#
sub RunPostSlurmAggregation
{

	my $splash_results_dir_name = $OPTIONS{'post-slurm'};

	my $rbg_dir_name = $splash_results_dir_name.'Results-by-Gene/';
	die "\n  ERROR:  Failed to locate '$rbg_dir_name'\n\n"
		if (!(-d $rbg_dir_name));

	$OPTIONS{'output-dir'} = $splash_results_dir_name;

	foreach my $species (keys %SPECIES_TO_GENOME)
	{
		AggregateAllResults($rbg_dir_name,$species);
	}

}










###########################################################################
#
#  Subroutine: CompileBasicResults
#
sub CompileBasicResults
{

	my $out_dir_name    = shift;
	my $base_names_ref  = shift;
	my $panther_species = shift;

	
	$out_dir_name =~ /\/([^\/]+)\/[^\/]+$/;
	my $family = $1;

	
	my $summary_file_name;
	if ($panther_species)
	{
		$summary_file_name = $out_dir_name.$panther_species.'-summary.out';
	}
	else
	{ 
		$summary_file_name = $out_dir_name.'summary.out';
	}

	open(my $SummaryFile,'>',$summary_file_name)
		|| die "\n  ERROR:  Failed to open summary file '$summary_file_name'\n\n";

	
	my @BaseNames = @{$base_names_ref};


	# We'll need to reduce to just the hits for the species
	# of interest if we're working with PANTHER
	if ($panther_species)
	{
		my @SpeciesBaseNames;
		foreach my $base_name (@BaseNames)
		{

			$base_name =~ /\.([^\.]+)$/;
			my $species = $1;

			push(@SpeciesBaseNames,$base_name)
				if ($species eq $panther_species);		

		}
		@BaseNames = @SpeciesBaseNames;
	}


	my $num_seqs = scalar(@BaseNames);


	print $SummaryFile "GENE FAMILY  : $family\n";
	print $SummaryFile "NUM SPLASHED : $num_seqs\n";



	# Now we'll iterate over all of the (applicable) outputs
	# and compile some super rad info!
	foreach my $base_name (sort @BaseNames)
	{

		my $spl_out_file_name = $out_dir_name.$base_name.'.out';
		my $spl_err_file_name = $out_dir_name.$base_name.'.err';

		# Just to get it out of the way, let's
		# grab the timing data from the error file
		my $splash_runtime = 'Not timed';
		my $bath_runtime   = 'Not timed';
		if (-e $spl_err_file_name)
		{
	
			open(my $SplErrFile,'<',$spl_err_file_name)
				|| die "\n  ERROR:  Failed to open error file '$spl_err_file_name' (looking for timing info.)\n\n";
			
			while (my $line = <$SplErrFile>) 
			{
				if ($line =~ /Time spent splashing/)
				{
					$line = <$SplErrFile>;
					if ($line =~ /Elapsed:\s+(\S+)/) 
					{
						$splash_runtime = $1;
					}
				}
				elsif ($line =~ /Elapsed \(wall clock\) time.*:\s+(\S+)\s*$/)
				{
					$bath_runtime = $1;
					if ($bath_runtime =~ /^(\d+):([^:]+)$/)
					{

						my $minutes = $1;
						my $seconds = $2;
						$minutes = '0'.$minutes if ($minutes < 10);

						$bath_runtime = '00:'.$minutes.':'.$seconds;

					}
				}
			}
			
			close($SplErrFile);
		
		}



		# Get this stuff taken care of
		print $SummaryFile "\n";
		print $SummaryFile "  Input ID       : $base_name\n";
		print $SummaryFile "  BATH Runtime   : $bath_runtime\n";
		print $SummaryFile "  Splash Runtime : $splash_runtime\n";


		# Now, it's time to look at the actual output!
		open(my $SplOutFile,'<',$spl_out_file_name)
			|| die "\n  ERROR:  Failed to open output file '$spl_out_file_name' (BAD NEWS!)\n\n";

		while (my $line = <$SplOutFile>) {
			next if ($line !~ /^\s*Query:\s+\S+\s+\[M=(\d+)\]\s*$/);
			print $SummaryFile "  Model Length   : $1\n";
			last;
		}


		# This is a bit inefficient, but it's useful to put this at the
		# top for each query
		my $num_exon_sets = 0;
		my $has_full_hit  = 0;
		my $chromosome;
		my $chr_start;
		while (my $line = <$SplOutFile>) 
		{
			
			$num_exon_sets++  if ($line =~ /\| = Exon Set/);
			$has_full_hit = 1 if ($line =~ /\(\* Full Model\)/);

			if ($line =~ /\| = Target Seq Name\s+(\S+)/)
			{
				$chromosome = $1;
			}

			# This should capture a nucleotide sequence row,
			# if we were searching within a target sequence
			# that was esl-sfetch'd from a larger sequence
			if (!$chr_start && $line =~ /^\s+\S+\/(\d+)-\d+\s+\d+/)
			{
				$chr_start  = $2;
			}

		}
		print $SummaryFile "  Num Exon Sets  : $num_exon_sets";


		# If it doesn't look like we were searching within a
		# sub-range of the chromosome, we were starting at 1, baby!
		$chr_start = 1 if (!$chr_start);


		# We can show our excitement *just a bit* in this case:
		if ($has_full_hit && $num_exon_sets == 1)
		{
			print $SummaryFile " (Full model coverage!)";
		}
		print $SummaryFile "\n";


		# Close and reopen the file :p
		close($SplOutFile);
		open($SplOutFile,'<',$spl_out_file_name)
			|| die "\n  ERROR:  Failed to open output file '$spl_out_file_name' (BAD NEWS!)\n\n";
		for (my $exon_set_id = 1; $exon_set_id <= $num_exon_sets; $exon_set_id++) 
		{

			print $SummaryFile "  - Exon Set     : $exon_set_id\n";


			# Advance to the start of the next exon set
			my $num_exons;
			while (my $line = <$SplOutFile>)
			{
				if ($line =~ /= Exon Set \d+ \((\d+) exons\)/) 
				{
					$num_exons = $1;
					last;
				}
			}
			print $SummaryFile "  - Num Exons    : $num_exons\n";


			# What range of model / nucleotide positions are covered?
			# NOTE: We'll need to correct these!
			my $model_pos_line = <$SplOutFile>;
			$model_pos_line =~ /(\d+)\.\.(\d+)/;
			my $model_start = $1;
			my $model_end   = $2;

			print $SummaryFile "  - Model Range  : $model_start..$model_end\n";


			my $target_seq_line = <$SplOutFile>;


			my $nucl_coord_line = <$SplOutFile>;
			$nucl_coord_line =~ /(\d+)\.\.(\d+)/;
			my $nucl_start = $1 + ($chr_start-1);
			my $nucl_end   = $2 + ($chr_start-1);

			# We *always* need this.
			my $revcomp = 0;
			$revcomp = 1 if ($nucl_start > $nucl_end);


			print $SummaryFile "  - Chromosome   : $chromosome\n";
			print $SummaryFile "  - Nucl. Range  : $nucl_start..$nucl_end\n";


			# Did this hit incorporate "missed" exons?
			# Is it significantly divergent from the
			# exon set supplied for search?
			my $has_missed_exons = 0;
			my $exon_set_discord = 0;
			while (my $extra_info_line = <$SplOutFile>) 
			{
				last if ($extra_info_line !~ /^\|/);

				if ($extra_info_line =~ /\| \+ Includes Missed Exons/)
				{
					$has_missed_exons = 1;
				}
				elsif ($extra_info_line =~ /\| \+ WARNING/)
				{
					$exon_set_discord = 1;
				}
			}

			print $SummaryFile "  - Missed Exons ? ";
			if ($has_missed_exons) { print $SummaryFile "Yes\n"; }
			else                   { print $SummaryFile "No\n";  }


			# If there was discord, we'll signal that this was
			# a problematic hit and skip to the next exon set.
			if ($exon_set_discord)
			{
				print $SummaryFile "  - EXON SET DISCORD -- Spliced alignment abandoned\n";
				next;
			}


			# Now we can iterate over the exons that make up this
			# "exon set" and report their particular datums!
			my $exon_id = 0;
			my $exon_model_start = 0;
			my $exon_model_end;
			my $exon_nucl_start = 0;
			my $exon_nucl_end;
			my $exon_nucl_seq = ''; # Stupid lazy way to get the splice dinucleotides...
			while ($exon_id <= $num_exons && !eof($SplOutFile))
			{

				my $line = <$SplOutFile>;

				$line =~ s/\n|\r//g;
				next if (!$line);


				# Are we starting off a new exon? Are we at *the end*?!
				if ($line =~ /\[ Exon Set/ || $line =~ /^\+=========/)
				{

					# Don't sweat it, special baby zero
					if ($exon_id == 0)
					{
						$exon_id = 1;
						next;
					}


					# We need to account for the fact that the splice signal
					# eats two of our precious nucleotides
					#
					# BUT: Don't trim the start of the first exon or the end
					#      of the last exon!
					#
					if ($revcomp) {
						$exon_nucl_start -= 2 if ($exon_id > 1         );
						$exon_nucl_end   += 2 if ($exon_id < $num_exons);
					} else {
						$exon_nucl_start += 2 if ($exon_id > 1         );
						$exon_nucl_end   -= 2 if ($exon_id < $num_exons);
					}

					
					# What's the splice signal, Kenneth?
					my $exon_3prime = '--';
					if ($exon_id > 1)
					{
						$exon_nucl_seq =~ /^(\S\S)/;
						$exon_3prime = $1;
					}

					my $exon_5prime = '--';
					if ($exon_id < $num_exons)
					{
						$exon_nucl_seq =~ /(\S\S)$/;
						$exon_5prime = $1;
					}


					print $SummaryFile "    + Exon $exon_id\n";
					print $SummaryFile "      Model Range   : $exon_model_start..$exon_model_end\n";
					print $SummaryFile "      Nucl. Range   : $exon_nucl_start..$exon_nucl_end\n";
					print $SummaryFile "      Splice Signal : $exon_3prime|$exon_5prime\n";

					$exon_model_start = 0;
					$exon_nucl_start  = 0;
					$exon_nucl_seq    = '';

					$exon_id++;
					next;

				}


				# Are we hitting a *chill* line?
				if ($exon_id && $line =~ /^\s+\S+\s+(\d+).*\s+(\d+)\s*$/)
				{
					
					my $ali_line_model_start = $1;
					my $ali_line_model_end   = $2;

					$exon_model_start = $ali_line_model_start if (!$exon_model_start);
					$exon_model_end   = $ali_line_model_end;


					$line = <$SplOutFile>; # Quality line
					$line = <$SplOutFile>; # Translation line
					$line = <$SplOutFile>; # Nucleotides!


					$line =~ /^\s+\S+\s+(\d+)\s+(\S+)\s+(\d+)\s*$/;
					my $ali_line_nucl_start = $1;
					my $ali_line_nucl_seq   = $2;
					my $ali_line_nucl_end   = $3;


					# DEBUGGING
					if (!$ali_line_nucl_end) {
						$line =~ s/\n|\r//g;
						print "\n  Warning:  Nucleotide line for $base_name is missing ali_line_nucl_end (line:'$line')\n\n";
					}


					# NOTE: We'll deal with the add/sub-2 for splice signal *later*
					$exon_nucl_start = $ali_line_nucl_start + ($chr_start-1) if (!$exon_nucl_start);
					$exon_nucl_end   = $ali_line_nucl_end   + ($chr_start-1);
					$line = <$SplOutFile>; # PP line (I'm a child)

					$exon_nucl_seq = $exon_nucl_seq.$ali_line_nucl_seq;

				}

			}

		}

		close($SplOutFile);

	}

	close($SummaryFile);

}









###########################################################################
#
#  Subroutine: GetSumPctCoverage
#
sub GetSumPctCoverage
{

	my $model_ranges_ref = shift;
	my $model_length = shift;

	my @ModelRanges = @{$model_ranges_ref};


	# This is the stupid AND lazy way of doing this,
	# so it's perfect for me!
	my @Covered;
	for (my $i=1; $i<=$model_length; $i++) 
	{
		$Covered[$i] = 0;
	}


	foreach my $model_range (@ModelRanges)
	{
		
		$model_range =~ /^(\d+)\.\.(\d+)$/;
		my $model_start = $1;
		my $model_end   = $2;

		for (my $i=$model_start; $i<=$model_end; $i++)
		{
			$Covered[$i] = 1;
		}

	}


	my $sum_coverage = 0;
	for (my $i=1; $i<=$model_length; $i++)
	{
		$sum_coverage += $Covered[$i];
	}


	return GetPct($sum_coverage,$model_length);

}









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


	# ... but I can prepare for bad presumptions!
	if (!$gene)
	{
		my $no_gene_message = "Error: Failed to pull gene name from directory path $family_dir_name";
		system("echo \"$no_gene_message\" >> $ERROR_FILE");
		return (0,0,0);
	}


	# Read through the family's directory, figuring
	# out what the situation is with our query data
	opendir(my $FamilyDir,$family_dir_name) 
		|| die "\n  ERROR:  Failed to open directory '$family_dir_name' (FamilySplash)\n\n";
	my @InputHMMs;
	my @FilesToBATHCONVERT;
	my @FilesToBATHBUILD;
	my @InputConsSeqs; # For miniprot tests
	while (my $file_name = readdir($FamilyDir)) 
	{

		next if ($file_name =~ /^\./);

		if ($OPTIONS{'miniprot'})
		{

			if (lc($file_name) =~ /\.cons\.fa$/)
			{
				# Too easy!
				push(@InputConsSeqs,$family_dir_name.$file_name);
			}

		}
		else
		{

			# Have we stumbled upon a real-life HMM?!
			if (lc($file_name) =~ /\.hmm$/) 
			{
				# Yes, but is it... a BATH HMM?!
				if (lc($file_name) =~ /\.bath\.hmm$/) 
				{
				
					# Very nice!
					push(@InputHMMs,$family_dir_name.$file_name);
				
				}
				else
				{
					# THIS isn't a BATH HMM... but it could be!
					# If there isn't already a BATH HMM that looks like
					# this one, queue it up for conversion.
					my $bath_file_name = $family_dir_name.$file_name;
					$bath_file_name =~ s/\.[^\.]+$/\.bath\.hmm/;
					
					push(@FilesToBATHCONVERT,$family_dir_name.$file_name)
						if (!(-e $bath_file_name));
				
				}
			}
			elsif (lc($file_name) =~ /^(\S+)[\.cons]?\.a?fa[sta]?$/)
			{
				my $file_base_name = $1;

				# We don't want to build an HMM on this if it's a
				# consensus FASTA and we're not miniprot-ing
				if (!$OPTIONS{'miniprot'} && lc($file_name) =~ /\.cons\.[^\.]+$/) { next; }

				# In case runtime was killed while there was a 'target' file,
				# we don't want to treat that as a query
				if ($file_base_name !~ /\.target$/ && !(-e $family_dir_name.$file_base_name.'.bath.hmm')) 
				{
					push(@FilesToBATHBUILD,$family_dir_name.$file_name) 
				}
			}

		}
	}
	closedir($FamilyDir);



	# Build HMMs for any sequence files that didn't have one
	foreach my $file_to_bathbuild (@FilesToBATHBUILD)
	{
		
		$file_to_bathbuild =~ /^(\S+)\.([^\.]+)$/;
		my $hmm_file_name  = $1.'.bath.hmm';
		my $in_file_ext    = $2;


		my $bathbuild_cmd  = "--amino \"$hmm_file_name\" \"$file_to_bathbuild\" 1>/dev/null";

		$bathbuild_cmd = "--unali ".$bathbuild_cmd 
			if (lc($in_file_ext) !~ /^a/);

		$bathbuild_cmd = $BATHBUILD.' '.$bathbuild_cmd;
		if (system($bathbuild_cmd)) 
		{
			die "\n  ERROR:  Failed to build HMM on '$file_to_bathbuild' (command:'$bathbuild_cmd')\n\n";
		}

		push(@InputHMMs,$hmm_file_name);

	}



	# You will become one with the BATHBORG!
	foreach my $file_to_bathconvert (@FilesToBATHCONVERT)
	{

		$file_to_bathconvert   =~ /^(\S+)\.[^\.]+$/;
		my $bath_hmm_file_name = $1.'.bath.hmm';

		my $bathconvert_cmd = "$BATHCONVERT \"$bath_hmm_file_name\" \"$file_to_bathconvert\" 1>/dev/null";
		if (system($bathconvert_cmd))
		{
			die "\n  ERROR:  Failed to convert HMM file '$file_to_bathconvert' to BATH HMM (command:'$bathconvert_cmd')\n\n";
		}

		push(@InputHMMs,$bath_hmm_file_name);

	}



	# If we didn't find any HMMs (or HMM-ables, or consensus seqs) shout, scream, and cry
	if (scalar(@InputHMMs) == 0 && scalar(@InputConsSeqs))
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
	# target sequence and running bathsearch!
	# Note that in the case of PANTHER testing we're going
	# to have multiple targets (each of the genomes)
	my @Inputs;
	if ($OPTIONS{'miniprot'}) { @Inputs = @InputConsSeqs; }
	else                      { @Inputs = @InputHMMs;     }

	my @FamilySuccesses;
	my @FamilyErrors;
	foreach my $input_file_name (@Inputs) 
	{

		$input_file_name =~ /\/([^\/]+)$/;
		my $query_file_name = $1;

		$query_file_name =~ /^([^\.]+)/;
		my $query_base_name = $1;


		my @TargetFileNames;
		my @QueryIDs;
		if ($OPTIONS{'panther'})
		{
			foreach my $species (keys %SPECIES_TO_GENOME)
			{
				my $genome   = $SPECIES_TO_GENOME{$species};
				my $query_id = $query_base_name.'.'.$species;
				push(@TargetFileNames,$genome);
				push(@QueryIDs,$query_id);
			}
		}
		else
		{
			# Single target
			push(@TargetFileNames,DetermineTargetSeq($input_file_name));
			push(@QueryIDs,$query_base_name);
		}


		for (my $target_id=0; $target_id<scalar(@TargetFileNames); $target_id++)
		{

			my $target_file_name = $TargetFileNames[$target_id];
			my $query_id         = $QueryIDs[$target_id];
			my $out_file_name    = $fam_out_dir_name.$query_id.'.out';
			my $err_file_name    = $fam_out_dir_name.$query_id.'.err';

			my $bathsearch_cmd = "/usr/bin/time -v $BATHSEARCH --cpu $OPTIONS{'threads-per-job'}";
			$bathsearch_cmd = $bathsearch_cmd." --qformat fasta" if ($OPTIONS{'miniprot'});
			$bathsearch_cmd = $bathsearch_cmd." -o $out_file_name $input_file_name $target_file_name 2>$err_file_name";
			
			if (system($bathsearch_cmd)) 
			{

				push(@FamilyErrors,$query_id);
				system("echo \"\% $bathsearch_cmd\" >> $ERROR_FILE");

				if ($OPTIONS{'err-kills'})
				{
					unless ($GENOME_LIST{$target_file_name}) {
						system("rm \"$target_file_name\"");
						system("rm \"$target_file_name.ssi\"");
					}
					die "\n  ERROR: Failed during execution of command: $bathsearch_cmd\n\n";
				}

			} 
			else 
			{
				push(@FamilySuccesses,$query_id);
			}


			# If we created a target sequence file, kill it!
			unless ($GENOME_LIST{$target_file_name}) 
			{
				system("rm \"$target_file_name\"");
				system("rm \"$target_file_name.ssi\"");
			}

		}



		# If we're not running miniprot, we're done with this family!
		next if (!$OPTIONS{'miniprot'});


		# MINIPROT TIME!
		for (my $target_id=0; $target_id<scalar(@TargetFileNames); $target_id++)
		{

			my $mp_out_dir_name = $fam_out_dir_name.'miniprot/';
			die "\n  ERROR:  Failed to create miniprot output directory '$mp_out_dir_name'\n\n"
				if (system("mkdir \"$mp_out_dir_name\""));

			my $target_file_name = $TargetFileNames[$target_id];
			$target_file_name    =~ s/\.[^\.]+$/\.mpi/;
			my $query_id         = $QueryIDs[$target_id];
			my $out_file_name    = $mp_out_dir_name.$query_id.'.out';
			my $err_file_name    = $mp_out_dir_name.$query_id.'.err';

			my $miniprot_cmd = "/usr/bin/time -v $MINIPROT --aln --trans -t$OPTIONS{'threads-per-job'} $target_file_name $input_file_name 1>$out_file_name 2>$err_file_name";

			system("echo \"\% $miniprot_cmd\" >> $ERROR_FILE")
				if (system($miniprot_cmd));

		}

	}


	# Compile some basic metadata about the results of
	# splashing around
	if ($OPTIONS{'panther'})
	{
		foreach my $species (keys %SPECIES_TO_GENOME)
		{
			CompileBasicResults($fam_out_dir_name,\@FamilySuccesses,$species);
		}
	}
	else
	{
		CompileBasicResults($fam_out_dir_name,\@FamilySuccesses,0);		
	}


	# If this is being run as part of a larger test,
	# we'll want to collect the number of error-full and
	# error-free runs.
	return ($fam_out_dir_name,\@FamilySuccesses,\@FamilyErrors);


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
		next if ($family =~ /^\./);

		$family =~ s/\/$//;
		my $family_dir_name = $protein_meta_dir_name.$family.'/';
		next if (!(-d $family_dir_name));

		# Confirm that there's at least one sequence in this directory
		opendir(my $FamilyDir,$family_dir_name)
			|| die "\n  ERROR:  Failed to open family directory '$family_dir_name' (looking to confirm at least one sequence)\n\n";

		my $contains_seqs = 0;
		while (my $file_check = readdir($FamilyDir)) 
		{
			if (lc($file_check) =~ /\.fa[sta]?$/ || lc($file_check) =~ /\.hmm$/)
			{
				$contains_seqs = 1;
				last;
			}
		}
		closedir($FamilyDir);

		push(@FamilyDirs,$family_dir_name) if ($contains_seqs);

	}
	closedir($ProteinMetaDir);


	# We need to make sure that we have consistency in our accesses of
	# family directories across runs / tasks, so sort 'em!
	@FamilyDirs  = sort @FamilyDirs;
	my $num_fams = scalar(@FamilyDirs);


	# Make a directory to hold the results
	# (unless we're slurmin', in which case it already exists!)
	my $rbg_dir_name = $OPTIONS{'output-dir'}.'Results-by-Gene/';
	if (!$OPTIONS{'slurm'})
	{
		die "\n  ERROR:  Failed to create 'Results-by-Gene' directory ($rbg_dir_name)\n\n"
			if (system("mkdir \"$rbg_dir_name\""));
	}



	# The way we handle process splitting depends on whether or not
	# we're slurming
	my $num_cpus = $OPTIONS{'num-cpus'};
	$num_cpus    = $num_fams if ($num_cpus > $num_fams);


	my $master_thread  = 1;
	my $thread_id      = 0;
	my $thread_portion = int($num_fams/$num_cpus);


	my $start_fam_id;
	my $end_fam_id;
	if ($OPTIONS{'slurm'})
	{

		# If we're using slurm, each task gets split into 4 threads,
		# so we'll need to adjust our threads per job and thread portion
		$OPTIONS{'threads-per-job'} = 2;
		$thread_portion = int($num_fams/(4 * $num_cpus));

		my $task_id = $OPTIONS{'slurm'} - 1;

		# This task gets four threads.  What a treat!
		# But where do they *generally* start?
		my $pid = 0;
		my $active_threads = 1;
		while ($active_threads < 4) 
		{
			if ($pid = fork)
			{
				die "\n  ERROR: Fork failed\n\n" if (not defined $pid);
				$active_threads++;
			}
			else
			{
				$thread_id = (4 * $task_id) + $active_threads;
				$master_thread = 0;
				last;
			}
		}

		$thread_id = 4 * $task_id if ($master_thread);

		$start_fam_id =  $thread_id    * $thread_portion;
		$end_fam_id   = ($thread_id+1) * $thread_portion;
		$end_fam_id   = $num_fams if ($thread_id+1 == 4*$num_cpus)

	}
	else
	{

		# This is just our run-of-the-mill process splitting
		my $pid = 0;
		my $active_threads = 1;
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
				$master_thread = 0;
				last;
			}
		}

		$start_fam_id =  $thread_id    * $thread_portion;
		$end_fam_id   = ($thread_id+1) * $thread_portion;
		$end_fam_id   = $num_fams if ($thread_id == $num_cpus-1);

	}


	# I don't think the system calls that write to
	# the error file would race, but let's play it safe	 
	$ERROR_FILE =~ s/\.err$/\.$thread_id\.err/ if ($thread_id);


	# DEBUGGING
	#my $check_in_str = "Thread $thread_id ready to tackle families $start_fam_id..$end_fam_id-1";
	#system("echo \"$check_in_str\" >> $ERROR_FILE");


	# Iterate over the families and splash 'em up!
	for (my $fam_id = $start_fam_id; $fam_id < $end_fam_id; $fam_id++) 
	{

		my $fam_in_dir_name = $FamilyDirs[$fam_id];

		$fam_in_dir_name =~ /\/([^\/]+)\/$/;
		my $family = $1;

		# NOTE: The 'successes' and 'errors' are the base names of the files
		#       (e.g., 'file.out' minus '.out').  Currently, we aren't really
		#       using these.  But we could.  And that's what matters.
		#
		my ($fam_out_dir_name,$fam_successes_ref,$fam_errors_ref) 
			= FamilySplash($fam_in_dir_name);

		if ($fam_out_dir_name && -d $fam_out_dir_name)
		{
			my $rbg_fam_dir_name = $rbg_dir_name.$family.'/';
			die "\n  ERROR:  Failed to move family directory '$fam_out_dir_name' to '$rbg_fam_dir_name'\n\n"
				if (system("mv \"$fam_out_dir_name\" \"$rbg_fam_dir_name\""));
		}
		else
		{
			my $fam_err_str = "Error: No output directory seems to have been created for family '$family'";
			system("echo \"$fam_err_str\" >> $ERROR_FILE");
		}

	}


	# End of the real work!
	exit(0) if (!$master_thread);
	while (wait() != -1) {}


	# If we weren't slurming, thread 0 can start aggregating
	# data across all of our processes!
	if (!$OPTIONS{'slurm'})
	{

		# Pull in the error output from each of the helpers
		for ($thread_id=1; $thread_id<$num_cpus; $thread_id++) 
		{
			my $thread_err_file = $ERROR_FILE;
			$thread_err_file =~ s/\.err$/\.$thread_id\.err/;

			if (-e $thread_err_file)
			{
				system("cat \"$thread_err_file\" >> \"$ERROR_FILE\"");
				system("rm \"$thread_err_file\"");
			}
		}


		# Now that the ugliness is over, the FUN!!!
		# Combine all of our family-specific output data into
		# one nice big output file.
		#
		# If we're doing a PANTHER run, we'll want to break
		# these up across species.
		#
		if ($OPTIONS{'panther'})
		{
			foreach my $species (keys %SPECIES_TO_GENOME)
			{
				AggregateAllResults($rbg_dir_name,$species);
			}
		}
		else
		{
			AggregateAllResults($rbg_dir_name,0);
		}

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
	$hmm_file_name =~ /^(.*\/)([^\/]+)\.bath\.hmm$/;
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
		$target_file_name =~ s/\.bath\.hmm$/\.target\.fa/;

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
	print "  USAGE: ./Run-Splash-Test.pl {OPT.S} [Gene-Super-Directory/] [Species-Guide.txt]\n";
	print "\n";
	print "  OPT.S: --full-genome       : Force use of full genome as target sequence.\n";
	print "         --panther           : Assume that input HMMs are from PANTHER and that\n";
	print "                                we're searching each HMM against all genomes.\n";
	print "         --miniprot          : Run miniprot alongside Splash.\n";
	print "                                (This requires having pre-run 'GenPANTHERConsSeqs.pl')\n";
	print "         --summarize   [str] : Given the name of a SLURM-based output directory,\n";
	print "                                run 'AggregateAllResults' (Still requires input info!)\n";
	print "         --err-kills         : If a bathsearch run fails, kill the script\n";
	print "                                (by default, the error is simply logged).\n";
	print "         --cpus/-n     [int] : Set the number of threads to use (default:1).\n";
	print "         --slurm/-s    [int] : Use srun to guide bathsearch (and miniprot) execution.\n";
	print "                                The provided integer is the job ID, from [1..cpus],\n";
	print "         --out-dir/-o  [str] : Name the output directory (default:'Splash-Results').\n";
	print "         --gene/-g     [str] : Run on a single gene within the input directory.\n";
	die   "\n\n";
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
			if ($line =~ /^=\s+(\S+)\s+(\d+)\s*/) {
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

	# SET THE DEFAULTS!

	# Output directory name
	$OPTIONS{'output-dir'} = 'Splash-Results';

	# How many CPUs do we want? (this should probably be communicated as 'jobs' or something...)
	$OPTIONS{'num-cpus'} = 1;

	# Are we also running miniprot?
	$OPTIONS{'miniprot'} = 0;

	# Are we using SLURM?  If so, use an 'srun' with supporting data.
	$OPTIONS{'slurm'} = 0;

	# Are we aggregating output data from a SLURM-based run?
	$OPTIONS{'post-slurm'} = 0;

	# How much sequence do we want to pull in around any tightly
	# defined nucleotide ranges?
	$OPTIONS{'num-extra-nucls'} = 5000;

	# Do we terminate this script's execution if an error occurs?
	$OPTIONS{'err-kills'} = 0;

	# Not really an option, but we'll call it one...
	# What does '~' mean?
	$OPTIONS{'tilde-dir'} = GetTildeDirName();


	# Time for non-stop fun! (parsing commandline arguments)
	my $num_args = scalar(@ARGV);
	for (my $arg_id=0; $arg_id<$num_args; $arg_id++) {

		my $Arg = $ARGV[$arg_id];

		if (lc($Arg) =~ /^-?-?full-genome$/) 
		{
			$OPTIONS{'full-genome'} = 1; # Overrides existence of '.genome-range.out' files
		}
		elsif (lc($Arg =~ /^-?-?panther$/))
		{
			$OPTIONS{'panther'} = 1; # PANTHER time, baby!
		}
		elsif (lc($Arg =~ /^-?-?miniprot$/))
		{
			ConfirmMiniprot();
			$OPTIONS{'miniprot'} = 1;
			$OPTIONS{'panther'}  = 1; # Miniprot forces PANTHER
		}
		elsif (lc($Arg =~ /^-?-?slurm$|^-s$/))
		{
			$arg_id++;
			my $task_id = $ARGV[$arg_id];
			die "\n  ERROR:  Task ID $task_id does not appear to be a positive integer...\n\n"
				if ($task_id =~ /\D/);
			$OPTIONS{'slurm'} = int($task_id); # Use srun on each system call

		}
		elsif (lc($Arg =~ /^-?-?summarize$/))
		{
			$arg_id++;
			my $result_dir_name = $ARGV[$arg_id];
			die "\n  ERROR:  Failed to locate result directory '$result_dir_name'\n\n"
				if (!(-e $result_dir_name));
			$result_dir_name = $result_dir_name.'/' if ($result_dir_name !~ /\/$/);
			$OPTIONS{'post-slurm'} = $result_dir_name;
		}
		elsif (lc($Arg =~ /^-?-?err-kills$/))
		{
			$OPTIONS{'err-kills'} = 1; # First failed run terminates script
		}
		elsif (lc($Arg) =~ /^-?-?cpus$/ || lc($Arg) =~ /^-?n$/)
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
		elsif (lc($Arg) =~ /^-?-?out-dir$/ || lc($Arg) =~ /^-?o$/)
		{
			$arg_id++;
			$OPTIONS{'output-dir'} = $ARGV[$arg_id];
		}
		elsif (lc($Arg) =~ /^-?-?gene$/ || lc($Arg) =~ /^-?g$/)
		{
			$arg_id++;
			$OPTIONS{'single-gene'} = lc($ARGV[$arg_id]);
		}
		elsif ($arg_id == $num_args-2)
		{

			if (!(-d $Arg))
			{
				die "\n  ERROR:  Failed to make sense of (inferred) query instruction '$Arg'\n\n";
			}


			# Whatever's going on, this had better be a directory!
			$Arg = $Arg.'/' if ($Arg !~ /\/$/);
			die "\n  ERROR:  Failed to locate gene super-directory '$Arg'\n\n"
				if (!(-d $Arg));

			if ($OPTIONS{'single-gene'})
			{
				$OPTIONS{'protein-input'} = $Arg.$OPTIONS{'single-gene'}.'/';
				die "\n  ERROR:  No directory for gene '$OPTIONS{'single-gene'}' found in directory '$Arg'\n\n"
					if (!(-d $OPTIONS{'protein-input'}));
			}
			else
			{
				$OPTIONS{'protein-input'} = $Arg;
			}

		}
		elsif ($arg_id == $num_args-1)
		{
			$OPTIONS{'genome-map'} = ValidateSpeciesGuide($Arg);
			BuildSpeciesToGenomeMap();
		}
		else
		{
			die "\n  ERROR:  Unrecognized (assumed) protein input argument: '$Arg'\n\n";
		}

	}



	# Determine the actual name of the output directory
	# If we're using slurm, the output directory and the
	# 'Results-by-Gene' subdirectory should be available
	# ahead of time.
	if ($OPTIONS{'post-slurm'})
	{
		RunPostSlurmAggregation();
		exit(0);
	}
	elsif (!$OPTIONS{'slurm'})
	{

		my $base_out_dir_name = $OPTIONS{'output-dir'};
		$base_out_dir_name =~ s/\/$//;
	
		my $out_dir_name = $base_out_dir_name;
		my $attempt = 1;
		while (-d $out_dir_name) {
			$attempt++;
			$out_dir_name = $base_out_dir_name.'-'.$attempt;
		}
	
		$OPTIONS{'output-dir'} = $out_dir_name.'/';
	
	}
	else
	{

		my $out_dir_name = $OPTIONS{'output-dir'};
		$out_dir_name = $out_dir_name.'/' if ($out_dir_name !~ /\/$/);

		die "\n  ERROR:  Slurm-ified running requires that the output directory already exists (failed to locate '$out_dir_name')\n\n"
			if (!(-d $out_dir_name));

		my $rbg_dir_name = $out_dir_name.'Results-by-Gene/';
		die "\n  ERROR:  Slurm-ified running requires that the Results-by-Gene directory already exists (failed to locate '$rbg_dir_name')\n\n"
			if (!(-d $rbg_dir_name));

		$OPTIONS{'output-dir'} = $out_dir_name;

	}

}







###########################################################################
#
#  Subroutine: GetTildeDirName
#
sub GetTildeDirName
{
	open(my $TildeCheck,'echo ~ |') 
		|| die "\n  ERROR:  Failed to run tilde dir check\n\n";

	my $tilde_dir_name = <$TildeCheck>;
	$tilde_dir_name =~ s/\n|\r//g;
	$tilde_dir_name = $tilde_dir_name.'/' if ($tilde_dir_name !~ /\/$/);

	close($TildeCheck);

	return $tilde_dir_name;
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
#  Subroutine: ValidateSpeciesGuide
#
sub ValidateSpeciesGuide
{

	my $species_guide_name = shift;

	open(my $SpeciesGuide,'<',$species_guide_name)
		|| die "\n  ERROR:  Failed to open species guide '$species_guide_name'\n\n";
	
	my @SpeciesGenomeMaps;
	while (my $line = <$SpeciesGuide>) {
		if ($line =~ /^\s*(\S+)\s+(\S+)\s/) {
			my $species = $1;
			my $genome  = ConfirmGenome($2);
			push(@SpeciesGenomeMaps,$species.'|'.$genome);
		}
	}
	close($SpeciesGuide);


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


	if ($genome =~ /^\~\//) 
	{
		$genome =~ s/^\~\///;
		$genome = $OPTIONS{'tilde-dir'}.$genome;
	} 
	elsif ($genome =~ /^\$HOME\//)
	{
		$genome =~ s/^\$HOME\///;
		$genome = $OPTIONS{'tilde-dir'}.$genome;
	}


	if (!(-e $genome)) 
	{
		die "\n  ERROR:  Failed to locate genomic file '$genome'\n\n";
	}

	my $ssi_file_name = $genome.'.ssi';
	if (!(-e $ssi_file_name) && system("$SFETCH --index \"$genome\"")) 
	{
		die "\n  ERROR:  Failed to produce an ssi index for genomic file '$genome'\n\n";
	}

	my $seqstat_file_name = $genome.'.seqstat.out';
	if (!(-e $seqstat_file_name) && system("$SEQSTAT -a \"$genome\" > \"$seqstat_file_name\"")) 
	{
		die "\n  ERROR:  Failed to aggregate sequence statistics for genomic file '$genome'\n\n";
	}


	# If we're using miniprot, be sure we have the '.mpi' index
	# pre-computed
	my $mpi_file_name = $genome;
	$mpi_file_name =~ s/\.[^\.]+$/\.mpi/;
	if ($OPTIONS{'miniprot'}) 
	{
		if (!(-e $mpi_file_name) && system("$MINIPROT -t2 -d \"$mpi_file_name\" \"$genome\""))
		{
			die "\n  ERROR:  Failed to produce an mpi index for genomic file '$genome'\n\n";
		}
	}

	return $genome; # Looks like a genome to me!

}








###########################################################################
#
#  Subroutine: ConfirmRequiredTools
#
sub ConfirmRequiredTools
{

	my $bath_base_dir = 'BATH/';
	die "\n  ERROR:  Failed to locate base BATH directory '$bath_base_dir'\n          Try running './Build.pl'\n\n"
		if (!(-d $bath_base_dir));


	my $easel_dir = $bath_base_dir.'easel/';
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


	my $bath_src_dir = $bath_base_dir.'src/';
	die "\n  ERROR:  Failed to locate BATH source directory '$bath_src_dir'\n\n"
		if (!(-d $bath_src_dir));


	$BATHBUILD = $bath_src_dir.'bathbuild';
	die "\n  ERROR:  Failed to locate bathbuild in BATH source directory (looking for '$BATHBUILD')\n\n"
		if (!(-x $BATHBUILD));


	$BATHCONVERT = $bath_src_dir.'bathconvert';
	die "\n  ERROR:  Failed to locate bathconvert in BATH source directory (looking for '$BATHCONVERT')\n\n"
		if (!(-x $BATHCONVERT));


	$BATHSEARCH = $bath_src_dir.'bathsearch';
	die "\n  ERROR:  Failed to locate bathhsearch in BATH source directory (looking for '$BATHSEARCH')\n\n"
		if (!(-x $BATHSEARCH));


}








###########################################################################
#
#  Subroutine: AnticipateMiniprot
#
sub AnticipateMiniprot
{
	my $miniprot = $0;
	$miniprot =~ s/[^\/]+$/miniprot\/miniprot/;
	return $miniprot;
}

###########################################################################
#
#  Subroutine: ConfirmMiniprot
#
sub ConfirmMiniprot
{
	if (!(-e $MINIPROT))
	{
		print "\n";
		print "  ERROR:  Failed to locate miniprot executable (looking for '$MINIPROT')\n";
		print "          Try: ./Build.pl --miniprot\n";
		die   "\n";
	}
}









###########################################################################
#
#  Subroutine: GetPct
#
sub GetPct
{
	my $numerator   = shift;
	my $denominator = shift;
	my $pct = int(1000.0 * $numerator / $denominator) / 10.0;
	return $pct;
}

