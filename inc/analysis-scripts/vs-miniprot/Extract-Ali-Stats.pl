#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub GetFamConsSeqLen;
sub GatherSplashStats;
sub GatherMiniprotStats;
sub GetTimeAndMaxRSS;


if (@ARGV != 1) { die "\n  USAGE:  ./Extract-Ali-Stats.pl [MP-vs-Splash-Output]\n\n"; }


my $result_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate input directory '$result_dir_name'\n\n"
	if (!(-d $result_dir_name));
$result_dir_name = $result_dir_name.'/' if ($result_dir_name !~ /\/$/);


my $rbg_dir_name = $result_dir_name.'Results-by-Gene/';
die "\n  ERROR:  Failed to locate 'Results-by-Gene' directory '$rbg_dir_name'\n\n"
	if (!(-d $rbg_dir_name));



opendir(my $RBG,$rbg_dir_name)
	|| die "\n  ERROR:  Failed to open directory '$rbg_dir_name'\n\n";

my %FamToDirName;
my %FamToSeqLen;
while (my $fam = readdir($RBG))
{
	next if ($fam =~ /^\./);
	$fam =~ s/\/$//;

	my $fam_dir_name = $rbg_dir_name.$fam.'/';
	next if (!(-d $fam_dir_name.'miniprot/'));

	$FamToDirName{$fam} = $fam_dir_name;
	$FamToSeqLen{$fam}  = GetFamConsSeqLen($fam_dir_name);

}

closedir($RBG);


my $num_fams = scalar(keys %FamToDirName);
die "\n  ERROR:  Failed to identify any suitable families...\n\n"
	if ($num_fams == 0);


# Hard-coding is best coding!
my @Species;
push(@Species,'human');
push(@Species,'chicken');
push(@Species,'stickleback');
push(@Species,'fruit_fly');


foreach my $species (@Species)
{

	my $species_out_file_name = $result_dir_name.$species.'-comparison.csv';
	open(my $SpeciesCSV,'>',$species_out_file_name)
		|| die "\n  ERROR:  Failed to open species output file '$species_out_file_name'\n\n";


	print $SpeciesCSV "Family,Seq-Length,";
	
	print $SpeciesCSV "Splash-Time(S),Splash-Max-Mem(MB),Splash-Chr,";
	print $SpeciesCSV "Splash-Ali-Start,Splash-Ali-End,Splash-Peak-Coverage,";
	print $SpeciesCSV "Splash-Num-Exons,Splash-Pct-ID,";
	
	print $SpeciesCSV "MP-Time(S),MP-Max-Mem(MB),MP-Chr,";
	print $SpeciesCSV "MP-Ali-Start,MP-Ali-End,MP-Peak-Coverage,";
	print $SpeciesCSV "MP-Num-Exons,MP-Pct-ID\n";
	
	print $SpeciesCSV "\n";


	my $only_splash  = 0;
	my $only_mp      = 0;
	my $both_succeed = 0;
	my $both_fail    = 0;

	foreach my $fam (sort keys %FamToDirName)
	{

		my $fam_dir_name = $FamToDirName{$fam};
		my $fam_seq_len  = $FamToSeqLen{$fam};


		my $splash_out_file_name = $fam_dir_name.$fam.'.'.$species.'.out';
		my $splash_err_file_name = $fam_dir_name.$fam.'.'.$species.'.err';

		my $mp_out_file_name = $fam_dir_name.'miniprot/'.$fam.'.'.$species.'.out';
		my $mp_err_file_name = $fam_dir_name.'miniprot/'.$fam.'.'.$species.'.err';


		my $cannot_compare = 0;
		if (!(-e $splash_out_file_name && $splash_err_file_name))
		{
			print "  WARNING: Missing splash out/err file for $species $fam\n";
			$cannot_compare = 1;
		}
		if (!(-e $mp_out_file_name && $mp_err_file_name))
		{
			print "  WARNING: Missing miniprot out/err file for $species $fam\n";
			$cannot_compare = 1;
		}
		next if ($cannot_compare);


		print $SpeciesCSV "$fam,$fam_seq_len,";


		my $splash_success = GatherSplashStats($splash_out_file_name,$splash_err_file_name,$SpeciesCSV,$fam_seq_len);
		my $mp_success = GatherMiniprotStats($mp_out_file_name,$mp_err_file_name,$SpeciesCSV,$fam_seq_len);


		if ($splash_success) 
		{
			if ($mp_success) { $both_succeed++; }
			else             { $only_splash++;  }
		}
		elsif ($mp_success)  { $only_mp++;      }
		else                 { $both_fail++;    }

	}

	close($SpeciesCSV);


	while (length("$both_succeed") < length("$num_fams")) { $both_succeed = ' '.$both_succeed; }
	while (length("$only_splash")  < length("$num_fams")) { $only_splash  = ' '.$only_splash;  }
	while (length("$only_mp")      < length("$num_fams")) { $only_mp      = ' '.$only_mp;      }
	while (length("$both_fail")    < length("$num_fams")) { $both_fail    = ' '.$both_fail;    }

	print "\n";
	print "$species\n";
	print "  Both Succeed  : $both_succeed\n";
	print "  Splash Only   : $only_splash\n";
	print "  MiniProt Only : $only_mp\n";
	print "  Both Fail     : $both_fail\n";
	print "\n";

}



1;






######################################################################################
#
#  Subroutine:  GetFamConsSeqLen
#
sub GetFamConsSeqLen
{
	my $fam_dir_name = shift;

	opendir(my $FamDir,$fam_dir_name)
		|| die "\n  ERROR:  Failed to open directory '$fam_dir_name'\n\n";

	my $seq_len = 0;
	while (my $file_name = readdir($FamDir))
	{

		next if ($file_name !~ /\.out$/);
		$file_name = $fam_dir_name.$file_name;

		open(my $HitFile,'<',$file_name)
			|| die "\n  ERROR:  Failed to open file '$file_name'\n\n";

		while (my $line = <$HitFile>)
		{
			if ($line =~ /\[M=(\d+)\]/) 
			{
				$seq_len = $1;
				last;
			}
		}

		close($HitFile);

		last if ($seq_len);
	
	}

	closedir($FamDir);

	return $seq_len;

}






######################################################################################
#
#  Subroutine:  GatherSplashStats
#
sub GatherSplashStats
{
	
	my $splash_out_file_name = shift;
	my $splash_err_file_name = shift;
	my $SpeciesCSV = shift;
	my $fam_seq_len = shift;


	my ($time,$rss) = GetTimeAndMaxRSS($splash_err_file_name);


	open(my $OutFile,'<',$splash_out_file_name)
		|| die "\n  ERROR:  Failed to open Splash output file '$splash_out_file_name'\n\n";

	my $num_exons  = '-';
	my $ali_start  = '-';
	my $ali_end    = '-';
	my $coverage   = '-';
	my $chr        = '-';
	my $pct_id     = '-';
	my $matches    =  0 ;
	my $mismatches =  0 ;

	# There are ~5 cases where it appears the program died,
	# so we need to catch those...
	my $peak_coverage = 0.0;

	while (my $line = <$OutFile>)
	{

		next if ($line !~ /\| = Exon Set \d+ \((\d+)/);
		my $es_num_exons = $1;

		$line = <$OutFile>;
		$line =~ /\| = Model Positions (\d+)\.\.(\d+)/;
		my $es_ali_start = $1;
		my $es_ali_end   = $2;

		my $es_coverage = int(1000.0 * ($es_ali_end - $es_ali_start + 1) / $fam_seq_len) / 10.0;

		# Is this the best exon set we've seen so far?
		next if ($es_coverage < $peak_coverage);

		$peak_coverage = $es_coverage;

		$num_exons = $es_num_exons;
		$ali_start = $es_ali_start;
		$ali_end   = $es_ali_end;
		$coverage  = $es_coverage.'%';


		$line = <$OutFile>;
		$line =~ /\| = Target Seq Name (\S+)/;
		$chr = $1;

		$line = <$OutFile> while ($line !~ /^:/);
		<$OutFile>;

		for (my $exon_id = 1; $exon_id <= $num_exons; $exon_id++)
		{
			
			<$OutFile>; # Exon ID Line;
			<$OutFile>; # Blank

			my $query_line = <$OutFile>; # Query
			if (!$query_line)
			{
				close($OutFile);
				print $SpeciesCSV "$time,$rss,-,-,-,-,-,-,";
				$splash_out_file_name =~ /\/([^\/]+)$/;
				print STDERR " [ SPLASH ERROR ]--> $1\n";
				return 0;
			}
			$query_line =~ s/\n|\r//g;

			while ($query_line && $query_line !~ /^:/)
			{

				<$OutFile>; # Quality line

				my $trans_line = <$OutFile>;
				if (!$trans_line)
				{
					close($OutFile);
					print $SpeciesCSV "$time,$rss,-,-,-,-,-,-,";
					$splash_out_file_name =~ /\/([^\/]+)$/;
					print STDERR " [ SPLASH ERROR ]--> $1\n";
					return 0;
				}
				$trans_line =~ s/\n|\r//g;

				my @Q = split(//,$query_line);
				my @T = split(//,$trans_line);

				# Walk up to the query sequence name
				my $pos = 0;
				$pos++ while ($Q[$pos] eq ' ');

				# Walk through the query sequence name
				$pos++ while ($Q[$pos] =~ /\S/);

				# NOW we can get countin'!
				# Note that our checks for whitespace should make this
				#   valid, even given that the query line will have indices
				while ($pos < scalar(@T))
				{

					my $q = uc($Q[$pos]);
					my $t = uc($T[$pos]);
					$pos++;

					# Whitespace skip
					next if ($q eq ' ');
					next if ($t eq ' ');

					# Gap skip
					next if ($q eq '.');
					next if ($t eq '-');

					if ($q eq $t) {    $matches++; }
					else          { $mismatches++; }

				}

				<$OutFile>; # Nucleotides
				<$OutFile>; # PP LINE!
				<$OutFile>; # Blank

				$query_line = <$OutFile>;
				if (!$query_line)
				{
					close($OutFile);
					$splash_out_file_name =~ /\/([^\/]+)$/;
					print STDERR " [ SPLASH ERROR ]--> $1\n";
					return 0;
				}
				$query_line =~ s/\n|\r//g;

			}

		}

		$pct_id = int(1000.0 * $matches / ($matches + $mismatches)) / 10.0;
		$pct_id = $pct_id.'%';

	}

	close($OutFile);


	print $SpeciesCSV "$time,$rss,$chr,";
	print $SpeciesCSV "$ali_start,$ali_end,$coverage,";
	print $SpeciesCSV "$num_exons,$pct_id,";


	if ($pct_id =~ /\%/)
	{
		# In order to be a "full success" (for fairness to miniprot)
		# we need at least 10% coverage
		$coverage =~ s/\%//;
		return 1 if ($coverage * 1.0 >= 10.0);
		return 0;
	}


	return 0;

}






######################################################################################
#
#  Subroutine:  GatherMiniprotStats
#
sub GatherMiniprotStats
{
	
	my $mp_out_file_name = shift;
	my $mp_err_file_name = shift;
	my $SpeciesCSV = shift;
	my $fam_seq_len = shift;


	my ($time,$rss) = GetTimeAndMaxRSS($mp_err_file_name);


	open(my $OutFile,'<',$mp_out_file_name)
		|| die "\n  ERROR:  Failed to open MiniProt output file '$mp_out_file_name'\n\n";

	my $num_exons  = '-';
	my $ali_start  = '-';
	my $ali_end    = '-';
	my $coverage   = '-';
	my $chr        = '-';
	my $pct_id     = '-';
	my $matches    =  0 ;
	my $mismatches =  0 ;

	my $summary_line = <$OutFile>;
	if ($summary_line)
	{
		
		$summary_line =~ /^\S+\s+\d+\s+(\d+)\s+(\d+)\s+\S\s+(\S+)/;
		$ali_start = $1 + 1; # For some reason miniprot seems to do [start,end)
		$ali_end   = $2;
		$chr       = $3;

		$coverage = int(1000.0 * ($ali_end - $ali_start + 1) / $fam_seq_len) / 10.0;
		$coverage = $coverage.'%';

		<$OutFile>; # Nucl line

		my $trans_line = <$OutFile>;
		$trans_line =~ s/\n|\r//g;
		$trans_line =~ s/^\S+\s+//;

		<$OutFile>; # Quality line

		my $query_line = <$OutFile>;
		$query_line =~ s/\n|\r//g;
		$query_line =~ s/^\S+\s+//;


		my $pos = 1;
		my $in_intron = 1;
		$num_exons = 0;
		
		my @Q = split(//,$query_line);
		my @T = split(//,$trans_line);

		while ($pos < scalar(@Q))
		{

			my $q = uc($Q[$pos]);
			my $t = uc($T[$pos]);

			if ($in_intron)
			{
				if ($t ne ' ')
				{
					$in_intron = 0;
					$num_exons++;
					next; # Don't advance -- we'll re-check this position!
				}

			}
			else
			{
				if ($t eq ' ')
				{
					$in_intron = 1;
				}
				elsif ($q =~ /[A-Z]/ && $t =~ /[A-Z]/)
				{
					if ($q eq $t) {    $matches++; }
					else          { $mismatches++; }
				}

			}

			$pos++;

		}

		$pct_id = int(1000.0 * $matches / ($matches + $mismatches)) / 10.0;
		$pct_id = $pct_id.'%';

	}

	close($OutFile);


	print $SpeciesCSV "$time,$rss,$chr,";
	print $SpeciesCSV "$ali_start,$ali_end,$coverage,";
	print $SpeciesCSV "$num_exons,$pct_id\n";
	

	if ($pct_id =~ /\%/)
	{
		# In order to be a "full success" (for fairness to splash)
		# we need at least 2 exons
		return 1 if ($num_exons > 1);
		return 0;
	}
	

	return 0;

}






######################################################################################
#
#  Subroutine:  GetTimeAndMaxRSS
#
sub GetTimeAndMaxRSS
{
	my $err_file_name = shift;
	open(my $ErrFile,'<',$err_file_name)
		|| die "\n  ERROR:  Failed to open error file '$err_file_name' (GetTimeAndMaxRSS)\n\n";

	my $time = '-';
	my $rss  = '-';
	while (my $line = <$ErrFile>)
	{

		$line =~ s/\n|\r//g;

		if ($line =~ /^\s*Elapsed \(wall clock\) time/)
		{
			$line =~ /\:\s+(\S+)\s*$/;
			$time = $1;

			my $hours   = 0;
			my $minutes = 0;
			my $seconds = 0;
			if ($time =~ /^([^:]+):([^:]+)$/)
			{
				$minutes = int($1);
				$seconds =     $2 ;
			}
			elsif ($time =~ /^([^:]+):([^:]+):([^:]+)$/)
			{
				$hours   = int($1);
				$minutes = int($2);
				$seconds =     $3 ;
			}
			else
			{
				print "  WARNING:  Error file '$err_file_name' had unexpected time formatting: '$line'\n";
				$time = '-';
				next;
			}

			$minutes += 60 * $hours;

			$seconds =~ s/\.\d+$//;
			$minutes = int(100.0 * $minutes + int($seconds) / 60.0) / 100.0;

			$time = $minutes;

		}
		elsif ($line =~ /^\s*Maximum resident set size \(kbytes\): (\d+)/)
		{
			$rss = int($1 / 1024);
		}

	}

	close($ErrFile);

	return ($time,$rss);

}









