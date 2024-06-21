#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV != 1) { die "\n  USAGE:  ./Reduce-to-Diffs.pl [MP-vs-Splash-Output]\n\n"; }


my $results_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate results directory '$results_dir_name'\n\n"
	if (!(-d $results_dir_name));
$results_dir_name = $results_dir_name.'/' if ($results_dir_name !~ /\/$/);


opendir(my $ResultsDir,$results_dir_name)
	|| die "\n  ERROR:  Failed to open results directory '$results_dir_name'\n\n";

while (my $csv_name = readdir($ResultsDir))
{

	next if ($csv_name !~ /^(\S+)-comparison\.csv$/);
	my $species = $1;

	$csv_name = $results_dir_name.$csv_name;

	open(my $InCSV,'<',$csv_name)
		|| die "\n  ERROR:  Failed to open input file '$csv_name'\n\n";

	my $out_csv_name = $csv_name;
	$out_csv_name =~ s/-comparison\.csv$/-differences\.csv/;

	open(my $OutCSV,'>',$out_csv_name)
		|| die "\n  ERROR:  Failed to open output file '$out_csv_name'\n\n";


	my $header = <$InCSV>;
	$header =~ s/\n|\r//g;

	my %ElementToColumn;

	my $col_id = 0;
	foreach my $element (split(/\,/,$header))
	{
		$col_id++;
		if ($element eq 'Family')
		{
			$ElementToColumn{$element} = $col_id;
		}
		elsif ($element eq 'Splash-Peak-Coverage')
		{
			$ElementToColumn{$element} = $col_id;
		}
		elsif ($element eq 'Splash-Pct-ID')
		{
			$ElementToColumn{$element} = $col_id;
		}
		elsif ($element eq 'Splash-Num-Exons')
		{
			$ElementToColumn{$element} = $col_id;
		}
		elsif ($element eq 'MP-Peak-Coverage')
		{
			$ElementToColumn{$element} = $col_id;
		}
		elsif ($element eq 'MP-Pct-ID')
		{
			$ElementToColumn{$element} = $col_id;
		}
		elsif ($element eq 'MP-Num-Exons')
		{
			$ElementToColumn{$element} = $col_id;
		}

	}

	print $OutCSV "Family,Peak-Coverage-Diff,Pct-ID-Diff\n";

	while (my $line = <$InCSV>)
	{
		$line =~ s/\n|\r//g;
		next if (!$line);

		my @Elements = split(/\,/,$line);
		
		my $family        = $Elements[$ElementToColumn{'Family'}               - 1];

		my $spl_coverage  = $Elements[$ElementToColumn{'Splash-Peak-Coverage'} - 1];
		my $spl_pct_id    = $Elements[$ElementToColumn{'Splash-Pct-ID'}        - 1];
		my $spl_num_exons = $Elements[$ElementToColumn{'Splash-Num-Exons'}     - 1];

		my $mp_coverage   = $Elements[$ElementToColumn{'MP-Peak-Coverage'}     - 1];
		my $mp_pct_id     = $Elements[$ElementToColumn{'MP-Pct-ID'}            - 1];
		my $mp_num_exons  = $Elements[$ElementToColumn{'MP-Num-Exons'}         - 1];


		if (!$mp_coverage)
		{
			die "$species / $family : $line\n";
		}

		if ($spl_coverage eq '-' || $mp_coverage eq '-') { next; }

		$spl_coverage =~ s/\%//;
		$mp_coverage  =~ s/\%//;


		if ($spl_coverage < 10.0 || $mp_coverage < 10.0) { next; }

		if ($spl_num_exons < 2   || $mp_num_exons < 2  ) { next; }


		my $coverage_diff = int(10.0 * ($spl_coverage - $mp_coverage))/10.0;
		$coverage_diff = $coverage_diff.'%';

		$spl_pct_id =~ s/\%//;
		$mp_pct_id  =~ s/\%//;
		my $pct_id_diff = int(10.0 * ($spl_pct_id - $mp_pct_id))/10.0;
		$pct_id_diff = $pct_id_diff.'%';

		print $OutCSV "$family,$coverage_diff,$pct_id_diff\n";

	}

	close($InCSV);
	close($OutCSV);

}

closedir($ResultsDir);



1;

