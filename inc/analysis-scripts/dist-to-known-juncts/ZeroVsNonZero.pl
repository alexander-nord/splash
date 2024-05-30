#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV != 1) { die "\n  USAGE:  ./ZeroVsNonZero.pl [junction-nucl-diff-counts.csv]\n\n"; }


open(my $CSV,'<',$ARGV[0])
	|| die "\n  ERROR:  Failed to open input file '$ARGV[0]'\n\n";

my $zero       = 0;
my $one_or_two = 0;
my $three_plus = 0;

while (my $line = <$CSV>)
{
	if ($line =~ /^(\d+),(\d+)/)
	{
		my $dist  = $1;
		my $count = $2;

		if    ($dist == 0) {    $zero     = $count; }
		elsif ($dist == 1) { $one_or_two += $count; }
		elsif ($dist == 2) { $one_or_two += $count; }
		else               { $three_plus += $count; }
	}
}

close($CSV);


my $num_juncts = $zero + $one_or_two + $three_plus;


my $zero_pct       = int(10000.0 * $zero       / $num_juncts) / 100.0;
my $one_or_two_pct = int(10000.0 * $one_or_two / $num_juncts) / 100.0;
my $three_plus_pct = int(10000.0 * $three_plus / $num_juncts) / 100.0;


# Assuming 0 is the majority
$one_or_two = ' '.$one_or_two while (length("$one_or_two") < length("$zero"));
$three_plus = ' '.$three_plus while (length("$three_plus") < length("$zero"));


print "\n";
print "  Zero := $zero => $zero_pct\%\n";
print "  1,2  := $one_or_two => $one_or_two_pct\%\n";
print "  3+   := $three_plus =>  $three_plus_pct\%\n";
print "\n";


1;
