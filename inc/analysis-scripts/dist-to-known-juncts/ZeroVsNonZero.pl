#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV != 1) { die "\n  USAGE:  ./ZeroVsNonZero.pl [junction-nucl-diff-counts.csv]\n\n"; }


open(my $CSV,'<',$ARGV[0])
	|| die "\n  ERROR:  Failed to open input file '$ARGV[0]'\n\n";

my    $zero = 0;
my $nonzero = 0; # lol

while (my $line = <$CSV>)
{
	if ($line =~ /^(\d+),(\d+)/)
	{
		my $dist  = $1;
		my $count = $2;

		if ($dist == 0) {    $zero  = $count; }
		else            { $nonzero += $count; }
	}
}

close($CSV);


my $zero_pct = int(1000.0 * $zero / ($zero + $nonzero)) / 10.0;


# Am I dumb or what?
$zero = ' '.$zero while (length("$zero") < length("$nonzero"));
$nonzero = ' '.$nonzero while (length("$nonzero") < length("$zero"));


print "\n";
print " Zero := $zero  ==>  ($zero_pct\%)\n";
print "!Zero := $nonzero\n";
print "\n";


1;
