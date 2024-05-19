#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV != 3) 
{
	die "\n  USAGE:  ./ExtractPANTHERFams.pl [path/to/PANTHER/] [Fam-List.txt] [out/dir/name/]\n\n";
}


my $in_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate input PANTHER directory '$in_dir_name'\n\n"
	if (!(-d $in_dir_name));
$in_dir_name = $in_dir_name.'/' if ($in_dir_name !~ /\/$/);


my $out_dir_name = $ARGV[2];
die "\n  ERROR:  Output directory '$out_dir_name' already exists\n\n"
	if (-d $out_dir_name);
$out_dir_name = $out_dir_name.'/' if ($out_dir_name !~ /\/$/);


my $fam_list_file_name = $ARGV[1];
die "\n  ERROR:  Failed to locate family list file '$fam_list_file_name'\n\n"
	if (!(-e $fam_list_file_name));


# First thing: What fams do we want?
open(my $FamListFile,'<',$fam_list_file_name)
	|| die "\n  ERROR:  Failed to open family list file '$fam_list_file_name'\n\n";

my %Fams;
while (my $line = <$FamListFile>) 
{

	$line =~ s/\n|\r//g;
	next if (!$line);

	foreach my $fam (split(/\,/,$line))
	{
		$Fams{uc($fam)} = 1;
	}

}
close($FamListFile);


# Uhhhhhh, weird?
if (scalar(keys %Fams) == 0)
{
	die "\n  ERROR:  Failed to discern any families from family list file '$fam_list_file_name'\n\n";
}


# We've got families -- time to find a place to put 'em!
die "\n  ERROR:  Failed	to create output directory '$out_dir_name'\n\n"
	if (system("mkdir \"$out_dir_name\""));


# Great! Now, let's pull 'em into the fold!
my $InputDir;
if (!opendir($InputDir,$in_dir_name))
{
	system("rmdir \"$out_dir_name\"");
	die "\n  ERROR:  Failed to open input directory '$in_dir_name'\n\n";
}


# We'll want to know if we didn't have any matches to the provided list
my $num_fams_matched = 0;
while (my $fam = readdir($InputDir))
{

	next if ($fam =~ /^\./);

	$fam =~ s/\/$//;
	my $fam_in_dir_name = $in_dir_name.$fam.'/';
	next if (!(-d $fam_in_dir_name));

	next if (!$Fams{uc($fam)});

	# Gotcha!
	my $fam_out_dir_name = $out_dir_name.$fam.'/';
	die "\n  ERROR:  Failed to create family output directory '$fam_out_dir_name'\n\n"
		if (system("mkdir \"$fam_out_dir_name\""));

	die "\n  ERROR:  Failed to copy files from '$fam_in_dir_name' to '$fam_out_dir_name'\n\n"
		if (system("cp $fam_in_dir_name* $fam_out_dir_name"));

	$num_fams_matched++;

}
closedir($InputDir);



# If we didn't find any of the families listed in our file, that's weird
if ($num_fams_matched == 0)
{
	system("rmdir \"$out_dir_name\"");
	die "\n  Warning:  None of the families listed in '$fam_list_file_name' matched sub-directory names in '$in_dir_name'\n\n"
}


1;
