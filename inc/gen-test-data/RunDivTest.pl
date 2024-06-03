#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


# UNSAFE
my $BATHSEARCH = $0;
$BATHSEARCH =~ s/inc\/gen-test-data\/\S+$/BATH\/src\/bathsearch/;
die "\n  ERROR:  Failed to locate bathsearch executable '$BATHSEARCH'\n\n"
	if (!(-e $BATHSEARCH));



if (@ARGV != 4) { die "\n  USAGE:  ./RunDivTest.pl [num-cpus] [diviner-based-inputs/] [genome.fa] [out-dir-name]\n\n"; }



my $num_cpus = int($ARGV[0]);


my $div_inputs_dir_name = $ARGV[1];
die "\n  ERROR:  Failed to locate Diviner-based input directory '$div_inputs_dir_name'\n\n"
	if (!(-d $div_inputs_dir_name));
$div_inputs_dir_name = $div_inputs_dir_name.'/' if ($div_inputs_dir_name !~ /\/$/);


my $genome = $ARGV[2];
die "\n  ERROR:  Failed to locate genome file '$genome'\n\n"
	if (!(-e $genome));


my $out_dir_name = $ARGV[3];
die "\n  ERROR:  Desired output directory '$out_dir_name' already exists\n\n"
	if (-d $out_dir_name);
$out_dir_name = $out_dir_name.'/' if ($out_dir_name !~ /\/$/);

die "\n  ERROR:  Failed to create output directory '$out_dir_name'\n\n"
	if (system("mkdir \"$out_dir_name\""));


opendir(my $DivDir,$div_inputs_dir_name)
	|| die "\n  ERROR:  Failed to open input directory '$div_inputs_dir_name'\n\n";

my @Cmds;
while (my $gene = readdir($DivDir))
{
	$gene =~ s/\/$//;
	next if ($gene =~ /^\./);

	my $gene_ali_file_name = $div_inputs_dir_name.$gene.'/'.$gene.'.afa';
	next if (!(-e $gene_ali_file_name));

	my $out_gene_dir_name = $out_dir_name.$gene.'/';
	die "\n  ERROR:  Failed to create gene output directory '$out_gene_dir_name'\n\n"
		if (system("mkdir \"$out_gene_dir_name\""));

	my $out_file_name = $out_gene_dir_name.$gene.'.out';
	my $err_file_name = $out_gene_dir_name.$gene.'.err';

	my $cmd = '/usr/bin/time -v '.$BATHSEARCH." --qformat afa -o \"$out_file_name\" \"$gene_ali_file_name\" \"$genome\" 2>\"$err_file_name\"";
	push(@Cmds,$cmd);

}

closedir($DivDir);



if (scalar(@Cmds) == 0)
{
	print "\n  No query alignment files found in '$div_inputs_dir_name'\n\n";
	system("rmdir \"$out_dir_name\"");
	exit(0);
}



my $num_cmds = scalar(@Cmds);
$num_cpus = $num_cmds if ($num_cpus > $num_cmds);

my $thread_id = 0;
my $active_threads = 1;
my $pid = 0;
while ($active_threads < $num_cpus)
{
	if ($pid = fork)
	{
		die "\n  ERROR:  Fork failed\n\n" if (not defined $pid);
		$active_threads++;
	}
	else
	{
		$thread_id = $active_threads;
		last;
	}
}


my $cmd_portion  = int($num_cmds/$num_cpus);
my $start_cmd_id =  $thread_id    * $cmd_portion;
my $end_cmd_id   = ($thread_id+1) * $cmd_portion;
$end_cmd_id = $num_cmds if ($thread_id == $num_cpus-1);


for (my $cmd_id = $start_cmd_id; $cmd_id < $end_cmd_id; $cmd_id++)
{
	if (system($Cmds[$cmd_id]))
	{
		print "  WARNING:  Command '$Cmds[$cmd_id]' failed!\n";
	}
}


exit(0) if ($thread_id);
while (wait() != -1) {}


1;
