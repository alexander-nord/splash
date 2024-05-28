#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


my $pbhmms_dir_name = $0;
$pbhmms_dir_name = $pbhmms_dir_name.'/' if ($pbhmms_dir_name !~ /\/$/);

my $BATHBUILD = $pbhmms_dir_name.'../../BATH/src/bathbuild';
die "\n  ERROR:  Failed to locate 'bathbuild' (looking for '$BATHBUILD')\n\n"
	if (!(-e $BATHBUILD));

my $BATHCONVERT = $pbhmms_dir_name.'../../BATH/src/bathconvert';
die "\n  ERROR:  Failed to locate 'bathconvert' (looking for '$BATHCONVERT')\n\n"
	if (!(-e $BATHCONVERT));


if (@ARGV != 1) { die "\n  USAGE:  ./PrepBathHMMs.pl [/path/to/gene/super-dir/]\n\n"; }


my $all_genes_dir_name = $ARGV[0];
$all_genes_dir_name = $all_genes_dir_name.'/' if ($all_genes_dir_name !~ /\/$/);
die "\n  ERROR:  Failed to locate input directory '$all_genes_dir_name'\n\n"
	if (!(-d $all_genes_dir_name));


opendir(my $AllGenesDir,$all_genes_dir_name)
	|| die "\n  ERROR:  Failed to open gene super-directory '$all_genes_dir_name'\n\n";

while (my $gene = readdir($AllGenesDir))
{
	
	next if ($gene =~ /^\./);
	$gene =~ s/\/$//;

	my $gene_dir_name = $all_genes_dir_name.$gene.'/';
	next if (!(-d $gene_dir_name));

	opendir(my $GeneDir,$gene_dir_name)
		|| die "\n  ERROR:  Failed to open gene directory '$gene_dir_name'\n\n";

	while (my $file_name = readdir($GeneDir))
	{
		
		next if ($file_name =~ /\.bath\.hmm/);
		next if (lc($file_name) !~ /\.a?fa[sta]?$|\.hmm$/)

		$file_name =~ /^(\S+)\.([^\.]+)$/;
		my $file_base_name = $1;
		my $file_extension = $2;

		next if (-e $gene_dir_name.$file_base_name.'.bath.hmm');

		$file_name        = $gene_dir_name.$file_name;
		my $out_file_name = $gene_dir_name.$file_base_name.'.bath.hmm';

		if (lc($file_extension) eq 'hmm')
		{

			my $convert_cmd = "$BATHCONVERT \"$out_file_name\" \"$file_name\" 1>/dev/null";

			if (system($convert_cmd))
			{
				die "\n  ERROR:  HMM conversion command '$convert_cmd' failed\n\n";
			}

		}
		else
		{

			my $build_cmd = " --amino \"$out_file_name\" \"$file_name\" 1>/dev/null";
			$build_cmd = ' --unali'.$build_cmd
				if (lc($file_extension) !~ /^a/);

			$build_cmd = $BATHBUILD.$build_cmd;

			if (system($build_cmd))
			{
				die "\n  ERROR:  HMM build command '$build_cmd' failed\n\n";
			}

		}

	}

	closedir($GeneDir);

}

closedir($AllGenesDir);


1;
