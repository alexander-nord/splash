#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;



my $BATHBUILD = $0;
$BATHBUILD =~ s/[^\/]+$//;
$BATHBUILD = './' if (!$BATHBUILD);
$BATHBUILD = $BATHBUILD.'../../BATH/src/bathbuild';
die "\n  ERROR:  Failed to locate bathbuild (looking for '$BATHBUILD')\n\n"
	if (!(-e $BATHBUILD));



if (@ARGV != 2) { die "\n  USAGE:  ./ExtractDivinerQueries.pl [Diviner-Results/] [out-dir-name/]\n\n"; }



my $diviner_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate Diviner results directory '$diviner_dir_name'\n\n"
	if (!(-d $diviner_dir_name));
$diviner_dir_name = $diviner_dir_name.'/' if ($diviner_dir_name !~ /\/$/);


my $out_dir_name = $ARGV[1];
die "\n  ERROR:  Output directory '$out_dir_name' already exists\n\n"
	if (-d $out_dir_name);
$out_dir_name = $out_dir_name.'/';


my $hbp_file_name = $diviner_dir_name.'Hits-by-Pct-ID.out';
open(my $HBP,'<',$hbp_file_name)
	|| die "\n  ERROR:  Failed to open Hits-by-Pct-ID file '$hbp_file_name'\n\n";

my %GenesToSourceSpecies;
my %GenesToTargetRanges;
while (my $line = <$HBP>)
{
	if ($line =~ /:\[ Novel \] (\S+)\s+(\S+)\s+\-\>\s+homo_sapiens\s+\(([^\)]+)\)/)
	{

		my $gene = $1;
		my $source_species = $2;
		my $target_range = $3;

		if ($GenesToSourceSpecies{$gene})
		{ 
			$GenesToSourceSpecies{$gene} = $GenesToSourceSpecies{$gene}.'/'.$source_species; 
			$GenesToTargetRanges{ $gene} = $GenesToTargetRanges{ $gene}.'/'.$target_range;
		}
		else
		{
			$GenesToSourceSpecies{$gene} = $source_species;
			$GenesToTargetRanges{ $gene} = $target_range;
		}
	
	}
}

close($HBP);


if (scalar(keys %GenesToSourceSpecies) == 0)
{
	print "\n  No novel hits to the human genome observed in Hits-by-Pct-ID file '$hbp_file_name'\n\n";
	exit(0);
}


# Looks like we're going to be making some dank BATH HMMs!
die "\n  ERROR:  Failed to create output directory '$out_dir_name'\n\n"
	if (system("mkdir \"$out_dir_name\""));


foreach my $gene (sort keys %GenesToSourceSpecies)
{

	my %ValidSources;
	foreach my $species (split(/\//,$GenesToSourceSpecies{$gene}))
	{
		$ValidSources{lc($species)} = 1;
	}


	my $gene_ali_file_name = $diviner_dir_name.'Results-by-Gene/'.$gene.'/mapped-seqs.afa';
	open(my $GeneAliFile,'<',$gene_ali_file_name)
		|| die "\n  ERROR:  Failed to open input alignment file '$gene_ali_file_name'\n\n";

	my @InitMSA;
	my @SeqNames;
	my $num_seqs = -1;
	my $init_msa_len;
	my $valid_source = 0;
	while (my $line = <$GeneAliFile>)
	{
		
		$line =~ s/\n|\r//g;
		next if (!$line);

		if ($line =~ /\>([^\|]+)\|/)
		{
			my $next_source_species = lc($1);
			if ($ValidSources{$next_source_species})
			{

				$valid_source = 1;
				$num_seqs++;
				$init_msa_len=0;
			
				$line =~ /\/>(\S+)/;
				$SeqNames[$num_seqs] = $1;

			}
			else
			{
				$valid_source = 0;
			}

		}
		elsif ($valid_source)
		{
			foreach my $char (split(//,$line))
			{
				$InitMSA[$num_seqs][$init_msa_len] = $char;
				$init_msa_len++;
			}
		}

	}
	close($GeneAliFile);
	$num_seqs++;


	# Reduce our initial MSA to columns with amino acid content
	my @MSA;
	my $msa_len = 0;
	for (my $init_col = 0; $init_col < $init_msa_len; $init_col++)
	{
		
		my $has_amino_content = 0;
		for (my $seq_id=0; $seq_id<$num_seqs; $seq_id++)
		{
			if ($InitMSA[$seq_id][$init_col] =~ /[A-Za-z]/)
			{
				$has_amino_content = 1;
			}
			$MSA[$seq_id][$msa_len] = $InitMSA[$seq_id][$init_col];
		}

		$msa_len++ if ($has_amino_content);

	}


	my $gene_out_dir_name = $out_dir_name.$gene.'/';
	die "\n  ERROR:  Failed to create gene output directory '$gene_out_dir_name'\n\n"
		if (system("mkdir \"$gene_out_dir_name\""));


	my $ali_out_file_name = $gene_out_dir_name.$gene.'.afa';
	open($GeneAliFile,'>',$ali_out_file_name)
		|| die "\n  ERROR:  Failed to open alignment output file '$ali_out_file_name'\n\n";

	for (my $seq_id=0; $seq_id<$num_seqs; $seq_id++)
	{
		print $GeneAliFile ">$SeqNames[$seq_id]";
		for (my $col=0; $col<$msa_len; $col++)
		{
			print $GeneAliFile "\n" if ($col % 60 == 0);
			print $GeneAliFile "$MSA[$seq_id][$col]";
		}
		print $GeneAliFile "\n";
	}

	close($GeneAliFile);


	my $hmm_out_file_name = $gene_out_dir_name.$gene.'.bath.hmm';
	my $bathbuild_cmd = " --amino \"$hmm_out_file_name\" \"$ali_out_file_name\"";
	$bathbuild_cmd = " --unali".$bathbuild_cmd if ($num_seqs == 1);
	$bathbuild_cmd = $BATHBUILD.$bathbuild_cmd;

	if (system($bathbuild_cmd))
	{
		die "\n  ERROR:  bathbuild command '$bathbuild_cmd' failed\n\n";
	}


	# As one last thing, we'll want to record the target range we're
	# hoping to "discover" with splashBATH
	my @TargetRanges = split(/\//,$GenesToTargetRanges{$gene});
	
	$TargetRanges[0] =~ /^(\S+):(\d+)\.\.(\d+)$/;
	my $target_chr   = $1;
	my $target_start = $2;
	my $target_end   = $3;

	my $target_revcomp = 0;
	$target_revcomp = 1 if ($target_chr =~ /\[revcomp\]/);

	for (my $range_id=1; $range_id<scalar(@TargetRanges); $range_id++)
	{
		$TargetRanges[$range_id] =~ /:(\d+)\.\.(\d+)$/;
		my $range_start = $1;
		my $range_end   = $2;

		if ($target_revcomp)
		{
			$target_start = $range_start if ($range_start > $target_start);
			$target_end   = $range_end   if ($range_end   < $target_end);
		}
		else
		{
			$target_start = $range_start if ($range_start < $target_start);
			$target_end   = $range_end   if ($range_end   > $target_end);			
		}
	}

	my $target_out_file_name = $gene_out_dir_name.'new-exon-range.out';
	system("echo \"CHR  : $target_chr\nRANGE: $target_start..$target_end\" > $target_out_file_name");


}



1;
