#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub RecordMappings;
sub RecordSeqs;
sub WriteSeqToFile;
sub CleanUpEmptyGeneDirs;



if (@ARGV != 1) {
	die "\n  USAGE:  ./MirageOutToSplashIn.pl [Mirage-Results]\n\n";
}


my $in_dir_name = $ARGV[0];
$in_dir_name = $in_dir_name.'/' if ($in_dir_name !~ /\/$/);
if (!(-d $in_dir_name)) {
	die "\n  ERROR:  Failed to locate input directory '$in_dir_name'\n\n";
}


my $out_dir_name = 'mirage-based-data/';
if (-d $out_dir_name) {
	die "\n  ERROR:  Output directory '$out_dir_name' already exists\n\n";
} elsif (system("mkdir \"$out_dir_name\"")) {
	die "\n  ERROR:  Failed to create directory '$out_dir_name'\n\n";
}


my $all_species_dir_name = $in_dir_name.'Species-MSAs/';
if (!(-d $all_species_dir_name)) {
	die "\n  ERROR:  Failed to locate species-specific head directory '$all_species_dir_name'\n\n";
}


opendir(my $AllSpeciesDir,$all_species_dir_name) || die "\n  ERROR:  Failed to open input directory '$all_species_dir_name'\n\n";
while (my $species = readdir($AllSpeciesDir)) {

	next if ($species =~ /^\./);
	
	$species =~ s/\/$//;
	my $species_dir_name     = $all_species_dir_name.$species.'/';
	my $species_ali_dir_name = $species_dir_name.'alignments/';
	my $species_map_dir_name = $species_dir_name.'mappings/';

	next if (!(-d $species_dir_name && -d $species_ali_dir_name && -d $species_map_dir_name));

	opendir(my $SpeciesMapDir,$species_map_dir_name);
	while (my $map_file_name = readdir($SpeciesMapDir)) {

		next if ($map_file_name !~ /^(\S+)\.out$/);
		my $gene = $1;

		# Titan is really slow, so for now we're just going to cut it out
		next if ($gene eq 'ttn');

		$map_file_name = $species_map_dir_name.$map_file_name;

		my $mapped_seqs_ref = RecordMappings($map_file_name,$out_dir_name);

		my $ali_file_name = $species_ali_dir_name.$gene.'.afa';

		RecordSeqs($mapped_seqs_ref,$ali_file_name,$out_dir_name);

	}
	closedir($SpeciesMapDir);

}
closedir($AllSpeciesDir);



# If there were genes where all mappings were derived from
# non-FastMap2 methods, we'll want to be sure to clear out those
# directories
CleanUpEmptyGeneDirs($out_dir_name);



1;







############################################################################
#
#  Subroutine: RecordMappings
#
sub RecordMappings
{

	my $map_file_name = shift;
	my $out_dir_name  = shift;


	$map_file_name =~ /\/([^\/]+)\/mappings\/([^\/]+)\.out$/;
	my $species = $1;
	my $gene    = $2;


	my $gene_out_dir_name = $out_dir_name.$gene.'/';
	if (!(-d $gene_out_dir_name) && system("mkdir \"$gene_out_dir_name\"")) {
		die "\n  ERROR:  Failed to create gene output directory '$gene_out_dir_name'\n\n"
	}


	open(my $MapFile,'<',$map_file_name) || die "\n  ERROR:  Failed to open input mapping file '$map_file_name'\n\n";


	my $canon_chr = <$MapFile>;
	$canon_chr =~ /Canonical Chromosome: (\S+)/;
	$canon_chr = $1;
	
	my $straight_chr = $canon_chr;
	$straight_chr =~ s/\[revcomp\]//;

	my %MappedSeqs;

	my $min_coord = -1;
	my $max_coord = -1;

	while (my $line = <$MapFile>) {

		next if ($line !~ /Sequence ID: (\S+)/);
		my $full_seq_name = $1;

		my $map_method = <$MapFile>;
		my $seq_chr    = <$MapFile>;
		my $num_exons  = <$MapFile>;

		next if ($seq_chr !~ /Chromosome : (\S+)/);
		$seq_chr = $1;

		next if ($seq_chr ne $canon_chr);

		$map_method =~ /Map Method : (\S+)/;
		$map_method = $1;

		# For now I'm going to require that the mapping came
		# from FastMap2, since I want to test with mapping
		# coordinates that I trust *very much*
		#
		# UPDATE: We'll catch these on the back-end
		#
		#next if ($map_method ne 'FastMap2');
		#

		$num_exons =~ /Num Exons  : (\d+)/;
		$num_exons = $1;


		# If there's only one exon, this isn't something we
		# want to search with Splash (duh)
		next if ($num_exons == 1);


		$full_seq_name =~ /\|([^\|]+)$/;
		my $seq_id = $1;

		my $seq_map_out_file_name = $gene_out_dir_name.$species.'.'.$seq_id.'.map.out';
		open(my $SeqMapOutFile,'>',$seq_map_out_file_name) || die "\n  ERROR:  Failed to open map output file '$seq_map_out_file_name'\n\n";


		print $SeqMapOutFile "SEQ   : $full_seq_name\n";
		print $SeqMapOutFile "CHR   : $seq_chr\n";
		print $SeqMapOutFile "EXONS : $num_exons\n";


		for (my $exon_id=0; $exon_id<$num_exons; $exon_id++) {

			my $map_line = <$MapFile>;
			<$MapFile>; # Eat the actual codon mapping

			$map_line =~ /:(\d+)\.\.(\d+)/;
			my $exon_start_nucl = $1;
			my $exon_end_nucl   = $2;

			print $SeqMapOutFile "$exon_start_nucl..$exon_end_nucl\n";

			if ($min_coord == -1 || $exon_start_nucl < $min_coord) {
				$min_coord = $exon_start_nucl;
			}
			if ($max_coord == -1 || $exon_start_nucl > $max_coord) {
				$max_coord = $exon_start_nucl;
			}

			if ($exon_end_nucl < $min_coord) {
				$min_coord = $exon_end_nucl;
			}
			if ($exon_end_nucl > $max_coord) {
				$max_coord = $exon_end_nucl;
			}

		}

		close($SeqMapOutFile);

		$MappedSeqs{$full_seq_name} = 1;

	}
	
	close($MapFile);


	# If we had mappings but none were FastMap2-derived
	# then we don't have anything to spit out here
	if ($min_coord != $max_coord) # (Both == -1)
	{
		my $species_range_file_name = $gene_out_dir_name.$species.'.genome-range.out';
		open(my $SpeciesRangeFile,'>',$species_range_file_name) || die "\n  ERROR:  Failed to open range output file '$species_range_file_name'\n\n";
		print $SpeciesRangeFile "SPECIES : $species\n";
		print $SpeciesRangeFile "CHR     : $straight_chr\n";
		print $SpeciesRangeFile "RANGE   : $min_coord..$max_coord\n";
		close($SpeciesRangeFile);
	}

	return \%MappedSeqs;

}







############################################################################
#
#  Subroutine: RecordSeqs
#
sub RecordSeqs
{
	my $mapped_seqs_ref = shift;
	my $ali_file_name   = shift;
	my $out_dir_name    = shift;

	my %MappedSeqs = %{$mapped_seqs_ref};

	# Sanity check...
	return if (scalar(keys %MappedSeqs) == 0);

	$ali_file_name =~ /\/([^\/]+)\/alignments\/([^\/]+)\.afa$/;
	my $species = $1;
	my $gene    = $2;

	my $gene_out_dir_name = $out_dir_name.$gene.'/';


	open(my $AliFile,'<',$ali_file_name) || die "\n  ERROR:  Failed to open input alignment file '$ali_file_name'\n\n";
	
	my $in_mapped_seq = 0;
	my $next_seq_name = '';
	my $next_seq      = '';
	while (my $line = <$AliFile>) {

		$line =~ s/\n|\r//g;
		next if (!$line);

		if ($line =~ /^\>/) {

			if ($in_mapped_seq) {
				$next_seq_name =~ /\|([^\|]+)$/;
				my $seq_id = $1;
				my $seq_file_name = $gene_out_dir_name.$species.'.'.$seq_id.'.fa';
				WriteSeqToFile($next_seq_name,$next_seq,$seq_file_name);
			}

			$line =~ /\>(\S+)/;
			$next_seq_name = $1;
			$next_seq      = '';

			if  ($MappedSeqs{$next_seq_name}) { $in_mapped_seq = 1; }
			else                              { $in_mapped_seq = 0; }

		} elsif ($in_mapped_seq) {

			$line =~ s/\s//g;
			$line =~ s/\-//g;
			$line =~ s/\///g;
			if ($line =~ /\S/) {
				$next_seq = $next_seq.uc($line);
			}

		}

	}
	close($AliFile);


	if ($in_mapped_seq) {
		$next_seq_name =~ /\|([^\|]+)$/;
		my $seq_id = $1;
		my $seq_file_name = $gene_out_dir_name.$species.'.'.$seq_id.'.fa';
		WriteSeqToFile($next_seq_name,$next_seq,$seq_file_name);
	}


}







############################################################################
#
#  Subroutine: WriteSeqToFile
#
sub WriteSeqToFile
{

	my $seq_name  = shift;
	my $seq_str   = shift;
	my $file_name = shift;

	my @Seq = split(//,$seq_str);

	open(my $SeqFile,'>',$file_name) || die "\n  ERROR:  Failed to open sequence output file '$file_name'\n\n";

	print $SeqFile ">$seq_name";
	my $line_pos =  0;
	my $line_len = 60;
	for (my $seq_pos = 0; $seq_pos < scalar(@Seq); $seq_pos++) {
		print $SeqFile "\n" if ($seq_pos % 60 == 0);
		print $SeqFile "$Seq[$seq_pos]";
	}
	print $SeqFile "\n";

	close($SeqFile);

}







############################################################################
#
#  Subroutine: CleanUpEmptyGeneDirs
#
sub CleanUpEmptyGeneDirs
{
	my $out_dir_name = shift;

	opendir(my $OutDir,$out_dir_name)
		|| die "\n  ERROR:  Failed to open output directory $out_dir_name for cleanup check\n\n";

	while (my $gene_dir_name = readdir($OutDir)) {

		next if ($gene_dir_name =~ /^\./);
		$gene_dir_name = $gene_dir_name.'/' if ($gene_dir_name !~ /\/$/);
		$gene_dir_name = $out_dir_name.$gene_dir_name;

		# I think everything in the output directory should be
		# a directory, but just in case...
		next if (!(-d $gene_dir_name));

		opendir(my $GeneDir,$gene_dir_name)
			|| die "\n  ERROR:  Failed to open gene directory $gene_dir_name for cleanup check\n\n";

		my $dir_is_empty = 1;
		while (my $file = readdir($GeneDir)) {
			next if ($file !~ /\.fa$/);
			$dir_is_empty = 0;
			last;
		}

		closedir($GeneDir);

		system("rm -rf \"$gene_dir_name\"") if ($dir_is_empty);

	}

	closedir($OutDir);

}
















