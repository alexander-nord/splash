#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

if (@ARGV != 4) {
    die "\n  USAGE:  ./GenHMMSearchData.pl [Mirage-Results] [Species-Guide] [HSI] [Out-Dir-Name]\n\n";
}


my $out_fa_line_len = 60;


my $hsi_dir_name = $ARGV[2];
$hsi_dir_name =~ s/\/$//;
my $sfetch = $hsi_dir_name.'/build/sfetch';
my $sstat  = $hsi_dir_name.'/build/sstat';

if (!(-e $sfetch && -e $sstat)) {
    die "\n  ERROR:  Failed to locate HSI tools (looking for '$sfetch' and '$sstat')\n\n";
}


open(my $SpeciesGuide,'<',$ARGV[1]) || die "\n  ERROR:  Failed to open species guide file '$ARGV[1]'\n\n";
my %SpeciesToGenome;
my %SpeciesChrToLen;
while (my $line = <$SpeciesGuide>) {

    next if ($line !~ /^(\S+)\s+(\S+)\s+\S+\s*$/);
    my $species = lc($1);
    my $genome  = $2;

    if (!(-e $genome)) {
	print "\n  WARNING:  $species genome location ($genome) appears invalid\n\n";
	next;
    }

    $SpeciesToGenome{$species} = $genome;
    
    open(my $SpeciesSstat,"$sstat $genome |") || die "\n  ERROR:  Failed to run sstat on '$genome'\n\n";
    while (my $line = <$SpeciesSstat>) {
	if ($line =~ /^\:\s+(\S+)\s+(\d+)/) {
	    $SpeciesChrToLen{$species.'|'.$1} = $2;
	}
    }
    close($SpeciesSstat);
   
}
close($SpeciesGuide);


my $mirage_results_dir_name = $ARGV[0];
$mirage_results_dir_name =~ s/\/$//;
$mirage_results_dir_name = $mirage_results_dir_name.'/Species-MSAs/';

if (!(-d $mirage_results_dir_name)) {
    die "\n  ERROR:  Failed to locate Mirage species-specific results directory '$mirage_results_dir_name'\n\n";
}


my $out_dir_name = $ARGV[3];
$out_dir_name =~ s/\/$//;
$out_dir_name = $out_dir_name.'/';
if (-d $out_dir_name) {
    system("rm -rf $out_dir_name");
}
if (system("mkdir $out_dir_name")) {
    die "\n  ERROR:  Failed to create output directory '$out_dir_name'\n\n";
}


opendir(my $MirageResultsDir, $mirage_results_dir_name) || die "\n  ERROR:  Failed to open directory '$mirage_results_dir_name'\n\n";
while (my $species = readdir($MirageResultsDir)) {

    next if ($species =~ /^\./);

    $species =~ s/\/$//;

    my $mapping_dir_name   = $mirage_results_dir_name.$species.'/mappings/';
    my $alignment_dir_name = $mirage_results_dir_name.$species.'/alignments/';

    next if (!(-d $mapping_dir_name && -d $alignment_dir_name));

    my $species_out_dir_name = $out_dir_name.$species.'/';
    if (system("mkdir $species_out_dir_name")) {
	die "\n  ERROR:  Failed to create species output directory '$species_out_dir_name'\n\n";
    }

    my $genome = $SpeciesToGenome{$species};

    opendir(my $MappingDir, $mapping_dir_name) || die "\n  ERROR:  Failed to open directory '$mapping_dir_name'\n\n";

    while (my $gene_map_file_name = readdir($MappingDir)) {

	next if ($gene_map_file_name !~ /^(\S+)\.out$/);
	my $gene = $1;

	$gene_map_file_name = $mapping_dir_name.$gene_map_file_name;
	my $gene_ali_file_name = $alignment_dir_name.$gene.'.afa';

	next if (!(-e $gene_ali_file_name));

	
	my @Seqs;
	my @SeqNames;
	my $num_seqs = -1;
	open(my $AliFile,'<',$gene_ali_file_name);
	while (my $line = <$AliFile>) {

	    $line =~ s/\n|\r//g;
	    next if (!$line);

	    if ($line =~ /\>(\S+)/) {
		my $seq_name = $1;
		$num_seqs++;
		$SeqNames[$num_seqs] = $seq_name;
		$Seqs[$num_seqs] = '';
	    } else {
		#$line =~ s/\-//g; # Input to hmmbuild is an alignment (DUH!!!)
		$line =~ s/\///g;
		$line =~ s/\s//g;
		next if (!$line);
		$Seqs[$num_seqs] = $Seqs[$num_seqs].lc($line);
	    }

	}
	close($AliFile);
	$num_seqs++;


	my $gene_out_dir_name = $species_out_dir_name.$gene.'/';
	if (system("mkdir $gene_out_dir_name")) {
	    die "\n  ERROR:  Failed to create gene output directory '$gene_out_dir_name'\n\n";
	}
	

	my $all_prots_fname = $gene_out_dir_name.'isos.afa';
	open(my $AllProts,'>',$all_prots_fname) || die "\n  ERROR:  Failed to open output file '$all_prots_fname'\n\n";

	for (my $i=0; $i<$num_seqs; $i++) {

	    $SeqNames[$i] =~ /\|([^\|]+)$/;
	    my $iso_id = $1;

	    #my $iso_out_file_name = $gene_out_dir_name.$iso_id.'.fa';
	    #open(my $IsoOutFile,'>',$iso_out_file_name) || die "\n  ERROR:  Failed to open output file '$iso_out_file_name'\n\n";
	    #print $IsoOutFile ">$gene\_$i\n";

	    print $AllProts ">$gene\_$i\n";

	    my @SeqChars = split(//,$Seqs[$i]);
	    my $char_id = 0;
	    while ($char_id < scalar(@SeqChars)) {

		print $AllProts "$SeqChars[$char_id++]";
		print $AllProts "\n" if ($char_id % $out_fa_line_len == 0);

		#print $IsoOutFile "$SeqChars[$char_id++]";
		#print $IsoOutFile "\n" if ($char_id % $out_fa_line_len == 0);
	    }

	    print $AllProts "\n" if ($char_id % $out_fa_line_len != 0);

	    #print $IsoOutFile "\n" if ($char_id % $out_fa_line_len != 0);
	    #print $IsoOutFile "\n";
	    #close($IsoOutFile);

	}
	close($AllProts);

	
	open(my $GeneMapFile,'<',$gene_map_file_name) || die "\n  ERROR:  Failed to open input file '$gene_map_file_name'\n\n";

	my $canon_chr = <$GeneMapFile>;
	$canon_chr =~ s/^Canonical Chromosome\: |\n|\r|\s//g;

	my $revcomp = 0;
	if ($canon_chr =~ /^(\S+)\[revcomp\]/) {
	    $canon_chr = $1;
	    $revcomp = 1;
	}

	# In case the 'canon_chr' is Undetermined we don't want
	# to worry about this gene
	if (!$SpeciesChrToLen{$species.'|'.$canon_chr}) {
	    system("rm -rf $gene_out_dir_name");
	    next;
	}

	my $min_coord = -1;
	my $max_coord = -1;

	while (my $line = <$GeneMapFile>) {

	    next if ($line !~ /Chromosome \: $canon_chr/);

	    $line = <$GeneMapFile>; # Num Exons

	    $line = <$GeneMapFile>;
	    while (!eof($GeneMapFile) && $line =~ /\* Aminos \S+ \S+\:(\d+)\.\.(\d+)/) {

		my $exon_start = $1;
		my $exon_end = $2;

		if ($revcomp) {
		    if ($min_coord == -1) {
			$min_coord = $exon_end;
			$max_coord = $exon_start;
		    } else {
			if ($min_coord > $exon_end) {
			    $min_coord = $exon_end;
			}
			if ($max_coord < $exon_start) {
			    $max_coord = $exon_start;
			}
		    }
		} else {
		    if ($min_coord == -1) {
			$min_coord = $exon_start;
			$max_coord = $exon_end;
		    } else {
			if ($min_coord > $exon_start) {
			    $min_coord = $exon_start;
			}
			if ($max_coord < $exon_end) {
			    $max_coord = $exon_end;
			}
		    }
		}

		$line = <$GeneMapFile>;
		$line = <$GeneMapFile>;

	    }

	}

	close($GeneMapFile);

	if ($min_coord  - 4000 > 0) {
	    $min_coord -= 4000;
	} else {
	    $min_coord = 1;
	}

	if ($max_coord  + 4000 < $SpeciesChrToLen{$species.'|'.$canon_chr}) {
	    $max_coord += 4000;
	} else {
	    $max_coord = $SpeciesChrToLen{$species.'|'.$canon_chr} - 1;
	}

	my $range = $min_coord.'..'.$max_coord;
	if ($revcomp) {
	    $range = $max_coord.'..'.$min_coord;
	}

	my $sfetch_cmd = "$sfetch -range $range $genome $canon_chr";
	open(my $NuclSfetch,"$sfetch_cmd |") || die "\n  ERROR:  Failed to run sfetch command '$sfetch_cmd'\n\n";

	my $gene_nucls_file_name = $gene_out_dir_name.'dna.fa';
	open (my $GeneNuclsFile,'>',$gene_nucls_file_name) || die "\n  ERROR:  Failed to open output file '$gene_nucls_file_name'\n\n";


	<$NuclSfetch>; # Eat the header
	print $GeneNuclsFile ">$gene\_coding\_region\n";
	while (my $line = <$NuclSfetch>) {
	    $line =~ s/\n|\r//g;
	    print $GeneNuclsFile "$line\n" if ($line);
	}

	close($NuclSfetch);
	close($GeneNuclsFile);


	open(my $SfetchRecord,'>',$gene_out_dir_name.'sfetch_cmd.txt');
	print $SfetchRecord "$sfetch_cmd\n";
	close($SfetchRecord);

    }
    closedir($MappingDir);

}
closedir($MirageResultsDir);


1;
