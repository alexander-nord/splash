#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


#  FUNCTIONAL SUBROUTINES
#
sub PrintUsage;
sub PullConsensusSeq;
sub EslSfetch;
sub VisAli;
sub TranslateCodon;
sub FormatForOutput;


#  BEUREAUCRATIC SUBROUTINES
#
sub MIN;
sub MAX;
sub RunSystemCommand;
sub OpenSystemCommand;
sub OpenInputFile;
sub NameOutputFile;
sub OpenOutputFile;
sub OpenFileToAppendTo;
sub DeleteFile;
sub ConfirmDirectory;
sub OpenDirectory;
sub NameDirectory;
sub CreateDirectory;



############
#          #
#  SCRIPT  #
#          #
############


if (@ARGV != 1) { PrintUsage(); }

my $linewidth = 66;

my @ProtFiles;
my @DNAFiles;

my $inf = OpenInputFile($ARGV[0]);
while (my $line = <$inf>) {

    $line =~ s/\n|\r//g;
    next if (!$line);

    $line =~ /^\s*(\S+)\s+(\S+)\s*$/;
    my $protname = $1;
    my $dnaname  = $2;

    if (!(-e $protname)) { die "\n  Failed to locate HMM file '$protname'\n\n"; }
    if (!(-e $dnaname )) { die "\n  Failed to locate DNA file '$dnaname'\n\n";  }

    if (!(-e $dnaname.'.ssi')) { RunSystemCommand("esl-sfetch --index \"$dnaname\""); }

    push(@ProtFiles,$protname);
    push(@DNAFiles,$dnaname);
    
}
close($inf);

if (!scalar(@ProtFiles)) { die "\n  No pairs found in input file '$ARGV[0]'?\n\n"; }

# Actually iterate over the hmm/sequence pairs, running hmmsearcht and parsing the
# outputs.
for (my $filenum = 0; $filenum < scalar(@ProtFiles); $filenum++) {

    my $protfile = $ProtFiles[$filenum];
    my $dnafile  = $DNAFiles[$filenum];

    my $consensus = PullConsensusSeq($protfile);
    my @ConsChars = split(//,$consensus);

    my $searchcmd = "hmmsearcht \"$protfile\" \"$dnafile\"";
    my $searchout = OpenSystemCommand($searchcmd);    

    my $num_hits = 0;
    while (my $line = <$searchout>) {
	$line =~ s/\n|\r//g;
	next if ($line !~ /RESULTS TIME\! \((\d+)\)/);
	$num_hits = $1;
	last;
    }

    my %SplicedHits;
    my %OrigHitDoms;

    for (my $hitnum = 0; $hitnum < $num_hits; $hitnum++) {

	# First up, the line with all of the specifics for this spliced hit
	my $line = <$searchout>;

	$line =~ /\((\d+)\:(\d+)\,(\d+)\:(\d+)\)/;
	my $hit1id = $1;
	my $dom1id = $2;
	my $hit2id = $3;
	my $dom2id = $4;

	$line =~ /(\d+)\s*\-\>\s*(\d+)\s*\/\s*(\d+)\s*\-\>\s*(\d+)/;
	my $hmmstart1 = $1;
	my $hmmend1   = $2;
	my $hmmstart2 = $3;
	my $hmmend2   = $4;

	$line =~ /\[(\S+)\: (\d+)\.\.(\d+) \/ (\d+)\.\.(\d+)\]/;
	my $chrname   = $1;
	my $dnastart1 = $2;
	my $dnaend1   = $3;
	my $dnastart2 = $4;
	my $dnaend2   = $5;

	$line =~ /\=\>\s+(\S+)/;
	my $score = $1;

	my $key   = "$hit1id\:$dom1id\|$hit2id\:$dom2id";
	my $entry = "$hmmstart1|$hmmend1|$dnastart1|$dnaend1";
	$entry    = $entry."|$hmmstart2|$hmmend2|$dnastart2|$dnaend2";
	$entry    = $entry."|$chrname|$score";
	$SplicedHits{$key} = $entry;

	# Second, some identifying info. so we can pair each side of the spliced
	# alignment with an original translated hit.
	$line = <$searchout>;

	$line =~ /\[HMMRange\:(\d+)\-(\d+)\|DNARange\:(\d+)\-(\d+)\] \[/;
	my $orig_hmmstart1 = $1;
	my $orig_hmmend1   = $2;
	my $orig_dnastart1 = $3;
	my $orig_dnaend1   = $4;

	$key   = "$orig_hmmstart1:$orig_hmmend1:$orig_dnastart1:$orig_dnaend1";
	$entry = "$hit1id:$dom1id";
	$OrigHitDoms{$key} = $entry;

	$line =~ /\] \[HMMRange\:(\d+)\-(\d+)\|DNARange\:(\d+)\-(\d+)\]/;
	my $orig_hmmstart2 = $1;
	my $orig_hmmend2   = $2;
	my $orig_dnastart2 = $3;
	my $orig_dnaend2   = $4;

	$key   = "$orig_hmmstart2:$orig_hmmend2:$orig_dnastart2:$orig_dnaend2";
	$entry = "$hit2id:$dom2id";
	$OrigHitDoms{$key} = $entry;

    }
    
    # Once we're done with the spliced hits, we'll need to take a nice little
    # peek at the original hits.
    my %IDsToOutput;
    while (my $line = <$searchout>) {

	$line =~ s/\n|\r//g;
	next if ($line !~ /\>\>/);

	$line = <$searchout>; # header
	$line = <$searchout>; # ------

	# Catch the locations of all the domains for this hit
	my @HitCoords;
	$line = <$searchout>;
	while ($line =~ /\S+/) {

	    my @LineData = split(/\s+/,$line);

	    my $hmmfrom = $LineData[7];
	    my $hmmto   = $LineData[8];
	    my $seqfrom = $LineData[10];
	    my $seqto   = $LineData[11];

	    my $key = "$hmmfrom:$hmmto:$seqfrom:$seqto";
	    push(@HitCoords,$key);
	    
	    $line = <$searchout>;
	}

	# Iterate over the domains, pulling in the original alignments
	# and storing them if they're applicable in any of our spliced
	# alignments
	foreach my $key (@HitCoords) {

	    while ($line !~ /\=\= domain/) { $line = <$searchout>; }

	    my $entry = $line;     # Score information
	    $line  = <$searchout>; # HMM Alignment
	    $entry = $entry.$line;
	    $line  = <$searchout>; # Alignment match indicators
	    $entry = $entry.$line;
	    $line  = <$searchout>; # Translation
	    $entry = $entry.$line;
	    $line  = <$searchout>; # DNA
	    $entry = $entry.$line;

	    if ($OrigHitDoms{$key}) {
		my $id = $OrigHitDoms{$key};
		$IDsToOutput{$id} = $entry;
	    }

	}

    }

    close($searchout);

    
    #
    # CHECKPOINT!
    #

    
    # Cool! Now we have everything we need to gather from the search output
    # itself, so it's time to make all this content presentable
    foreach my $hitkey (keys %SplicedHits) {

	# Grab the original hit outputs
	$hitkey =~ /(\d+\:\d+)\|(\d+\:\d+)/;
	my $id1 = $1;
	my $id2 = $2;
	my $origoutput1 = $IDsToOutput{$id1};
	my $origoutput2 = $IDsToOutput{$id2};

	# Grab our spliced alignment data
	my $hitentry  = $SplicedHits{$hitkey};
	my @EntryData = split(/\|/,$hitentry);

	my $hmm1from  = $EntryData[0];
	my $hmm1to    = $EntryData[1];
	my $dna1from  = $EntryData[2];
	my $dna1to    = $EntryData[3];

	my $hmm2from  = $EntryData[4];
	my $hmm2to    = $EntryData[5];
	my $dna2from  = $EntryData[6];
	my $dna2to    = $EntryData[7];

	my $chromosome = $EntryData[8];
	my $score      = $EntryData[9];

	my $revcomp = 0;
	$revcomp    = 1 if ($dna1from > $dna1to);

	# Pull in the DNA sequences, with a little bit extra (i.e., the SS dinucls)
	my $eslsfetch;

	my $adjusted_to = $dna1to;
	if ($revcomp) { $adjusted_to -= 2; }
	else          { $adjusted_to += 2; }
	$eslsfetch = "esl-sfetch -c $dna1from\.\.$adjusted_to";
	$eslsfetch = $eslsfetch." \"$dnafile\" \"$chromosome\"";
	my $dna1   = EslSfetch($eslsfetch);

	my @DNAChars1 = split(//,$dna1);
	my $dnalen1 = scalar(@DNAChars1);

	# We'll make the splice signal characters lower-case
	$DNAChars1[$dnalen1-2] = lc($DNAChars1[$dnalen1-2]);
	$DNAChars1[$dnalen1-1] = lc($DNAChars1[$dnalen1-1]);
	
	my $adjusted_from = $dna2from;
	if ($revcomp) { $adjusted_from += 2; }
	else          { $adjusted_from -= 2; }
	$eslsfetch = "esl-sfetch -c $adjusted_from\.\.$dna2to";
	$eslsfetch = $eslsfetch." \"$dnafile\" \"$chromosome\"";
	my $dna2   = EslSfetch($eslsfetch);
	
	my @DNAChars2 = split(//,$dna2);
	my $dnalen2 = scalar(@DNAChars2);

	# We'll make the splice signal characters lower-case
	$DNAChars2[0] = lc($DNAChars2[0]);
	$DNAChars2[1] = lc($DNAChars2[1]);


	# Great! Now that we have our DNA sequences, we can think about aligning
	# things in a visually approachable way!
	my ($model1,$match1,$trans1,$nucls1,$num_spare_nucls1) 
	    = VisAli($origoutput1,\@ConsChars,$hmm1from,$hmm1to,\@DNAChars1,0);
	my ($model2,$match2,$trans2,$nucls2,$num_spare_nucls2) 
	    = VisAli($origoutput2,\@ConsChars,$hmm2from,$hmm2to,\@DNAChars2,2);

	# Split everything into arrays
	my @Model1 = split(//,$model1);
	my @Match1 = split(//,$match1);
	my @Trans1 = split(//,$trans1);
	my @Nucls1 = split(//,$nucls1);

	my @Model2 = split(//,$model2);
	my @Match2 = split(//,$match2);
	my @Trans2 = split(//,$trans2);
	my @Nucls2 = split(//,$nucls2);

	
	# Our final little task is identifying whether an amino has been split
	# over the splice boundary and recording it.
	#
	# We'll split things between the two ways this can happen (aa|b) vs (a|bb)
	# for simplified coding.
	#
	if ($num_spare_nucls1 == 2) {

	    my $codon = $Nucls1[scalar(@Nucls1)-4].$Nucls1[scalar(@Nucls1)-3];
	    $codon    = $codon.$Nucls2[2];
	    my $trans = TranslateCodon($codon);

	    $Trans1[scalar(@Nucls1)-3] = $trans;
	    $Model1[scalar(@Nucls1)-3] = $ConsChars[$hmm1to];
	    
	} elsif ($num_spare_nucls2 == 2) {

	    my $codon = $Nucls1[scalar(@Nucls1)-3];
	    $codon    = $codon.$Nucls2[2].$Nucls2[3];
	    my $trans = TranslateCodon($codon);

	    $Trans2[2] = $trans;
	    $Model2[2] = $ConsChars[$hmm2from-1];

	}

	# Join 'em
	push(@Model1,'|');
	push(@Model1,@Model2);

	push(@Match1,'|');
	push(@Match1,@Match2);

	push(@Trans1,'|');
	push(@Trans1,@Trans2);

	push(@Nucls1,'|');
	push(@Nucls1,@Nucls2);


	# Format things nicely
	my ($h1f,$h2f) = FormatForOutput($hmm1from,$hmm2from);
	my ($h1t,$h2t) = FormatForOutput($hmm1to,$hmm2to);
	my ($d1f,$d2f) = FormatForOutput($dna1from,$dna2from);
	my ($d1t,$d2t) = FormatForOutput($dna1to,$dna2to);

	
	# PRINT PRINT PRINT!
	print "\n";
	for (my $i=0; $i<$linewidth+11; $i++) { print "-"; }
	print "\n\n\n";
	print "  HMM  :  $protfile\n";
	print "  DNA  :  $dnafile\n";
	print "  Chr. :  $chromosome\n";
	print "  Hit 1:  $h1f-$h1t  ($d1f..$d1t)\n";
	print "  Hit 2:  $h2f-$h2t  ($d2f..$d2t)\n";
	print "  Score:  $score\n";
	print "\n\n";
	for (my $i=0; $i<scalar(@Nucls1); $i+=$linewidth) {

	    print "  Model:   ";
	    for (my $j=0; $j<$linewidth && $i+$j<scalar(@Nucls1); $j++) {
		print "$Model1[$i+$j]";
	    }
	    print "\n";

	    print "           ";
	    for (my $j=0; $j<$linewidth && $i+$j<scalar(@Nucls1); $j++) {
		print "$Match1[$i+$j]";
	    }
	    print "\n";

	    print "  Trans:   ";
	    for (my $j=0; $j<$linewidth && $i+$j<scalar(@Nucls1); $j++) {
		print "$Trans1[$i+$j]";
	    }
	    print "\n";

	    print "  Nucls:   ";
	    for (my $j=0; $j<$linewidth && $i+$j<scalar(@Nucls1); $j++) {
		print "$Nucls1[$i+$j]";
	    }
	    print "\n";
	    print "\n\n";

	}

	print "\n" if (scalar(@Nucls1) % $linewidth);
	
	
    }

    # Pizza!
    
}

    

1;  ##   END OF SCRIPT   ##





############################
#                          #
#  FUNCTIONAL SUBROUTINES  #
#                          #
############################





#################################################################
#
#  FUNCTION:  PrintUsage
#
sub PrintUsage
{
    print "\n";
    print "  USAGE:  ./RunAndParseHmmsearcht.pl [search-pairs-file]\n";
    print "\n";
    print "  WHERE:  The \"search-pairs-file\" contains whitespace-separated pairs\n";
    print "          of protein HMM file locations and DNA target sequence files\n";
    print "          to use as inputs to hmmsearcht.\n";
    die   "\n";
}




#################################################################
#
#  FUNCTION:  PullConsensusSeq
#
sub PullConsensusSeq
{
    my $hmmfname = shift;
    my $inf = OpenInputFile($hmmfname);

    # Eat the header
    my $line;
    while ($line = <$inf>) { last if ($line =~ /\s+1\s+\d+\.\d+/); }

    # Run through the sequence
    my $consensus = '';
    while (!eof($inf)) {
	$line =~ /\d+\s+(\S)\s+\S+\s+\S+\s+\S+\s*$/;
	$consensus = $consensus.uc($1);
	$line = <$inf>;
	$line = <$inf>;
	$line = <$inf>;
	last if ($line =~ /\/\//);
    }

    close($inf);
    
    return $consensus;
}




#################################################################
#
#  FUNCTION:  EslSfetch
#
sub EslSfetch
{
    my $cmd  = shift;
    my $inf  = OpenSystemCommand($cmd);
    my $line = <$inf>;
    my $seq  = "";
    while ($line = <$inf>) {
	$line =~ s/\n|\r|\s//g;
	$seq = $seq.uc($line);
    }
    close($inf);
    return $seq;
}




#################################################################
#
#  FUNCTION:  VisAli
#
sub VisAli
{
    my $orig_output = shift;

    my $cons_chars_ref = shift;
    my @ConsChars = @{$cons_chars_ref};

    my $hmmfrom = shift;
    my $hmmto   = shift;

    my $dna_chars_ref = shift;
    my @DNAChars = @{$dna_chars_ref};

    # This is used to indicate whether we're working with the upstream exon (0) or
    # downstream exon (2).
    my $starts_at = shift;

    # Note that because we say the final amino acid is whichever amino the final coding
    # nucleotides belong to, we may have an "hmmto - hmmfrom" that gives us one more
    # amino than can wholly be attributed to this sequence.
    my $num_trans_aminos = int((scalar(@DNAChars)-2) / 3);
    my $num_spare_nucls  = scalar(@DNAChars) - 2 - (3 * $num_trans_aminos);

    # Is there an uglier way to make sure we don't grab the splice signal? I doubt it.
    my @TransChars;
    if ($starts_at == 0) {
	# UPSTREAM
	for (my $i=0; $i<3*$num_trans_aminos; $i+=3) {
	    my $codon = $DNAChars[$i].$DNAChars[$i+1].$DNAChars[$i+2];
	    my $amino = TranslateCodon($codon);
	    push(@TransChars,$amino);
	}
    } else {
	# DOWNSTREAM
	for (my $i=2+$num_spare_nucls; $i<scalar(@DNAChars); $i+=3) {
	    my $codon = $DNAChars[$i].$DNAChars[$i+1].$DNAChars[$i+2];
	    my $amino = TranslateCodon($codon);
	    push(@TransChars,$amino);
	}
    }

    # Now the "fun" task of parsing the original HMMER output!
    # These are: 0. Scoring info.
    #            1. HMM
    #            2. Quality/Matches
    #            3. Translation
    #            4. DNA
    my @OrigLines = split(/\n/,$orig_output);
    $OrigLines[1] =~ /^\s+\S+\s+(\d+)(\D+)(\d+)/;
    my $orig_hmmfrom = $1;
    my $orig_hmmali = uc($2);
    my $orig_hmmto = $3;
    $orig_hmmali =~ s/\s+//g;
    my @OrigAli = split(//,$orig_hmmali);

    my $Model = '';
    my $Trans = '';
    my $Nucls = '';
    my $Match = '';

    # We're going to need to take pretty different approaches if we're working with the
    # upstream or downstream exon, due to the orientation of the model relative to the
    # rest of the spliced hit...

    # Note that we'll fill in with consensus characters anywhere that we extend past
    # what's given by the alignment.
    
    # UPSTREAM
    if ($starts_at == 0) {

	my $aa_pos = 0;
	my $dna_pos = 0;
	
	while ($aa_pos < MIN(scalar(@TransChars),scalar(@OrigAli))) {
	    $Model = $Model.' '.$OrigAli[$aa_pos].' ';
	    $Trans = $Trans.' '.$TransChars[$aa_pos].' ';
	    if ($OrigAli[$aa_pos] eq $TransChars[$aa_pos]) { 
		$Match = $Match.' | ';
	    } else {
		$Match = $Match.'   ';
	    }
	    $aa_pos++;
	    my $codon=$DNAChars[$dna_pos].$DNAChars[$dna_pos+1].$DNAChars[$dna_pos+2];
	    $Nucls = $Nucls.$codon;
	    $dna_pos += 3;
	}

	while ($aa_pos < $num_trans_aminos) {
	    $Model = $Model.' '.$ConsChars[$hmmto-($num_trans_aminos-$aa_pos)].' ';
	    $Trans = $Trans.' '.$TransChars[$aa_pos].' ';
	    if ($ConsChars[$hmmto-($num_trans_aminos-$aa_pos)] eq $TransChars[$aa_pos]) { 
		$Match = $Match.' | ';
	    } else {
		$Match = $Match.'   ';
	    }
	    $aa_pos++;
	    my $codon=$DNAChars[$dna_pos].$DNAChars[$dna_pos+1].$DNAChars[$dna_pos+2];
	    $Nucls = $Nucls.$codon;
	    $dna_pos += 3;
	}

	while ($dna_pos < scalar(@DNAChars)) {
	    $Model = $Model.' ';
	    $Trans = $Trans.' ';
	    $Match = $Match.' ';
	    $Nucls = $Nucls.$DNAChars[$dna_pos++];
	}


    }

    # DOWNSTREAM
    if ($starts_at == 2) {

	# We may need to add characters from the consensus sequence to fill things out
	if (scalar(@TransChars) > scalar(@OrigAli)) {
	    my @ConsFill;
	    for (my $i=0; $i<scalar(@TransChars)-scalar(@OrigAli); $i++) {
		push(@ConsFill,$ConsChars[$hmmfrom-1+$i]);
	    }
	    push(@ConsFill,@OrigAli);
	    @OrigAli = @ConsFill;
	} 

	if (scalar(@TransChars) < scalar(@OrigAli)) {
	    my $offset = scalar(@OrigAli)-scalar(@TransChars);
	    for (my $i=0; $i<scalar(@TransChars); $i++) {
		$OrigAli[$i] = $OrigAli[$i+$offset];
	    }
	}

	my $aa_pos = scalar(@TransChars);
	my $dna_pos = scalar(@DNAChars)-1;

	while ($aa_pos) {
	    $aa_pos--;
	    $Model = ' '.$OrigAli[$aa_pos].' '.$Model;
	    $Trans = ' '.$TransChars[$aa_pos].' '.$Trans;
	    if ($OrigAli[$aa_pos] eq $TransChars[$aa_pos]) {
		$Match = ' | '.$Match;
	    } else {
		$Match = '   '.$Match;
	    }
	    my $codon=$DNAChars[$dna_pos-2].$DNAChars[$dna_pos-1].$DNAChars[$dna_pos];
	    $dna_pos -= 3;
	    $Nucls = $codon.$Nucls;
	}

	while ($dna_pos >= 0) {
	    $Model = ' '.$Model;
	    $Trans = ' '.$Trans;
	    $Match = ' '.$Match;
	    $Nucls = $DNAChars[$dna_pos--].$Nucls;
	}
	
	
    }
    
    return ($Model,$Match,$Trans,$Nucls,$num_spare_nucls);
    
}





#########################################################################
#
#  Function Name: TranslateCodon
#
#  About: Convert a DNA triple to an amino acid. ('x' for 'stop')
#
sub TranslateCodon
{
    
    my $codon = shift;
    $codon    = uc($codon);
    $codon    =~ s/U/T/g;
    return 'X' if ($codon =~ /[^ACGT]/);
    
    my @codonArray = split('',$codon);
    
    if ($codonArray[0] eq 'A') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "AAA") { return 'K'; }
	    if ($codon eq "AAC") { return 'N'; }
	    if ($codon eq "AAG") { return 'K'; }
	    if ($codon eq "AAT") { return 'N'; }
	}
	if ($codonArray[1] eq 'C') {
	    if ($codon eq "ACA") { return 'T'; }
	    if ($codon eq "ACC") { return 'T'; }
	    if ($codon eq "ACG") { return 'T'; }
	    if ($codon eq "ACT") { return 'T'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "AGA") { return 'R'; }
	    if ($codon eq "AGC") { return 'S'; }
	    if ($codon eq "AGG") { return 'R'; }
	    if ($codon eq "AGT") { return 'S'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "ATA") { return 'I'; }
	    if ($codon eq "ATC") { return 'I'; }
	    if ($codon eq "ATG") { return 'M'; }
	    if ($codon eq "ATT") { return 'I'; }
	}
    }

    if ($codonArray[0] eq 'C') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "CAA") { return 'Q'; }
	    if ($codon eq "CAC") { return 'H'; }
	    if ($codon eq "CAG") { return 'Q'; }
	    if ($codon eq "CAT") { return 'H'; }
	}
	if ($codonArray[1] eq 'C') {
	    if ($codon eq "CCA") { return 'P'; }
	    if ($codon eq "CCC") { return 'P'; }
	    if ($codon eq "CCG") { return 'P'; }
	    if ($codon eq "CCT") { return 'P'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "CGA") { return 'R'; }
	    if ($codon eq "CGC") { return 'R'; }
	    if ($codon eq "CGG") { return 'R'; }
	    if ($codon eq "CGT") { return 'R'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "CTA") { return 'L'; }
	    if ($codon eq "CTC") { return 'L'; }
	    if ($codon eq "CTG") { return 'L'; }
	    if ($codon eq "CTT") { return 'L'; }
	}
    }
    
    if ($codonArray[0] eq 'G') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "GAA") { return 'E'; }
	    if ($codon eq "GAC") { return 'D'; }
	    if ($codon eq "GAG") { return 'E'; }
	    if ($codon eq "GAT") { return 'D'; }
	}
	if ($codonArray[1] eq 'C') {    
	    if ($codon eq "GCA") { return 'A'; }
	    if ($codon eq "GCC") { return 'A'; }
	    if ($codon eq "GCG") { return 'A'; }
	    if ($codon eq "GCT") { return 'A'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "GGA") { return 'G'; }
	    if ($codon eq "GGC") { return 'G'; }
	    if ($codon eq "GGG") { return 'G'; }
	    if ($codon eq "GGT") { return 'G'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "GTA") { return 'V'; }
	    if ($codon eq "GTC") { return 'V'; }
	    if ($codon eq "GTG") { return 'V'; }
	    if ($codon eq "GTT") { return 'V'; }
	}
    }
    if ($codonArray[0] eq 'T') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "TAA") { return 'x'; }
	    if ($codon eq "TAC") { return 'Y'; }
	    if ($codon eq "TAG") { return 'x'; }
	    if ($codon eq "TAT") { return 'Y'; }
	}
	if ($codonArray[1] eq 'C') {
	    if ($codon eq "TCA") { return 'S'; }
	    if ($codon eq "TCC") { return 'S'; }
	    if ($codon eq "TCG") { return 'S'; }
	    if ($codon eq "TCT") { return 'S'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "TGA") { return 'x'; }
	    if ($codon eq "TGC") { return 'C'; }
	    if ($codon eq "TGG") { return 'W'; }
	    if ($codon eq "TGT") { return 'C'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "TTA") { return 'L'; }
	    if ($codon eq "TTC") { return 'F'; }
	    if ($codon eq "TTG") { return 'L'; }
	    if ($codon eq "TTT") { return 'F'; }
	}
    }

    # Weird codon is weird. TO THE BIN WITH YOU!
    return 'X';
    
}





#################################################################
#
#  FUNCTION:  FormatForOutput
#
sub FormatForOutput
{
    my $a = shift;
    my $b = shift;
    while (length("$a")<length("$b")) { $a = ' '.$a; }
    while (length("$b")<length("$a")) { $b = ' '.$b; }
    return ($a,$b);
}








###############################
#                             #
#  BEUREAUCRATIC SUBROUTINES  #
#                             #
###############################



#################################################################
#
#  FUNCTION:  MIN
#
sub MIN
{
    my $a = shift;
    my $b = shift;
    return $a if ($a < $b);
    return $b;
}



#################################################################
#
#  FUNCTION:  MAX
#
sub MAX
{
    my $a = shift;
    my $b = shift;
    return $a if ($a > $b);
    return $b;
}



#################################################################
#
#  FUNCTION:  RunSystemCommand
#
sub RunSystemCommand
{
    my $command = shift;
    if (system($command)) { die "\n  ERROR:  System command '$command' failed during execution\n\n"; }
}



#################################################################
#
#  FUNCTION:  OpenSystemCommand
#
sub OpenSystemCommand
{
    my $command = shift;
    if ($command !~ /\s+\|\s*$/) { $command = $command.' |'; }
    open(my $command_output,$command) || die "\n  ERROR:  Failed to open output from system command '$command'\n\n";
    return $command_output;
}



#################################################################
#
#  FUNCTION:  OpenInputFile
#
sub OpenInputFile
{
    my $filename = shift;
    if (!(-e $filename)) { die "\n  ERROR:  Failed to locate input file '$filename'\n\n"; }
    open(my $filehandle,'<',$filename) || die "\n  ERROR:  Failed to open input file '$filename'\n\n";
    return $filehandle;
}



#################################################################
#
#  FUNCTION:  NameOutputFile
#
sub NameOutputFile
{
    my $intended_name = shift;
    my $basename;
    my $extension;
    if ($intended_name =~ /(\S+)(\.[^\.]+)$/) {
	$basename = $1;
	$extension = $2;
    } else {
	$basename = $intended_name;
	$extension = '';
    }
    my $filename = $basename.$extension;
    my $attempt = 1;
    while (-e $filename) {
	$attempt++;
	$filename = $basename.'_'.$attempt.$extension;
    }
    return $filename;
}



#################################################################
#
#  FUNCTION:  OpenOutputFile
#
sub OpenOutputFile
{
    my $filename = shift;
    $filename = NameOutputFile($filename);
    open(my $filehandle,'>',$filename) || die "\n  ERROR:  Failed to open output file '$filename'\n\n";
    return $filehandle;
}



#################################################################
#
#  FUNCTION:  OpenFileToAppendTo
#
sub OpenFileToAppendTo
{
    my $filename = shift;
    open(my $filehandle,'>>',$filename) || die "\n  ERROR:  Failed to open output file '$filename' (for appending)\n\n";
    return $filehandle;
}



#################################################################
#
#  FUNCTION:  ConfirmDirectory
#
sub ConfirmDirectory
{
    my $dirname = shift;
    if (!(-d $dirname)) { die "\n  ERROR:  Failed to locate directory '$dirname'\n\n"; }
    if ($dirname !~ /\/$/) { $dirname = $dirname.'/'; }
    return $dirname;
}



#################################################################
#
#  FUNCTION:  DeleteFile
#
sub DeleteFile
{
    my $filename;
    if (-e $filename) {
	my $rm_command = 'rm '.$filename;
	RunSystemCommand($rm_command);
    }
}



#################################################################
#
#  FUNCTION:  OpenDirectory
#
sub OpenDirectory
{
    my $dirname = shift;
    $dirname = ConfirmDirectory($dirname);
    opendir(my $dirhandle,$dirname) || die "\n  ERROR:  Failed to open directory '$dirname'\n\n";
    return $dirhandle;
}



#################################################################
#
#  FUNCTION:  NameDirectory
#
sub NameDirectory
{
    my $intended_name = shift;
    $intended_name =~ s/\/$//;
    my $dirname = $intended_name;
    my $attempt = 1;
    while (-d $dirname) {
	$attempt++;
	$dirname = $intended_name.'_'.$attempt;
    }
    $dirname = $dirname.'/';
    return $dirname;
}



#################################################################
#
#  FUNCTION:  CreateDirectory
#
sub CreateDirectory
{
    my $dirname = shift;
    $dirname = NameDirectory($dirname);
    RunSystemCommand("mkdir $dirname");
    if (!(-d $dirname)) { die "\n  ERROR:  Creation of directory '$dirname' failed\n\n"; }
    return $dirname;
}


##   END OF FILE   ##













