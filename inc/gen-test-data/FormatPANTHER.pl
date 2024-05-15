#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;



if (@ARGV != 2) { die "\n  USAGE:  ./FormatPANTHER.pl [path/to/PANTHER18.0/] [Out-Dir-Name]\n\n"; }



my $panther_in_dir_name = $ARGV[0];
$panther_in_dir_name = $panther_in_dir_name.'/' if ($panther_in_dir_name !~ /\//);

my $out_dir_name = $ARGV[1];
$out_dir_name = $out_dir_name.'/' if ($out_dir_name !~ /\//);



die "\n  ERROR:  Failed to locate input PANTHER dataset '$panther_in_dir_name'\n\n"
	if (!(-d $panther_in_dir_name));

my $panther_books_dir_name = $panther_in_dir_name.'books/';
die "\n  ERROR:  Failed to locate expected 'books' directory '$panther_books_dir_name'\n\n"
	if (!(-d $panther_books_dir_name));



die "\n  ERROR:  Output directory '$out_dir_name' already exists\n\n"
	if (-d $out_dir_name);

if (system("mkdir \"$out_dir_name\"")) {
	die "\n  ERROR:  Failed to make output directory '$out_dir_name'\n\n";
}



opendir(my $PantherBooksDir,$panther_books_dir_name)
	|| die "\n  ERROR:  Failed to open 'books' directory '$panther_books_dir_name'\n\n";

while (my $book_name = readdir($PantherBooksDir))
{

	next if ($book_name =~ /^\./);
	$book_name =~ s/\/$//;

	my $book_in_dir_name = $panther_books_dir_name.$book_name.'/';
	next if (!(-d $book_in_dir_name));

	my $book_in_hmm_name = $book_in_dir_name.'hmmer.hmm';
	next if (!(-e $book_in_hmm_name));


	my $book_out_dir_name = $out_dir_name.$book_name.'/';
	die "\n  ERROR:  Failed to create book output directory '$book_out_dir_name'\n\n"
		if (system("mkdir \"$book_out_dir_name\""));

	my $out_hmm_name = $book_out_dir_name.$book_name.'.hmm';
	die "\n  ERROR:  Failed to copy PANTHER HMM file '$book_in_hmm_name' to '$out_hmm_name'\n\n"
		if (system("cp \"$book_in_hmm_name\" \"$out_hmm_name\""));

}

closedir($PantherBooksDir);



# I assume this is always available, but it's not technically necessary
my $panther_names_file_name = $panther_in_dir_name.'globals/names.tab';
if (-e $panther_names_file_name)
{

	open(my $PantherNames,'<',$panther_names_file_name)
		|| die "\n  ERROR:  Failed to open PANTHER name guide '$panther_names_file_name'\n\n";

	my $all_names_list = $out_dir_name.'../PANTHER-Name-Guide.out';
	system("rm \"$all_names_list\"") if (-e $all_names_list);

	while (my $line = <$PantherNames>)
	{
		if ($line =~ /^(\S+)\.mag\.mod\s+(.+)\s*$/)
		{
			
			my $book_name = $1;
			my $book_desc = $2;

			my $book_out_dir_name = $out_dir_name.$book_name.'/';
			next if (!(-d $book_out_dir_name));

			my $book_out_desc_file = $book_out_dir_name.'description.txt';
			system("echo \"$book_desc\" > $book_out_desc_file");

			system("echo \"$book_name     $book_desc\" >> \"$all_names_list\"");

		}
	}

	close($PantherNames);

}





1;
