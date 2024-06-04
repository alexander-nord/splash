#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub SetupMiniprot;
sub GetLocalZlib;
sub CopyMiniprotMakefile;



if (@ARGV == 1 && lc($ARGV[0]) =~ /^\-?\-?miniprot/)
{
    SetupMiniprot();
    exit(0);
}



my $splash_bathsearch = 'bathsearch.c';
if (!(-e $splash_bathsearch)) {
    die "\n  ERROR: I can't find bathsearch.c!\n\n";
}


my $splash_pipeline = 'inc/p7_pipeline.c';
if (!(-e $splash_pipeline)) {
    die "\n  ERROR: There should be a 'p7_pipeline.c' file in 'inc/'? (mainly force-disables frameshift-aware pipeline for splash)\n\n";
}


my $splash_hmmer_h = 'inc/hmmer.h';
if (!(-e $splash_hmmer_h)) {
    die "\n  ERROR: There should be a 'hmmer.h' file in 'inc/'\n\n";
}


my $bath_dir_name = 'BATH/';


if (!(-d $bath_dir_name)) {

    my $bath_git_url = 'https://github.com/traviswheelerlab/BATH';
    if (system("git clone $bath_git_url")) {
	die "\n  ERROR: Failed to clone into $bath_git_url\n\n";
    }
    
    chdir($bath_dir_name);

    system("git clone https://github.com/TravisWheelerLab/easel");
    chdir("easel/");
    system("git checkout BATH");
    chdir("../");

    system("autoconf");
    system("./configure");

    system("cp ../$splash_bathsearch src/");
    system("cp ../$splash_pipeline   src/");
    system("cp ../$splash_hmmer_h    src/");
    system("make");
    
} else {

    my $bath_src_name = $bath_dir_name.'src/';
    system("cp $splash_bathsearch $bath_src_name");
    system("cp $splash_pipeline   $bath_src_name");
    system("cp $splash_hmmer_h    $bath_src_name");
    chdir($bath_src_name);
    system("make -B");

}

1;






sub SetupMiniprot
{

    if (-d 'miniprot/')
    {
        print "  Miniprot directory already exists -- remove to re-install\n";
        return;
    }


    # We're just going to do a local install of zlib
    my $zlib_location = GetLocalZlib();

    system("git clone https://github.com/lh3/miniprot");
    CopyMiniprotMakefile("inc/MiniprotMakefile","miniprot/Makefile",$zlib_location);
    
    chdir("miniprot");
    system("make");
    chdir("..");
    
}


sub GetLocalZlib
{
    system("wget http://zlib.net/current/zlib.tar.gz");
    system("gunzip zlib.tar.gz");
    system("tar -xf zlib.tar");
    system("rm zlib.tar");

    opendir(my $CurrentDir,'.');
    while (my $sub_dir_name = readdir($CurrentDir))
    {
        if (lc($sub_dir_name) =~ /^zlib/)
        {
            system("mv $sub_dir_name zlib") 
                if ($sub_dir_name ne 'zlib');
            last;
        }
    }
    closedir($CurrentDir);

    chdir('zlib');
    system('./configure');
    system('make');
    chdir('..');

    open(my $PWD,'pwd |');
    my $pwd = <$PWD>;
    $pwd =~ s/\n|\r//g;
    close($PWD);

    $pwd = $pwd.'/' if ($pwd !~ /\/$/);
    return $pwd.'zlib/';

}


sub CopyMiniprotMakefile
{
    my $file_to_copy  = shift;
    my $target_file   = shift;
    my $zlib_location = shift;

    open(my  $InFile,'<',$file_to_copy);
    open(my $OutFile,'>',$target_file );

    while (my $line = <$InFile>)
    {
        $line =~ s/\n|\r//g;
        $line =~ s/<<SPLASH_LOCAL_ZLIB>>/$zlib_location/g;
        print $OutFile "$line\n";
    }

    close( $InFile);
    close($OutFile);
}





