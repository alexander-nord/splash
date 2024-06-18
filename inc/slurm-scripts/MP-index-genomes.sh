#!/usr/bin/bash
#SBATCH --job-name=mpindex
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --partition=standard
#SBATCH --account=twheeler

cd $HOME/splash/miniprot

./miniprot -t4 -d /xdisk/twheeler/alexnord/inputs-to-splash/genomes/hg38.mpi    /xdisk/twheeler/alexnord/inputs-to-splash/genomes/hg38.fa
./miniprot -t4 -d /xdisk/twheeler/alexnord/inputs-to-splash/genomes/galGal6.mpi /xdisk/twheeler/alexnord/inputs-to-splash/genomes/galGal6.fa
./miniprot -t4 -d /xdisk/twheeler/alexnord/inputs-to-splash/genomes/gasAcu1.mpi /xdisk/twheeler/alexnord/inputs-to-splash/genomes/gasAcu1.fa
./miniprot -t4 -d /xdisk/twheeler/alexnord/inputs-to-splash/genomes/dm6.mpi     /xdisk/twheeler/alexnord/inputs-to-splash/genomes/dm6.fa

