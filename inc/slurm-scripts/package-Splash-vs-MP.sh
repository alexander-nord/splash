#!/usr/bin/bash
#SBATCH --job-name=spl-vs-mp-pack
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --partition=standard
#SBATCH --account=twheeler

cd /xdisk/twheeler/alexnord/results/
tar -czf Splash-vs-MP-Output.tgz Splash-vs-MP-Output
