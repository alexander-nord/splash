#!/usr/bin/bash
#SBATCH --job-name=splanther-pack
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --partition=standard
#SBATCH --account=twheeler

cd /xdisk/twheeler/alexnord/results/
tar -czf Splash-PANTHER-Output.tgz Splash-PANTHER-Output
