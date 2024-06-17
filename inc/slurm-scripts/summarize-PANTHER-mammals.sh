#!/usr/bin/bash
#SBATCH --job-name=splanther-summary
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --partition=standard
#SBATCH --account=twheeler


SPECIES_GUIDE="$HOME/data/PANTHER-species-guide"
IN_DIR_NAME="/xdisk/twheeler/alexnord/mammal-PANTHER"
OUT_DIR_NAME="/xdisk/twheeler/alexnord/results/Splash-PANTHER-Output"

cd $HOME/splash

time perl Run-Splash-Test.pl --summarize $OUT_DIR_NAME $IN_DIR_NAME $SPECIES_GUIDE 1>summarize.out 2>summarize.err
