#!/bin/bash
#SBATCH --job-name=splanther
#SBATCH --array=1-40
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --partition=standard
#SBATCH --account=twheeler

#
#  NOTE:  BE SURE TO INITIALIZE YOUR OUTPUT DIRECTORIES AHEAD OF TIME!!!!
#         ^ This includes 'Results-by-Gene'!
#

SPECIES_GUIDE="$HOME/data/PANTHER-species-guide"
IN_DIR_NAME="/xdisk/twheeler/alexnord/mammal-PANTHER"
OUT_DIR_NAME="/xdisk/twheeler/alexnord/results/Splash-PANTHER-Output"

cd $HOME/splash

time perl Run-Splash-Test.pl --slurm %a --panther -o $OUT_DIR_NAME -n $SLURM_ARRAY_TASK_MAX $IN_DIR_NAME $SPECIES_GUIDE 1>%a.fresh.out 2>%a.fresh.err

