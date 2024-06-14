#!/bin/bash
#SBATCH --job-name=splanther
#SBATCH --array=1-20
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
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

time perl Run-Splash-Test.pl --slurm $SLURM_ARRAY_TASK_ID --panther -o $OUT_DIR_NAME -n $SLURM_ARRAY_TASK_COUNT $IN_DIR_NAME $SPECIES_GUIDE 1>$SLURM_ARRAY_TASK_ID.fresh.out 2>$SLURM_ARRAY_TASK_ID.fresh.err

