#!/bin/bash -l
#SBATCH --mem=50G
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user rvnair@stanford.edu
#SBATCH -A baas_lab_mpsnyder
#SBATCH -o %x.o%j
date
hostname
module load cellranger/2.2.0
cellranger aggr --id=AGG_TS-D --csv=/BaaS/labs/mpsnyder/Rabinovitch/Toshie/docs/TS-D-AGGR.csv --normalize=mapped
echo 'Done'
