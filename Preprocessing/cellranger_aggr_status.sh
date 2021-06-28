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
cellranger aggr --id=AGG_STATUS_TS-D --csv=/BaaS/labs/marlener/Toshie/docs/TS-D-STATUS-AGGR.csv --normalize=mapped
echo 'Done'
