#!/bin/bash -l
#SBATCH --mem=50G
#SBATCH -t 48:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user rvnair@stanford.edu
#SBATCH -A baas_lab_mpsnyder
#SBATCH -o %x.o%j
date
hostname
module load bcl2fastq2/2.20.0
module load cellranger/2.2.0
cellranger mkfastq --run=/BaaS/labs/mpsnyder/Rabinovitch/Toshie/180924_COOPER_0219_AHW373BBXX_L8_TS-D/raw_data --csv=/BaaS/labs/mpsnyder/Rabinovitch/Toshie/docs/TS-D.csv --use-bases-mask=y26n*,i8,n*,y100n* --ignore-dual-index
echo 'Done'
