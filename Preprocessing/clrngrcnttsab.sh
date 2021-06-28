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
cellranger count --id=ST06 --transcriptome=$HG19 --fastqs=/BaaS/labs/mpsnyder/Rabinovitch/Toshie/cellranger_mkfastq/HW373BBXX/outs/fastq_path --sample=ST06 --expect-cells=2000
echo 'Done'
