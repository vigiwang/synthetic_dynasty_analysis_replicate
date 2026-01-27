#!/bin/bash
#SBATCH --job-name=vg_bootstrap
#SBATCH --output=boot_main.out
#SBATCH --error=boot_main.err
#SBATCH --account= rcc-uchicago
#SBATCH --partition=caslake
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=45
#SBATCH --mem=178G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=weiqiw@rcc.uchicago.edu

module load R/4.4.1
Rscript /home/weiqiw/synthetic_dynasty_analysis/master_main_results.R
