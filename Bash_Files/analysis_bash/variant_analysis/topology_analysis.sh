#!/bin/bash
#SBATCH	-p single
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --mem=5000
#SBATCH --job-name="topology_analysis"
#SBATCH --output=topology_analysis.out
#SBATCH --mail-user=p.burmedi@stud.uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_pburmedi/repos/ligrec_robustness


/net/data.isilon/ag-saez/bq_pburmedi/SOFTWARE/miniconda3/envs/liana_env/bin/Rscript Code/Analysis_Scripts/Supplementary_distribution_of_top_ranks.R