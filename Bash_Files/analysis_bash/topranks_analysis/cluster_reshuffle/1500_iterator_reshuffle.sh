#!/bin/bash
#SBATCH	-p single
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --mem=20000
#SBATCH --job-name="1500_reshuffle_iterator"
#SBATCH --output=1500_reshuffle_iterator.out
#SBATCH --mail-user=p.burmedi@stud.uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_pburmedi/repos/ligrec_robustness


# Add arguments to pass to Run_Iterator. 1 = reshuffle_or_subset, 2 = top_n, 3 = job_id
/net/data.isilon/ag-saez/bq_pburmedi/SOFTWARE/miniconda3/envs/liana_env/bin/Rscript Code/Analysis_Scripts/Clusters_Run_Iterator.R "reshuffle" 1500 $SLURM_JOBID  
