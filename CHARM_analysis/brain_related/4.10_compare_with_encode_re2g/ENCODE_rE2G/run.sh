# CHARM_pre/CHARM_analysis/brain_related/4.10_compare_with_encode_re2g/ENCODE_rE2G/run.sh
#! /bin/bash
#SBATCH -J re2g 
#SBATCH --partition=fatcomp    
#SBATCH --nodelist=node03    
#SBATCH --cpus-per-task=50                   
#SBATCH --output=/share/home/mwang/shared/mwang/CHARM-seq/ENCODE_rE2G/slurm_log/E2G.%j.out   
#SBATCH --error=/share/home/mwang/shared/mwang/CHARM-seq/ENCODE_rE2G/slurm_log/E2G.%j.err     

# source ~/.bashrc
conda activate abc-env
cd /share/home/mwang/shared/mwang/CHARM-seq/ENCODE_rE2G/
snakemake --unlock
snakemake --cores 50
