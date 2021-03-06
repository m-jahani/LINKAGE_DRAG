#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for BayPass as job arrays
# ---------------------------------------------------------------------
#SBATCH --account=def-rieseber
#SBATCH --job-name=GenomicPredictionModel
#SBATCH --cpus-per-task=48
#SBATCH --mem=192000M
#SBATCH --time=02:59:00
#SBATCH --output=%j_%x.out
#SBATCH --array=1-46
# ---------------------------------------------------------------------
echo "Current working directory: $(pwd)"
echo "Starting run at: $(date)"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""

#Load the conda environment
#source /home/mjahani/miniconda3/bin/activate rrblup

module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2

#input data
ROOT_DIR=/home/mjahani/scratch/LINKAGE_DRAG/new_data
SCRIPT=${ROOT_DIR}/CC_GP_Permutation.R
GENOTYPE=${ROOT_DIR}/SAM_introgression_donor_ANNUUS_maf0.03_rrBLUP.in
PHENOTYPE=${ROOT_DIR}/phenotype_common_georgia_corrected.csv
DONOR=WILD_ANNUUS
# Output directory
SAVE_DIR=${ROOT_DIR}/result

#run th escript
Rscript $SCRIPT $GENOTYPE $PHENOTYPE $SAVE_DIR $SLURM_ARRAY_TASK_ID $DONOR

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: $(date)"
# ---------------------------------------------------------------------
