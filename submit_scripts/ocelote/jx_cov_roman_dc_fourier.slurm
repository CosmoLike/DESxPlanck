#!/bin/bash
#SBATCH --job-name=RomanDCCov
#SBATCH --output=/xdisk/timeifler/jiachuanxu/job_logs/RomanDCCov_%A_%a.out
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1

### Roman 3x2pt DC1 (4950)
### Roman 3x2pt DC1 medium (5671, due to different ggl exclude)
#SBATCH --array=1-38

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiachuanxu@arizona.edu

###--- REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
#SBATCH --qos=qual_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=00:30:00


# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
### Load required modules/libraries if needed
#module load gnu/5.4.0
module load gsl
module load openmpi3

### change to your script submission directory
cd $SLURM_SUBMIT_DIR

### Run your work
### run your executable program with begin and end date and time output

### If run full covariance matrix
### =============================

#INI_FILE=ini_files/cov_Roman_DC1_fourier.ini
INI_FILE=ini_files/cov_Roman_DC1_fourier_medium.ini
for (( c=0; c<150; c++ ))
do
    hit=$(( ${SLURM_ARRAY_TASK_ID} + ${c} * 38 ))
    ./compute_covariances_fourier_3x2pt ${hit} ${INI_FILE}
done

### If only running gk and sk, hit = 1035-1539, which has 504 blocks
### ===============================================================
#hit=$(( ${SLURM_ARRAY_TASK_ID} + 1034 ))
#./compute_covariances_real_6x2pt $hit ini_files/cov_y1_mcal_mix.ini ./cmblensrec/plancksmica/cmb_lmax3000.txt >&/home/u17/jiachuanxu/output/job_output_${hit}.log

### In case some jobs failed, uncomment the following lines
### =======================================================
#TODO_list=(1445)
#hit=${TODO_list[${SLURM_ARRAY_TASK_ID}]}

#hit=1540
#echo ${hit}
#./compute_covariances_real_6x2pt ${hit} ini_files/cov_y1_mcal_mix.ini ./cmblensrec/plancksmica/cmb_lmax3000.txt >&/home/u17/jiachuanxu/output/job_output_${hit}.log

