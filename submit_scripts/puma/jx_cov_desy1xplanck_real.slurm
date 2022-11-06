#!/bin/bash
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------

###--- Optional. Set the job name
#SBATCH --job-name=Y1_6x2pt_cov

###--- Optional. Set the output filename.
###--- SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=/home/u17/jiachuanxu/output/%x-%j.out
#SBATCH --error=/home/u17/jiachuanxu/output/%x-%j.err

###--- Optionan. Array job submission
#SBATCH --array=1-770

###--- REQUIRED. Specify the PI group for this job
#SBATCH --account=timeifler

###--- Optional. Request email when job begins and ends
#SBATCH --mail-type=ALL

###--- Optional. Specify email address to use for notification
#SBATCH --mail-user=jiachuanxu@email.arizona.edu

###--- REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
###SBATCH --qos=user_qos_timeifler

###--- REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks-per-node=2

###--- REQUIRED. Set the number of nodes
#SBATCH --nodes=1

###--- REQUIRED. Set the memory-per-node required for this job.
#SBATCH --mem=16gb

###--- REQUIRED. Specify the walltime required for this job, hhh:mm:ss
#SBATCH --time=60:00:00


# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
### Load required modules/libraries if needed
#module load gnu/5.4.0
module load gsl/2.6
module load python/3.6
#module load mpich/3.3.1
module load openmpi3/3.1.4

### change to your script submission directory
cd $SLURM_SUBMIT_DIR

### Run your work
### run your executable program with begin and end date and time output

### If run full covariance matrix
### =============================
for (( c=0; c<2; c++ ))
do
    hit=$(( ${SLURM_ARRAY_TASK_ID} + $c * 770 ))
    ./compute_covariances_real_6x2pt $hit ini_files/cov_y1_mcal_mix_planck_DR4.ini ./cmblensrec/plancksmica/cmb_lmax3000.txt >&/home/u17/jiachuanxu/output/job_output_${hit}.log
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

