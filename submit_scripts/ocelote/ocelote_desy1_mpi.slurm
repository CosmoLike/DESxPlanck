#!/bin/bash
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------

###--- Optional. Set the job name
#SBATCH --job-name=Y1_3x2_test

###--- Optional. Set the output filename.
###--- SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=/home/u17/jiachuanxu/output/%x-%j.out
#SBATCH --error=/home/u17/jiachuanxu/output/%x-%j.err

###--- REQUIRED. Specify the PI group for this job
#SBATCH --account=timeifler

###--- Optional. Request email when job begins and ends
#SBATCH --mail-type=ALL

###--- Optional. Specify email address to use for notification
#SBATCH --mail-user=jiachuanxu@email.arizona.edu

###--- REQUIRED. Set the partition for your job.
#SBATCH --partition=standard

###--- REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks-per-node=28

###--- REQUIRED. Set the number of nodes
#SBATCH --nodes=20

###--- REQUIRED. Set the memory-per-node required for this job.
#SBATCH --mem-per-cpu=6gb

###--- REQUIRED. Specify the walltime required for this job, hhh:mm:ss
#SBATCH --time=24:00:00
 
 
# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
### Load required modules/libraries if needed
#module load gnu/5.4.0
module load gsl/2.6
module load python/3.6/3.6.5
#module load mpich/3.3.1
module load openmpi3/3.1.4

### change to your script submission directory
cd $SLURM_SUBMIT_DIR

### activate your python virtual environment
VIRTUAL_ENV="/home/u17/jiachuanxu/python3.6_virtualenv"
#export VIRTUAL_ENV
#_OLD_VIRTUAL_PATH="$PATH"
#PATH="$VIRTUAL_ENV/bin:$PATH"
#export PATH

### Run your work
### run your executable program with begin and end date and time output
export MPI_DSM_DISTRIBUTE
EXEC_DIR="/home/u17/jiachuanxu/CosmoLike/DESxPlanck"
cd $EXEC_DIR
source $VIRTUAL_ENV/bin/activate
date
/usr/bin/time mpiexec -n 560 python run_cosmolike_wrapper_3x2pt.py yaml/DESY1_3x2pt_test_ocelote.yaml
date

#PATH="$_OLD_VIRTUAL_PATH"
#export PATH
