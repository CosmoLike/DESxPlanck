#!/bin/bash
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------

###--- Optional. Set the job name
#SBATCH --job-name=ccoce

###--- Optional. Set the output filename.
###--- SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=/xdisk/timeifler/jiachuanxu/cosmolike_chains/outputs/%x-%j.out
#SBATCH --error=/xdisk/timeifler/jiachuanxu/cosmolike_chains/outputs/%x-%j.err

###--- REQUIRED. Specify the PI group for this job
#SBATCH --account=timeifler

###--- Optional. Request email when job begins and ends
#SBATCH --mail-type=ALL

###--- Optional. Specify email address to use for notification
#SBATCH --mail-user=jiachuanxu@email.arizona.edu

###--- REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
#SBATCH --qos=qual_qos_timeifler

###--- REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks-per-node=28

###--- REQUIRED. Set the number of nodes
#SBATCH --nodes=20

###--- REQUIRED. Set the memory-per-node required for this job.
#SBATCH --mem=168gb

###--- REQUIRED. Specify the walltime required for this job, hhh:mm:ss
#SBATCH --time=24:00:00
 
 
# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
### Load required modules/libraries if needed
module load gsl
module load anaconda
module load openmpi3

### change to your script submission directory
cd $SLURM_SUBMIT_DIR

### activate your python virtual environment
conda init bash
source ~/.bashrc 
conda activate kltools

### Run your work
### run your executable program with begin and end date and time output
export MPI_DSM_DISTRIBUTE
date
/usr/bin/time mpiexec -n 560 --mca btl tcp,self python run_cosmolike_wrapper_3x2pt.py yaml/DESY1_3x2pt_code_comparison.yaml
date
