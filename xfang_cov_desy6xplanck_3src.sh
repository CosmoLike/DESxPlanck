#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -W group_list=cosmo
#PBS -q standard
#PBS -J 1-630
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l place=pack:shared
#PBS -l walltime=7:00:00
#PBS -N des6xplanck_cov
#PBS -e /home/u1/xfang/output/
#PBS -o /home/u1/xfang/output/

module load gsl/2/2.1

cd $PBS_O_WORKDIR
/home/u1/xfang/CosmoLike/DESxPlanck/compute_covariances_fourier $PBS_ARRAY_INDEX ini_files/cov_y6_mcal_fourier_3src.ini >&/home/u1/xfang/output/job_output_$PBS_ARRAY_INDEX.log