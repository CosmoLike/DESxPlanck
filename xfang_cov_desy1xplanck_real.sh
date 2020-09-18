#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -W group_list=cosmo
#PBS -q high_pri
#PBS -J 1-1540
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l place=pack:shared
#PBS -l walltime=7:00:00
#PBS -N des1xplanck_cov
#PBS -e /home/u1/xfang/output/
#PBS -o /home/u1/xfang/output/

module load gsl/2/2.1

cd $PBS_O_WORKDIR
/home/u1/xfang/CosmoLike/DESxPlanck/compute_covariances_real_6x2pt $PBS_ARRAY_INDEX ini_files/cov_y1_mcal_mix.ini >&/home/u1/xfang/output/job_output_$PBS_ARRAY_INDEX.log