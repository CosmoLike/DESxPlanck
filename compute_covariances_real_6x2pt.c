#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <string.h>
 
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <fftw3.h>

#include "../cosmolike_core/theory/baryons.h"
#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/emu17/P_cb/emu.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift_spline.c"
#include "../cosmolike_core/theory/halo.c"
#include "../cosmolike_core/theory/HOD.c"
#include "../cosmolike_core/theory/pt.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/CMBxLSS_fourier.c"
#include "../cosmolike_core/theory/cluster.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "../cosmolike_core/theory/covariances_3D.c"
#include "../cosmolike_core/theory/covariances_fourier.c"
#include "../cosmolike_core/theory/covariances_CMBxLSS_fourier.c"
// #include "../cosmolike_core/theory/covariances_real_binned.c"
#include "../cosmolike_core/theory/covariances_binned_simple.c"
// #include "../cosmolike_core/theory/run_covariances_real.c"
#include "../cosmolike_core/theory/run_covariances_real_fullsky.c"
#include "../cosmolike_core/theory/run_covariances_real_fullsky_6x2pt.c"
#include "init_LSSxCMB.c"

// Usage example
// compute_covariances_real_6x2pt $PBS_ARRAY_INDEX ini_files/cov_y1_mcal_mix.ini >&/home/u1/xfang/output/job_output_$PBS_ARRAY_INDEX.log
int main(int argc, char** argv)
{
  int hit=atoi(argv[1]);
  FILE *F1,*F2;
  int i,l,m,n,o,s,p,output;
  double ktmp;
  char OUTFILE[400],filename[400];
  Ntable.N_a=100;
  set_cov_parameters_to_(argv[2],1);// ini_files -> covparams
  //here: setting values internally

  // set this to zero to quickly run Gaussian-only covariances for testing
  if (covparams.ng==1){
    NG = 1;
  }
  else {
    NG = 0;
  }

  // set this to one to output details about inputs for diagnostics
  output = 0;
  FILE *F;
  printf("running multi_covariance_real with NG = %d\n",NG);
  
  set_cosmological_parameters_to_(argv[2],1);// ini_files -> cosmoparams?

  set_survey_parameters_to_(argv[2],1);// ini_files -> surveyparams?
  //init_clusters();
  //init_IA("none", "GAMA");
  //printf("test values: %d, %d, %s",redshift.clustering_photoz,tomo.clustering_Nbin,redshift.clustering_REDSHIFT_FILE);
  // printf("end of setup in main\n");

  double theta_min=covparams.tmin;
  double theta_max=covparams.tmax;
  int Ntheta=covparams.ntheta; 

  double logdt=(log(theta_max)-log(theta_min))/Ntheta;
  double *theta,*thetamin,*thetamax, *dtheta;
  theta=create_double_vector(0,Ntheta-1);
  thetamin=create_double_vector(0,Ntheta);
  thetamax=create_double_vector(0,Ntheta-1);
  dtheta=create_double_vector(0,Ntheta-1);
  for(i=0; i<Ntheta ; i++){
    thetamin[i]=exp(log(theta_min)+(i+0.0)*logdt);
    thetamax[i]=exp(log(theta_min)+(i+1.0)*logdt);
    theta[i] = 2./3.*(pow(thetamax[i],3.)-pow(thetamin[i],3.))/(pow(thetamax[i],2.)-pow(thetamin[i],2.));
    //printf("%d %e %e\n", i, thetamin[i],theta[i]);
    dtheta[i]=thetamax[i]-thetamin[i];
  }
  thetamin[Ntheta] = thetamax[Ntheta-1];
  like.theta=theta; // like.theta ponts to memory of theta
  like.Ntheta = Ntheta;

  double lmin=covparams.lmin;
  double lmax=covparams.lmax;
  int Nell=covparams.ncl; 
  double logdl=(log(lmax)-log(lmin))/Nell;
  double *ellmin, *dell;
  ellmin=create_double_vector(0,Nell);
  dell=create_double_vector(0,Nell-1);
  double ellmax;
  for(i=0; i<Nell ; i++){
    ellmin[i]=exp(log(lmin)+(i+0.0)*logdl);
    ellmax = exp(log(lmin)+(i+1.0)*logdl);
    dell[i]=ellmax-ellmin[i];
  }
  ellmin[Nell] = ellmax;
  like.Ncl = Nell;
  like.lmin = covparams.lmin;
  like.lmax = covparams.lmax;


  init_source_sample_();
  init_lens_sample_();

  init_probes("6x2pt");
  init_cmb("planck");

  /* pre-Calculate galaxy bias for src (5 bins) and lens (5 bins) galaxies */
  // double zbins[10] = {0.318457,0.518719,0.724785,0.993135,1.595836,0.320976,0.508596,0.686747,0.882423,1.131005};
  // double grow_z;
  // printf("galaxy bias: ");
  // for(i=0; i<10; i++){
  //   grow_z = growfac(1./(zbins[i]+1.))/growfac(1.0);
  //   printf("%lg, ", 1.05/grow_z);
  // }
  // printf("\ngalaxy bias computed!\n");
  // exit(0);

  int k=1;
  if (strcmp(covparams.ss,"true")==0)
  {
    sprintf(OUTFILE,"%s_ssss_++_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);

    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=l;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_shear_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,1,1,k);
        }

        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_ssss_--_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=l;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_shear_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,0,0,k);
        } 
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_ssss_+-_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_shear_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,1,0,k);
        } 
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ll,"true")==0)
  {
    sprintf(OUTFILE,"%s_llll_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){ 
      for (m=l;m<tomo.clustering_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_clustering_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,k); 
        }        
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ls,"true")==0)
  {
    sprintf(OUTFILE,"%s_lsls_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
      for (m=l;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ggl_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,k);  
        }        
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ls,"true")==0 && strcmp(covparams.ss,"true")==0)
  {
    sprintf(OUTFILE,"%s_lsss_+_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ggl_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,1,k);  
        } 
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_lsss_-_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ggl_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,0,k);  
        } 
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ll,"true")==0 && strcmp(covparams.ss,"true")==0)
  {
    sprintf(OUTFILE,"%s_llss_+_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_clustering_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,1,k);  
        } 
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_llss_-_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_clustering_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,0,k);  
        } 
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ll,"true")==0 && strcmp(covparams.ls,"true")==0)
  {
    sprintf(OUTFILE,"%s_llls_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){
      for (m=0;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_clustering_ggl_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.lk,"true")==0 && strcmp(covparams.ss,"true")==0)
  {
    sprintf(OUTFILE,"%s_lkss_+_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Nbin; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_gk_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,1,k);  
        } 
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_lkss_-_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Nbin; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_gk_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,0,k);  
        } 
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ks,"true")==0 && strcmp(covparams.ss,"true")==0)
  {
    sprintf(OUTFILE,"%s_ksss_+_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,1,k);  
        } 
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_ksss_-_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,0,k);  
        } 
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.lk,"true")==0 && strcmp(covparams.ls,"true")==0)
  {
    sprintf(OUTFILE,"%s_lkls_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Nbin; l++){
      for (m=0;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_gk_ggl_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ks,"true")==0 && strcmp(covparams.ls,"true")==0)
  {
    sprintf(OUTFILE,"%s_ksls_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=0;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_ggl_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.lk,"true")==0 && strcmp(covparams.ll,"true")==0)
  {
    sprintf(OUTFILE,"%s_lkll_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Nbin; l++){
      for (m=0;m<tomo.clustering_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_gk_clustering_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ks,"true")==0 && strcmp(covparams.ll,"true")==0)
  {
    sprintf(OUTFILE,"%s_ksll_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=0;m<tomo.clustering_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_clustering_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.lk,"true")==0)
  {
    sprintf(OUTFILE,"%s_lklk_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Nbin; l++){
      for (m=l;m<tomo.clustering_Nbin; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_gk_gk_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ks,"true")==0 && strcmp(covparams.lk,"true")==0)
  {
    sprintf(OUTFILE,"%s_kslk_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=0;m<tomo.clustering_Nbin; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_gk_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ks,"true")==0)
  {
    sprintf(OUTFILE,"%s_ksks_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=l;m<tomo.shear_Nbin; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_ks_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  ////////

  // if (strcmp(covparams.kk,"true")==0 && strcmp(covparams.ss,"true")==0)
  // {
  //   sprintf(OUTFILE,"%s_kkss_+_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
  //   for (m=0;m<tomo.shear_Npowerspectra; m++){
  //     if(k==hit) {
  //       sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
  //       // if (fopen(filename, "r") != NULL){exit(1);}
  //       run_cov_kk_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,m,1,k);  
  //     } 
  //     k=k+1;
  //   }
  //   sprintf(OUTFILE,"%s_kkss_-_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
  //   for (m=0;m<tomo.shear_Npowerspectra; m++){
  //     if(k==hit) {
  //       sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
  //       // if (fopen(filename, "r") != NULL){exit(1);}
  //       run_cov_kk_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,m,0,k);  
  //     } 
  //     k=k+1;
  //   }
  // }
  // if (strcmp(covparams.kk,"true")==0 && strcmp(covparams.ls,"true")==0)
  // {
  //   sprintf(OUTFILE,"%s_kkls_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
  //   for (m=0;m<tomo.ggl_Npowerspectra; m++){
  //     if(k==hit) {
  //       sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
  //       // if (fopen(filename, "r") != NULL){exit(1);}
  //       run_cov_kk_ggl_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,m,k);  
  //     }
  //     k=k+1;
  //   }
  // }
  // if (strcmp(covparams.kk,"true")==0 && strcmp(covparams.ll,"true")==0)
  // {
  //   sprintf(OUTFILE,"%s_kkll_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
  //   for (m=0;m<tomo.clustering_Npowerspectra; m++){
  //     if(k==hit) {
  //       sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
  //       // if (fopen(filename, "r") != NULL){exit(1);}
  //       run_cov_kk_clustering_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,m,k);  
  //     }
  //     k=k+1;
  //   }
  // }
  // if (strcmp(covparams.kk,"true")==0 && strcmp(covparams.lk,"true")==0)
  // {
  //   sprintf(OUTFILE,"%s_kklk_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
  //   for (m=0;m<tomo.clustering_Nbin; m++){
  //     if(k==hit) {
  //       sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
  //       // if (fopen(filename, "r") != NULL){exit(1);}
  //       run_cov_kk_gk_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,m,k);  
  //     }
  //     k=k+1;
  //   }
  // }
  // if (strcmp(covparams.kk,"true")==0 && strcmp(covparams.ks,"true")==0)
  // {
  //   sprintf(OUTFILE,"%s_kkks_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
  //   for (m=0;m<tomo.shear_Nbin; m++){
  //     if(k==hit) {
  //       sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
  //       // if (fopen(filename, "r") != NULL){exit(1);}
  //       run_cov_kk_ks_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,m,k);  
  //     }
  //     k=k+1;
  //   }
  // }
  // // real kkkk
  // if (strcmp(covparams.kk,"true")==0)
  // {
  //   sprintf(OUTFILE,"%s_kkkk_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
  //   if(k==hit) {
  //     sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
  //     // if (fopen(filename, "r") != NULL){exit(1);}
  //     run_cov_kk_kk_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,k);  
  //   }
  //   k=k+1;
  // }



  // Mixed covs
  if (strcmp(covparams.kk,"true")==0 && strcmp(covparams.ss,"true")==0)
  {
    sprintf(OUTFILE,"%s_kkss_+_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (m=0;m<tomo.shear_Npowerspectra; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_shear_mix_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta, ellmin,Nell,m,1,k);  
      } 
      k=k+1;
    }
    sprintf(OUTFILE,"%s_kkss_-_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (m=0;m<tomo.shear_Npowerspectra; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_shear_mix_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta, ellmin,Nell,m,0,k);  
      } 
      k=k+1;
    }
  }
  if (strcmp(covparams.kk,"true")==0 && strcmp(covparams.ls,"true")==0)
  {
    sprintf(OUTFILE,"%s_kkls_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (m=0;m<tomo.ggl_Npowerspectra; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_ggl_mix_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta, ellmin,Nell,m,k);  
      }
      k=k+1;
    }
  }
  if (strcmp(covparams.kk,"true")==0 && strcmp(covparams.ll,"true")==0)
  {
    sprintf(OUTFILE,"%s_kkll_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (m=0;m<tomo.clustering_Npowerspectra; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_clustering_mix_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta, ellmin,Nell,m,k);  
      }
      k=k+1;
    }
  }
  if (strcmp(covparams.kk,"true")==0 && strcmp(covparams.lk,"true")==0)
  {
    sprintf(OUTFILE,"%s_kklk_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (m=0;m<tomo.clustering_Nbin; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_gk_mix_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta, ellmin,Nell,m,k);  
      }
      k=k+1;
    }
  }
  if (strcmp(covparams.kk,"true")==0 && strcmp(covparams.ks,"true")==0)
  {
    sprintf(OUTFILE,"%s_kkks_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (m=0;m<tomo.shear_Nbin; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_ks_mix_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta, ellmin,Nell,m,k);  
      }
      k=k+1;
    }
  }
  // fourier kkkk
  // This does not include CMB smoothing.
  if (strcmp(covparams.kk,"true")==0)
  {
    sprintf(OUTFILE,"%s_kkkk_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    if(k==hit) {
      sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
      // if (fopen(filename, "r") != NULL){exit(1);}
      run_cov_kk_kk(OUTFILE,covparams.outdir,ellmin, dell,k);  
    }
    k=k+1;
  }

  if (hit==0)
  {
    sprintf(OUTFILE,"%s",covparams.filename);
    write_gglensing_zbins(OUTFILE);

    sprintf(OUTFILE,"%s%s.blocks",covparams.outdir,covparams.filename);
    F1 = fopen(OUTFILE,"w");
    fprintf(F1,"%d\n",k-1);
    fclose(F1);
  }

  printf("number of cov blocks for parallelization: %d\n",k-1); 
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0;   
}





