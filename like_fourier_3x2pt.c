#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
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
#include "../cosmolike_core/cfftlog/cfftlog.h"
#include "../cosmolike_core/cfftlog/utils.h"
#include "../cosmolike_core/cfastpt/cfastpt.h"
#include "../cosmolike_core/cfastpt/utils.h"

#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/emu17/P_cb/emu.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift_spline.c"
//#include "../cosmolike_core/theory/halo.c"
#include "../cosmolike_core/theory/halo_fast.c"
#include "../cosmolike_core/theory/HOD.c"
//#include "../cosmolike_core/theory/pt.c"
#include "../cosmolike_core/theory/pt_cfastpt.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/reduced_shear.c"
#include "../cosmolike_core/theory/CMBxLSS_fourier.c"
#include "../cosmolike_core/theory/cosmo2D_exact.c"
#include "../cosmolike_core/theory/cosmo2D_real.c"
#include "../cosmolike_core/theory/cosmo2D_fullsky.c"
#include "../cosmolike_core/theory/cosmo2D_fullsky_TATT.c"
#include "../cosmolike_core/theory/cluster.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "../cosmolike_core/theory/init_baryon.c"
#include "init_LSSxCMB.c"
//#include "../cosmolike_core/theory/init.c"

#include "../cosmolike_core/theory/priors_mpp.c"
// Naming convention:
// g = galaxy positions ("g" as in "galaxy")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "shear")
// And alphabetical order

typedef double (*C_tomo_pointer)(double l, int n1, int n2);
void twopoint_via_hankel(double **xi, double *logthetamin, double *logthetamax, C_tomo_pointer C_tomo, int ni, int nj, int N_Bessel);

#include "../cosmolike_core/theory/CMBxLSS_real_fullsky.c"

typedef struct input_nuisance_params_y3 {
    double bias[10];
    // double bias2[10];
    double b_mag[10];
    double lens_z_bias[10];
    double source_z_bias[10];
    double shear_m[10];
    double p_ia[10];
    double Q1;
    double Q2;
    double Q3;
} input_nuisance_params_y3;

typedef struct input_cosmo_params_y3 {
    double omega_m;
    double sigma_8;
    double A_s;
    double n_s;
    double w0;
    double wa;
    double omega_b;
    double omega_nuh2;
    double h0;
    double MGSigma;
    double MGmu;
    double theta_s;
} input_cosmo_params_y3;

double C_shear_tomo_sys(double *ell, int nl, int z1, int z2);
double C_gl_tomo_sys(double *ell, int nl, int zl, int zs);
double C_cl_tomo_sys(double *ell, int nl, int zs);

void set_data_shear(double *ell, double *data, int start);
void set_data_ggl(double *ell, double *data, int start);
void set_data_clustering(double *ell, double *data, int start);
void set_data_gk(double *ell, double *data, int start);
void set_data_ks(double *ell, double *data, int start);
void set_data_kk_logbin(double *ell, double *data, int start);
void set_data_kk_bandpower(double *data, int start);
void compute_data_vector(char *filename, double OMM, double S8, double NS, double W0,double WA, double OMB, double OMNUh2, double H0, double MGSigma, double MGmu, double THETA_S, \
                      double *B, double *b_mag,\
                      double *SP, double *CP, double *M, \
                      double *p_ia);
double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double OMNUh2, double H0, double MGSigma, double MGmu, double THETA_S, \
                      double *B, double *b_mag,\
                      double *SP, double *CP, double *M, \
                      double *p_ia, double Q1, double Q2, double Q3);
void write_datavector_wrapper(char *filename, input_cosmo_params_y3 ic, input_nuisance_params_y3 in);
double log_like_wrapper(input_cosmo_params_y3 ic, input_nuisance_params_y3 in);

double get_sigma_8(input_cosmo_params_y3 ic);
double get_h0(input_cosmo_params_y3 ic);

double get_h0(input_cosmo_params_y3 ic){
  return cosmology.h0;
}

double get_sigma_8(input_cosmo_params_y3 ic){
  if (ic.A_s != cosmology.A_s){
    printf("cosmology changed before calling get_sigma_8\n");
    return -1.;
  } 
  return cosmology.sigma_8;
}


double C_shear_tomo_sys(double *ell, int nl, int z1, int z2)
{
  double Cl = 0.0;
  int l_min, l_max, l, Nl;
  l_min = (int)ceil(ell[nl]);
  l_max = (int)ceil(ell[nl+1])-1;
  Nl = l_max - l_min + 1;
  for(l=l_min; l<=l_max; l++){
    if(like.IA==0){
      Cl += C_shear_tomo((double)l, z1, z2);
    }
    else if (like.IA==1||like.IA==3||like.IA==4){
      Cl += C_shear_shear_IA((double)l, z1, z2);
    }
    else if (like.IA==5||like.IA==6){
      Cl += C_EE_TATT((double)l, z1, z2);
    }
    else{
      printf("like.IA = %d is not supported!\n", like.IA);
      exit(-1);
    }
  }
  Cl /= (double)Nl;
  if(like.shearcalib==1) Cl *=(1.0+nuisance.shear_calibration_m[z1])*(1.0+nuisance.shear_calibration_m[z2]);
  return Cl;
}

double C_gl_tomo_sys(double *ell, int nl, int zl, int zs)
{
  double Cl = 0.0;
  int l_min, l_max, l, Nl;
  l_min = (int)ceil(ell[nl]);
  l_max = (int)ceil(ell[nl+1])-1;
  Nl = l_max - l_min + 1;
  for(l=l_min; l<=l_max; l++){
    if(like.IA==0){
      Cl += C_gl_tomo((double)l, zl, zs);
    }
    else if (like.IA==1||like.IA==3||like.IA==4){
      // w/ RSD and magnification bias
      Cl += C_ggl_IA((double)l, zl, zs);
    }
    else if (like.IA==5||like.IA==6){
      // w/ magnification bias, RSD?
      Cl += C_ggl_TATT((double)l, zl, zs);
    }
    else{
      printf("like.IA = %d is not supported!\n", like.IA);
      exit(-1);
    }
  }
  Cl /= (double)Nl;
  if(like.shearcalib==1) Cl *=(1.0+nuisance.shear_calibration_m[zs]);
  return Cl;
}

double C_cl_tomo_sys(double *ell, int nl, int zs)
{
  double Cl = 0.0;
  int l_min, l_max, l, Nl;
  l_min = (int)ceil(ell[nl]);
  l_max = (int)ceil(ell[nl+1])-1;
  Nl = l_max - l_min + 1;
  for(l=l_min; l<=l_max; l++){
    // w/ magnification bias and RSD 
    Cl += C_cl_tomo((double)l, zs, zs);
  }
  Cl /= (double)Nl;
  return Cl;
}

void set_data_shear(double *ell, double *data, int start)
{
  int i,z1,z2,nz,j;
  for (nz = 0; nz < tomo.shear_Npowerspectra; nz++){
    z1 = Z1(nz); z2 = Z2(nz);
    for (i = 0; i < like.Ncl; i++){
      if (mask(start+like.Ncl*nz+i)){
        data[start+like.Ncl*nz+i] = C_shear_tomo_sys(ell, i, z1, z2);
      }
    }
  }
}

void set_data_ggl(double *ell, double *data, int start)
{
  int i, zl,zs,nz;  
  for (nz = 0; nz < tomo.ggl_Npowerspectra; nz++){
    zl = ZL(nz); zs = ZS(nz);
    //printf("ggl bin combos %d %d\n",zl,zs);
    for (i = 0; i < like.Ncl; i++){
      if (mask(start+(like.Ncl*nz)+i)){
        data[start+(like.Ncl*nz)+i] = C_gl_tomo_sys(ell, i, zl, zs);
      }
    }
  }
}

void set_data_clustering(double *ell, double *data, int start)
{
  int i,nz,j;
  for (nz = 0; nz < tomo.clustering_Npowerspectra; nz++){
    for (i = 0; i < like.Ncl; i++){
      if (mask(start+(like.Ncl*nz)+i)){
        data[start+(like.Ncl*nz)+i] = C_cl_tomo_sys(ell, i, nz);
      }
    }
  }
}

/*void set_data_gk(double *theta, double *data, int start)
{
   for (int nz=0; nz<tomo.clustering_Nbin; nz++){
      for (int i=0; i<like.Ncl; i++){
         if (mask(start+(like.Ncl*nz)+i)){
            data[start+(like.Ncl*nz)+i] = w_gk_fullsky(i,nz);
         }
         else{
            data[start+(like.Ncl*nz)+i] = 0.;
         }
      }
   }
}

void set_data_ks(double *theta, double *data, int start)
{
   for (int nz=0; nz<tomo.shear_Nbin; nz++){
      for (int i=0; i<like.Ncl; i++){
         if (mask(start+(like.Ncl*nz)+i)){
            data[start+(like.Ncl*nz)+i] = w_ks_sys(i,nz);
         }
         else{
            data[start+(like.Ncl*nz)+i] = 0.;
         }
      }
   }
}

void set_data_kk_logbin(double *ell, double *data, int start)
{
   for (int i=0; i<like.Ncl; i++){
      if (mask(start+i)){
         data[start+i] = C_kk_nointerp(ell[i]);
      }
      else{
         data[start+i] = 0.;
      }
   }
}

void set_data_kk_bandpower(double *data, int start)
{
  for(int L=like.lmin_bp_with_corr; L<like.lmax_bp_with_corr+1; L++)
  {
    //double _C_kk_limber = C_kk((double)L);// use the interpolate version
    double _C_kk_limber = C_kk_nointerp((double)L);
    for(int i=0; i<like.Nbp; i++)
    {
      if(mask(start+i))
      {
        data[start+i] += _C_kk_limber * binmat_read(1, i, L-like.lmin_bp_with_corr);
      }
    }
  }
  for(int i=0; i<like.Nbp; i++)
  {
    if(mask(start+i))
    {
      data[start+i] -= ckk_offset_read(1, i);
    }
  }
}*/

int set_cosmology_params(double OMM, double NORM, double NS, double W0,double WA, double OMB, double OMNUh2, double H0, double MGSigma, double MGmu, double THETA_S)
{
  cosmology.Omega_m=OMM;
  cosmology.Omega_v= 1.0-cosmology.Omega_m;
  if (NORM < 1.e-7){
    cosmology.A_s = NORM;
    cosmology.sigma_8 = 0.;
  }
  else{
    cosmology.sigma_8=NORM;
    cosmology.A_s = 0.;
  }
  cosmology.theta_s = THETA_S;
  cosmology.n_spec= NS;
  cosmology.w0=W0;
  cosmology.wa=WA;
  cosmology.omb=OMB;
  if (H0> 0){
    cosmology.Omega_nu=OMNUh2/H0/H0;
  }
  else{cosmology.Omega_nu =0.0;}
  cosmology.h0=H0;
  cosmology.MGSigma =  MGSigma;
  cosmology.MGmu =  MGmu;

   if (cosmology.Omega_m < 0.05 || cosmology.Omega_m > 0.6) return 0;
   if (cosmology.omb < 0.04 || cosmology.omb > 0.055) return 0;
   //if (cosmology.sigma_8 < 0.5 || cosmology.sigma_8 > 1.1) return 0;
   if (cosmology.n_spec < 0.84 || cosmology.n_spec > 1.06) return 0;
   if (cosmology.w0 < -2.1 || cosmology.w0 > -0.0) return 0;
   if (cosmology.wa < -2.6 || cosmology.wa > 2.6) return 0;
   if (cosmology.h0 < 0.4 || cosmology.h0 > 0.9) return 0;
  
  /*printf("cosmology.theta_s = %le \n", cosmology.theta_s);
  printf("cosmology.A_s = %le \n", cosmology.A_s);
  printf("cosmology.w0 = %le \n", cosmology.w0);
  printf("cosmology.wa = %le \n", cosmology.wa);
  printf("cosmology.sigma_8 = %le \n", cosmology.sigma_8);
  printf("cosmology.h0= %le \n", cosmology.h0);
  printf("cosmology.omb = %le \n", cosmology.omb);
  printf("cosmology.Omega_m = %le \n", cosmology.Omega_m);
  printf("cosmology.Omega_nu = %le \n", cosmology.Omega_nu);*/

  return 1;
}

void set_nuisance_shear_calib(double *M)
{
  int i;
  for(i=0;i<tomo.shear_Nbin;i++) {nuisance.shear_calibration_m[i] = M[i];}
}

int set_nuisance_shear_photoz(double *SP)
{
  int i;
  for(i=0;i<tomo.shear_Nbin;i++) {nuisance.bias_zphot_shear[i]=SP[i];}
  
  // for (i=0;i<tomo.shear_Nbin; i++){ 
  //   nuisance.sigma_zphot_shear[i]=SPS1;
  //   if (nuisance.sigma_zphot_shear[i]<0.0001) return 0;
  // }
  return 1;
}

int set_nuisance_clustering_photoz(double *CP)
{
  int i;
  for(i=0;i<tomo.clustering_Nbin;i++) {nuisance.bias_zphot_clustering[i]=CP[i];}
  
  // for (i=0;i<tomo.clustering_Nbin; i++){ 
  //   nuisance.sigma_zphot_clustering[i]=CPS1;
  //   if (nuisance.sigma_zphot_clustering[i]<0.0001) return 0;
  // }
  return 1;
}

int set_nuisance_ia(double *p_ia)
{
  switch (like.IA)
  {
    case 4: // NLA power-law
      nuisance.A_ia=p_ia[0];
      nuisance.eta_ia=p_ia[1];
      nuisance.oneplusz0_ia = 1.62;
      break;
    case 6: // TATT power-law
      nuisance.A_ia = p_ia[0];
      nuisance.eta_ia = p_ia[1];
      nuisance.b_ta_z[0] = p_ia[2];
      nuisance.A2_ia = p_ia[3];
      nuisance.eta_ia_tt = p_ia[4];
      nuisance.oneplusz0_ia = 1.62;
      break;
    default:
      printf("ERROR: like.IA can not be %d! Only support 4,6\n", like.IA);
      exit(-1);
  }
  // if (nuisance.A_ia < 0.0 || nuisance.A_ia > 10.0) return 0;
  // if (nuisance.eta_ia < -10.0 || nuisance.eta_ia> 10.0) return 0;
  return 1;
}

int set_nuisance_gbias(double *B)
{

  int i;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    gbias.b[i] = B[i];
  }
  return 1;
} 

int set_nuisance_bmag(double *b_mag)
{

  int i;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    gbias.b_mag[i] = b_mag[i];
  }
  return 1;
} 

double log_multi_like(double OMM, double NORM, double NS, double W0,double WA, double OMB, double OMNUh2, double H0, double MGSigma, double MGmu, double THETA_S, \
                      double *B, double *b_mag,\
                      double *SP, double *CP, double *M, \
                      double *p_ia, double Q1, double Q2, double Q3)
{
  int i,j,k,m=0,l;
  // printf("%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, \n", OMM,NORM,NS,W0,WA,OMB,OMNUh2,H0, MGSigma, MGmu, THETA_S);
  // for(i=0;i<10;i++){
  //   printf("%lg, ", B[i]);
  // }printf("\n");
  // for(i=0;i<10;i++){
  //   printf("%lg, ", b_mag[i]);
  // }printf("\n");
  // for(i=0;i<10;i++){
  //   printf("%lg, ", SP[i]);
  // }printf("\n");
  // for(i=0;i<10;i++){
  //   printf("%lg, ", CP[i]);
  // }printf("\n");
  // for(i=0;i<10;i++){
  //   printf("%lg, ", M[i]);
  // }printf("\n");
  // for(i=0;i<10;i++){
  //   printf("%lg, ", p_ia[i]);
  // }printf("\n");
  // printf("finish print parameters\n");

  static double *pred;
  static double *ellmin;
  static double darg;
  double chisqr,a,log_L_prior=0.0, log_L=0.0;;
  
  pred= create_double_vector(0, like.Ndata-1);

  ellmin = create_double_vector(0, like.Ncl);
  darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
  for (l=0; l<like.Ncl; l++){
    ellmin[l]=exp(log(like.lmin)+(l+0.0)*darg);
  }
  ellmin[like.Ncl] = like.lmax;

  set_cosmology_params(OMM,NORM,NS,W0,WA,OMB,OMNUh2,H0, MGSigma, MGmu, THETA_S);
  if (strcmp(pdeltaparams.runmode,"class")==0||strcmp(pdeltaparams.runmode,"CLASS")==0) {
    int status = 0;
    if (H0> 0 &&(OMB*H0*H0 >= 0.04 || OMB*H0*H0 <= 0.005)){printf("BBN\n"); return -1.e+15;}
    p_class(1.,1.,0,&status);
    if (status){printf("CLASS error\n"); return -1.e+15;}
  }
  set_nuisance_shear_calib(M);
  if (set_nuisance_shear_photoz(SP)==0){
    printf("Shear photo-z sigma too small\n");
    return -1.0e15;
  }
  if (set_nuisance_clustering_photoz(CP)==0){
    printf("Clustering photo-z sigma too small\n");
    return -1.0e15;
  }
  if (set_nuisance_ia(p_ia)==0){
    printf("IA parameters out of bounds\n");
    return -1.0e15; 
  }
  if (set_nuisance_gbias(B)==0){
    printf("Bias out of bounds\n");
    return -1.0e15;
  }
  if (set_nuisance_bmag(b_mag)==0){
    printf("b_mag out of bounds\n");
    return -1.0e15;
  }
  // printf("like %le %le %le %le %le %le %le %le\n",cosmology.Omega_m, cosmology.Omega_v,cosmology.sigma_8,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,cosmology.h0); 
  // printf("like %le %le %le %le\n",gbias.b[0][0], gbias.b[1][0], gbias.b[2][0], gbias.b[3][0]);    
  // for (i=0; i<10; i++){
  //   printf("nuisance %le %le %le\n",nuisance.shear_calibration_m[i],nuisance.bias_zphot_shear[i],nuisance.sigma_zphot_shear[i]);
  // }

  log_L_prior=0.0;
  // if(like.Aubourg_Planck_BAO_SN==1) log_L_prior+=log_L_Planck_BAO_SN();
  // if(like.SN==1) log_L_prior+=log_L_SN();
  //if(like.BAO==1) log_L_prior+=log_L_BAO();
  // if(like.Planck==1) log_L_prior+=log_L_Planck();
  // if(like.Planck15_BAO_w0wa==1) log_L_prior+=log_L_Planck15_BAO_w0wa();//CH
  //if(like.Planck15_BAO_H070p6_JLA_w0wa==1) log_L_prior+=log_L_Planck15_BAO_H070p6_JLA_w0wa();//CH
  // Flat priors are implemented in python wrapper
  // if(like.IA!=0) log_L_prior+=log_L_ia();
  // if(like.IA!=0) log_L_prior+=log_like_f_red();
  //test likelihood evaluation JX
  if(like.wlphotoz!=0) log_L_prior+=log_L_wlphotoz();
  if(like.clphotoz!=0) log_L_prior+=log_L_clphotoz();
  if(like.shearcalib==1) log_L_prior+=log_L_shear_calib();

  // Flat priors are implemented in python wrapper
  // if(like.IA!=0) {
  //   log_L = 0.0;
  //   log_L -= pow((nuisance.A_ia - prior.A_ia[0])/prior.A_ia[1],2.0);
  //   log_L -= pow((nuisance.eta_ia - prior.eta_ia[0])/prior.eta_ia[1],2.0);
  //   log_L_prior+=0.5*log_L;
  // }
  // if(like.baryons==1){
  //   log_L = 0.0;
  //   log_L -= pow((Q1 - prior.bary_Q1[0])/prior.bary_Q1[1],2.0);
  //   log_L -= pow((Q2 - prior.bary_Q2[0])/prior.bary_Q2[1],2.0);
  //   log_L -= pow((Q3 - prior.bary_Q3[0])/prior.bary_Q3[1],2.0);
  //   log_L_prior+=0.5*log_L;
  // }
 
  // if(like.clusterMobs==1) log_L_prior+=log_L_clusterMobs();
 
  // printf("%d %d %d %d\n",like.BAO,like.wlphotoz,like.clphotoz,like.shearcalib);
  // printf("logl %le %le %le %le\n",log_L_shear_calib(),log_L_wlphotoz(),log_L_clphotoz(),log_L_clusterMobs());
  int start=0;  

  if(like.shear_shear==1) {
    set_data_shear(ellmin, pred, start);
    start=start+like.Ncl*tomo.shear_Npowerspectra; 
  }
  if(like.shear_pos==1){
    set_data_ggl(ellmin, pred, start);
    start=start+like.Ncl*tomo.ggl_Npowerspectra;
  } 
  if(like.pos_pos==1){
    set_data_clustering(ellmin, pred, start);
    start=start+like.Ncl*tomo.clustering_Npowerspectra;
  }
  /*if(like.gk==1) {
    set_data_gk(theta, pred, start);
    start += like.Ntheta*tomo.clustering_Nbin;
  }
  if(like.ks==1) {
    set_data_ks(theta, pred, start);
    start += like.Ntheta*tomo.shear_Nbin;
  } 
  if(like.kk==1) {
    if(like.Ncl>0){
      set_data_kk_logbin(ell, pred, start);
      start += like.Ncl;
    }
    else{
      set_data_kk_bandpower(pred, start);
      start += like.Nbp;
    }
  }*/
  
  chisqr=0.0;
  for (i=0; i<like.Ndata; i++){
    for (j=0; j<like.Ndata; j++){
      //a=(pred[i]-data_read(1,i)+Q1*bary_read(1,0,i)+Q2*bary_read(1,1,i)+Q3*bary_read(1,2,i))*invcov_mask(1,i,j)*(pred[j]-data_read(1,j)+Q1*bary_read(1,0,j)+Q2*bary_read(1,1,j)+Q3*bary_read(1,2,j));
      a=(pred[i]-data_read(1,i))*invcov_mask(1,i,j)*(pred[j]-data_read(1,j));
      chisqr=chisqr+a;
    }
    // if (fabs(data_read(1,i)) < 1.e-25){
    //    printf("%d %le %le %le\n",i,data_read(1,i),pred[i],invcov_read(1,i,i));
    // }
  }
  if (chisqr<0.0){
    printf("error: chisqr = %le\n",chisqr);
    //exit(EXIT_FAILURE);
  }
  if (chisqr<-1.0) exit(EXIT_FAILURE);
  if (isnan(chisqr)){return -1.e+16;}
  //printf("%le\n",chisqr);
  return -0.5*chisqr+log_L_prior;
}

void compute_data_vector(char *filename, double OMM, double NORM, double NS, double W0,double WA, double OMB, double OMNUh2, double H0, double MGSigma, double MGmu, double THETA_S, \
                      double *B, double *b_mag,\
                      double *SP, double *CP, double *M, \
                      double *p_ia){

  int i,j,k,m=0,l;
  static double *pred;
  static double *ellmin;
  static double darg;
  double chisqr,a,log_L_prior=0.0;
  
  pred= create_double_vector(0, like.Ndata-1);

  ellmin = create_double_vector(0, like.Ncl);
  darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
  for (l=0; l<like.Ncl; l++){
    ellmin[l]=exp(log(like.lmin)+(l+0.0)*darg);
  }
  ellmin[like.Ncl] = like.lmax;

  if(set_cosmology_params(OMM,NORM,NS,W0,WA,OMB,OMNUh2,H0, MGSigma, MGmu, THETA_S)){printf("set cosmo params success\n");}
  else{printf("set cosmo params failed\n");exit(1);}
  set_nuisance_shear_calib(M);
  set_nuisance_shear_photoz(SP);
  set_nuisance_clustering_photoz(CP);
  set_nuisance_ia(p_ia);
  set_nuisance_gbias(B);
  set_nuisance_bmag(b_mag);
  printf("Setting model vector\n");
  int start=0;  
  if(like.shear_shear==1) {
    printf("Start with shear-shear\n");
    set_data_shear(ellmin, pred, start);
    start=start+like.Ncl*tomo.shear_Npowerspectra;
    printf("Done with shear-shear\n");
  }
  if(like.shear_pos==1){
    printf("Start with galaxy-galaxy lensing\n");
    set_data_ggl(ellmin, pred, start);
    start=start+like.Ncl*tomo.ggl_Npowerspectra;
    printf("Done with galaxy-galaxy lensing\n");
  } 
  if(like.pos_pos==1){
    printf("Start with clustering\n");
    set_data_clustering(ellmin, pred, start);
    start=start+like.Ncl*tomo.clustering_Npowerspectra;
    printf("Done with clustering\n");
  }
  /*if(like.gk==1) {
    printf("Start with galaxy-kappa\n");
    set_data_gk(theta, pred, start);
    start += like.Ntheta*tomo.clustering_Nbin;
    printf("Done with galaxy-kappa\n");
  }
  if(like.ks==1) {
    printf("Start with shear-kappa\n");
    set_data_ks(theta, pred, start);
    start += like.Ntheta*tomo.shear_Nbin;
    printf("Done with shear-kappa\n");
  } 
  if(like.kk==1) {
    if(like.Ncl > 0)
    {
      printf("Start with kappa-kappa (log-bin)\n");
      set_data_kk_logbin(ell, pred, start);
      start += like.Ncl;
      printf("Done with kappa-kappa (log-bin)\n");
    }
    else
    {
      printf("Start with kappa-kappa (band-power)\n");
      set_data_kk_bandpower(pred, start);
      start += like.Nbp;
      printf("Done with kappa-kappa (band-power)\n");
    }
  }*/

  FILE *F;
  F=fopen(filename,"w");
  if(F==NULL){
    printf("ERORR: Can not open file %s\nAborting...\n", filename);
    exit(1);
  }
  for (i=0;i<like.Ndata; i++){  
    fprintf(F,"%d %le\n",i,pred[i]);
    //printf("%d %le\n",i,pred[i]);
  }
  fclose(F);
  // printf("&gbias.b1_function %p\n",&gbias.b1_function);
  // printf("gbias.b1_function  %p\n",gbias.b1_function);
  // printf("bgal_z   %p\n",bgal_z);
  // printf("&bgal_z  %p\n",&bgal_z);
  // printf("b1_per_bin   %p\n",b1_per_bin);
  // printf("&b1_per_bin  %p\n",&b1_per_bin);

}


void write_datavector_wrapper(char *filename, input_cosmo_params_y3 ic, input_nuisance_params_y3 in)
{
  printf("write_datavector_wrapper: path to test data vector: %s\n",filename);

  double NORM;
  if (ic.A_s > 0. && ic.A_s < 1.e-5){NORM = ic.A_s;}
  else{NORM = ic.sigma_8;}
  if (NORM <= 0){
    printf("write_datavector_wrapper called with A_s = %e, sigma_8 =%e\nEXIT\n",ic.A_s,ic.sigma_8);
    exit(1);
  }
  compute_data_vector(filename, 
    ic.omega_m, NORM, ic.n_s, ic.w0, ic.wa, ic.omega_b, ic.omega_nuh2, 
    ic.h0, ic.MGSigma, ic.MGmu, ic.theta_s,
    in.bias, in.b_mag,
    in.source_z_bias, in.lens_z_bias, in.shear_m, 
    in.p_ia);
}

double log_like_wrapper(input_cosmo_params_y3 ic, input_nuisance_params_y3 in)
{
  double NORM;
  if (ic.A_s > 0. && ic.A_s < 1.e-5){NORM = ic.A_s;}
  else{NORM = ic.sigma_8;}
  if (NORM <= 0){
    printf("log_like_wrapper called with A_s = %e, sigma_8 =%e\nEXIT\n",ic.A_s,ic.sigma_8);
    exit(1);
  }
  double like = log_multi_like(ic.omega_m, ic.sigma_8, ic.n_s, ic.w0, ic.wa, ic.omega_b,ic.omega_nuh2, ic.h0, ic.MGSigma, ic.MGmu,ic.theta_s, 
    in.bias, in.b_mag,
    in.source_z_bias,in.lens_z_bias,in.shear_m, 
    in.p_ia, in.Q1, in.Q2, in.Q3);
  return like;
}



void save_zdistr_sources(int zs){
  double z,dz =(redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/300.0;
  printf("Printing redshift distribution n(z) for source redshift bin %d\n",zs);
  
   FILE *F1;
   char filename[300];
   sprintf(filename,"zdistris/zdist_sources_bin%d.txt",zs);
   F1 = fopen(filename,"w");
   for (z =redshift.shear_zdistrpar_zmin; z< redshift.shear_zdistrpar_zmax; z+= dz){
      fprintf(F1,"%e %e\n", z, zdistr_photoz(z,zs));
   }
}


void save_zdistr_lenses(int zl){
   double z,dz =(redshift.clustering_zdistrpar_zmax-redshift.clustering_zdistrpar_zmin)/300.0;
  printf("Printing redshift distribution n(z) and bias b(z) for lens redshift bin %d\n",zl);
   
   FILE *F1;
   char filename[300];
   sprintf(filename,"zdistris/zdist_lenses_bin%d.txt", zl);
   F1 = fopen(filename,"w");
   for (z =redshift.clustering_zdistrpar_zmin; z< redshift.clustering_zdistrpar_zmax; z+= dz){
      fprintf(F1,"%e %e\n", z, pf_photoz(z,zl));
   }
}



