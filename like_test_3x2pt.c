#include "like_fourier_6x2pt.c"

// void test_Cls_desy3_planck()
// {
//   clock_t begin, end;
//   double time_spent;

//   int i;

//   init_binning_fourier(20, 30., 3000.);
//   init_scalecuts(20., 3000.); // Rmin_bias = 20Mpc/h, lmax_shear=3000

//   sprintf(survey.name,"%s","DESY3xplanck");

//   // init_bary(argv[2]);

//   // init_priors(0.002,sigma_zphot_shear[sce],0.001,0.001,sigma_zphot_clustering[sce],0.001,0.001,3.0,1.2,3.8,2.0,16.0,5.0,0.8);


//   double NORM;
//   like.IA = 4; 
// /* sample parameter values */
//   // nuisance.oneplusz0_ia = 1.62;

// /* here, do your time-consuming job */
//   // init_cosmo_runmode("CLASS");
//   init_cosmo_runmode("halofit");
// //  init_bary("owls_AGN");
//   cosmology.Omega_m   = 0.3;
//   NORM    = 0.82355;
//   // NORM  = 2.19e-09;
//   cosmology.n_spec    = 0.97;
  
//   cosmology.w0=-1.;
//   cosmology.wa=0.;
//   cosmology.omb=0.048;
//   cosmology.h0=0.69;
//   double Omega_nuh2 = 0.0;//0.00083;//0.001743331232934258;
//   cosmology.Omega_v=1.0-cosmology.Omega_m;
// //  cosmology.theta_s  = -1.;
//   double b1[10]={0.,0.,0.,0.,0.,2.,2.,2.,2.,2.}, b2[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},b_mag[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//   b1[0] = 1.38927;
//   b1[1] = 1.68849;
//   b1[2] = 1.38347;
//   b1[3] = 1.51844;
//   b1[4] = 1.59457;
// //
// //  b2[0] = 0.23; b2[1] = 0.23; b2[2] = 0.23; b2[3] = 0.5; b2[4] = .5;
//   b_mag[0] = 0.0347082;
//   b_mag[1] = 0.0230357;
//   b_mag[2] = 0.0235965;
//   b_mag[3] = 0.0507362;
//   b_mag[4] = 0.0961916;

//   double A_ia=0.0591106, eta_ia=0.0641435;
//   double p_ia[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//   p_ia[0] = A_ia; p_ia[1] = eta_ia;

//   double mean_m[10]={-0.00289242 , 0.00135348,0.000818622,-0.00202546 ,0.0,0.0,0.0,0.0,0.0,0.0};
//   double sigma_m[10]={0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03};
//   double bias_photoz_s[10]={-0.0027565,-1.09488e-05,-0.00255251,-5.42763e-05,0.0,0.0,0.0,0.0,0.0,0.0};
//   // double sigma_b_photoz_s[10]={0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08};
//   double bias_photoz_l[10]={0.00187054,0.00241709,0.000165675,0.00325307,-0.000617269,0.0,0.0,0.0,0.0,0.0};
//   // double sigma_b_photoz_l[10]={0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};

//   // init_source_sample_mpp("./zdistris/nz_v0.16_smoothed.txt",4);
//   // init_lens_sample_mpp("./zdistris/nz_y3_redmagic_v0.5.1_wide_gold_2.2.1_combined_hd3_hl2_z_samp.txt",5,b1,b2,0.0);

//   init_source_sample_mpp("./zdistris/mcal_1101_source.nz",4);
//   init_lens_sample_mpp("./zdistris/mcal_1101_lens.nz",5,b1,b2,0.0);


//   // init_binning_mpp(20,2.5,250.);

//   init_probes("6x2pt");
//   init_cmb("planck");
  
//   // set_shear_priors_mpp(mean_m,sigma_m);
//   //set_wlphotoz_priors_mpp(bias_photoz_s,sigma_b_photoz_s);
//   //set_clphotoz_priors_mpp(bias_photoz_l,sigma_b_photoz_l);

//   sprintf(like.MASK_FILE,"%s","none");
//   printf("PATH TO MASK: %s\n",like.MASK_FILE);
//   begin = clock();
//   char datavfile[200];
//   sprintf(datavfile,"datav/xi_Y3_test_test");
// //  sprintf(datavfile,"datav/xi_Y3_baseline+b2_MICE+bary_owls_AGN");
//   compute_data_vector(datavfile,0.325498, 0.768139, 0.979752, -1, 0.0, 0.0472401, 0.0, 0.691452, 0.0, 0.0, 0.0,
//     b1, b_mag,
//     bias_photoz_s, //source photo-z bias
//     bias_photoz_l, //lens photo-z bias
//     mean_m, //shear calibration
//     p_ia); // IA
  
//   end = clock();
//   time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
//    printf("time spent %le\n",time_spent);
//    printf("like.IA = %d\n", like.IA);

// }

void test_Cls_desy3()
{
  clock_t begin, end;
  double time_spent;

  int i;

  init_binning_fourier(20, 30., 3000.);
  init_scalecuts(20., 3000.); // Rmin_bias = 20Mpc/h, lmax_shear=3000

  sprintf(survey.name,"%s","DESY3");

  // init_bary(argv[2]);

  // init_priors(0.002,sigma_zphot_shear[sce],0.001,0.001,sigma_zphot_clustering[sce],0.001,0.001,3.0,1.2,3.8,2.0,16.0,5.0,0.8);


  double NORM;
  like.IA = 4; 
/* sample parameter values */
  // nuisance.oneplusz0_ia = 1.62;

/* here, do your time-consuming job */
  // init_cosmo_runmode("CLASS");
  init_cosmo_runmode("halofit");
//  init_bary("owls_AGN");
  cosmology.Omega_m   = 0.3;
  NORM    = 0.82355;
  // NORM  = 2.19e-09;
  cosmology.n_spec    = 0.97;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.048;
  cosmology.h0=0.69;
  double Omega_nuh2 = 0.0;//0.00083;//0.001743331232934258;
  cosmology.Omega_v=1.0-cosmology.Omega_m;
//  cosmology.theta_s  = -1.;
  double b1[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, b2[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},b_mag[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  b1[0] = 1.7;
  b1[1] = 1.7;
  b1[2] = 1.7;
  b1[3] = 2.0;
  b1[4] = 2.0;
//
//  b2[0] = 0.23; b2[1] = 0.23; b2[2] = 0.23; b2[3] = 0.5; b2[4] = .5;
  b_mag[0] = -0.19375;
  b_mag[1] = -0.6285407;
  b_mag[2] = -0.69319886;
  b_mag[3] = 1.17735723;
  b_mag[4] = 1.87509758;

  double p_ia[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double A_ia=0.5, eta_ia=0.;
  p_ia[0] = A_ia; p_ia[1] = eta_ia;

  double mean_m[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double sigma_m[10]={0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03};
  double bias_photoz_s[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  // double sigma_b_photoz_s[10]={0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08};
  double bias_photoz_l[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  // double sigma_b_photoz_l[10]={0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};

  // init_source_sample_mpp("./zdistris/nz_v0.16_smoothed.txt",4);
  // init_lens_sample_mpp("./zdistris/nz_y3_redmagic_v0.5.1_wide_gold_2.2.1_combined_hd3_hl2_z_samp.txt",5,b1,b2,0.0);

  init_source_sample_mpp("./zdistris/mcal_1101_source.nz",4);
  init_lens_sample_mpp("./zdistris/mcal_1101_lens.nz",5,b1,b2,0.0);


  // init_binning_mpp(20,2.5,250.);

  init_probes("3x2pt");
  
  // set_shear_priors_mpp(mean_m,sigma_m);
  //set_wlphotoz_priors_mpp(bias_photoz_s,sigma_b_photoz_s);
  //set_clphotoz_priors_mpp(bias_photoz_l,sigma_b_photoz_l);

  sprintf(like.MASK_FILE,"%s","none");
  printf("PATH TO MASK: %s\n",like.MASK_FILE);
  begin = clock();
  char datavfile[200];
  sprintf(datavfile,"datav/xi_Y3_3x2pt");
//  sprintf(datavfile,"datav/xi_Y3_baseline+b2_MICE+bary_owls_AGN");
  compute_data_vector(datavfile,cosmology.Omega_m,NORM ,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,Omega_nuh2,cosmology.h0,0.0,0.0,cosmology.theta_s,
    b1, b_mag,
    bias_photoz_s, //source photo-z bias
    bias_photoz_l, //lens photo-z bias
    mean_m, //shear calibration
    p_ia); // IA
  
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
   printf("time spent %le\n",time_spent);
   printf("like.IA = %d\n", like.IA);

}

void test_Cls_desy6()
{
  clock_t begin, end;
  double time_spent;

  int i;

  init_binning_fourier(15, 40., 3000.);
  init_scalecuts(20., 3000.); // Rmin_bias = 20Mpc/h, lmax_shear=3000


  sprintf(survey.name,"%s","DESY6");

  // init_bary(argv[2]);

  // init_priors(0.002,sigma_zphot_shear[sce],0.001,0.001,sigma_zphot_clustering[sce],0.001,0.001,3.0,1.2,3.8,2.0,16.0,5.0,0.8);


  double NORM;
  like.IA = 4; 
/* sample parameter values */
  // nuisance.oneplusz0_ia = 1.62;

/* here, do your time-consuming job */
  // init_cosmo_runmode("CLASS");
  init_cosmo_runmode("halofit");
//  init_bary("owls_AGN");
  cosmology.Omega_m   = 0.3;
  NORM    = 0.82355;
  // NORM  = 2.19e-09;
  cosmology.n_spec    = 0.97;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.048;
  cosmology.h0=0.69;
  double Omega_nuh2 = 0.0;//0.00083;//0.001743331232934258;
  cosmology.Omega_v=1.0-cosmology.Omega_m;
//  cosmology.theta_s  = -1.;
  double b1[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, b2[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},b_mag[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  b1[0] = 1.44;
  b1[1] = 1.70;
  b1[2] = 1.698;
  b1[3] = 1.997;
  b1[4] = 2.058;
//
//  b2[0] = 0.23; b2[1] = 0.23; b2[2] = 0.23; b2[3] = 0.5; b2[4] = .5;
  b_mag[0] = -0.102;
  b_mag[1] = -0.102;
  b_mag[2] = -0.102;
  b_mag[3] = 1.06;
  b_mag[4] = 1.06;

  double A_ia=0.5, eta_ia=0.;
  double p_ia[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  p_ia[0] = A_ia; p_ia[1] = eta_ia;

  double mean_m[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double sigma_m[10]={0.005,0.005,0.005,0.005,0.005,0.0,0.0,0.0,0.0,0.0};
  double bias_photoz_s[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  // double sigma_b_photoz_s[10]={0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08};
  double bias_photoz_l[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  // double sigma_b_photoz_l[10]={0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};

  // init_source_sample_mpp("./zdistris/nz_v0.16_smoothed.txt",4);
  // init_lens_sample_mpp("./zdistris/nz_y3_redmagic_v0.5.1_wide_gold_2.2.1_combined_hd3_hl2_z_samp.txt",5,b1,b2,0.0);

  init_source_sample_mpp("./zdistris/source_DESY6.nz",5);
  init_lens_sample_mpp("./zdistris/mcal_1101_lens.nz",5,b1,b2,0.1);


  // init_binning_mpp(20,2.5,250.);

  init_probes("3x2pt");
  
  // set_shear_priors_mpp(mean_m,sigma_m);
  //set_wlphotoz_priors_mpp(bias_photoz_s,sigma_b_photoz_s);
  //set_clphotoz_priors_mpp(bias_photoz_l,sigma_b_photoz_l);

  sprintf(like.MASK_FILE,"%s","none");
  printf("PATH TO MASK: %s\n",like.MASK_FILE);
  begin = clock();
  char datavfile[200];
  sprintf(datavfile,"datav/xi_Y6_3x2pt");
//  sprintf(datavfile,"datav/xi_Y3_baseline+b2_MICE+bary_owls_AGN");
  compute_data_vector(datavfile,cosmology.Omega_m,NORM ,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,Omega_nuh2,cosmology.h0,0.0,0.0,cosmology.theta_s,
    b1, b_mag,
    bias_photoz_s, //source photo-z bias
    bias_photoz_l, //lens photo-z bias
    mean_m, //shear calibration
    p_ia); // IA
  
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
   printf("time spent %le\n",time_spent);
   printf("like.IA = %d\n", like.IA);

}

void test_Cls_desy6_3src()
{
  clock_t begin, end;
  double time_spent;

  int i;

  init_binning_fourier(20, 30., 3000.);
  init_scalecuts(20., 3000.); // Rmin_bias = 20Mpc/h, lmax_shear=3000


  sprintf(survey.name,"%s","DESY6_3src");

  // init_bary(argv[2]);

  // init_priors(0.002,sigma_zphot_shear[sce],0.001,0.001,sigma_zphot_clustering[sce],0.001,0.001,3.0,1.2,3.8,2.0,16.0,5.0,0.8);


  double NORM;
  like.IA = 4; 
/* sample parameter values */
  // nuisance.oneplusz0_ia = 1.62;

/* here, do your time-consuming job */
  // init_cosmo_runmode("CLASS");
  init_cosmo_runmode("halofit");
//  init_bary("owls_AGN");
  cosmology.Omega_m   = 0.3;
  NORM    = 0.82355;
  // NORM  = 2.19e-09;
  cosmology.n_spec    = 0.97;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.048;
  cosmology.h0=0.69;
  double Omega_nuh2 = 0.0;//0.00083;//0.001743331232934258;
  cosmology.Omega_v=1.0-cosmology.Omega_m;
//  cosmology.theta_s  = -1.;
  double b1[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, b2[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},b_mag[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  b1[0] = 1.44;
  b1[1] = 1.70;
  b1[2] = 1.698;
  b1[3] = 1.997;
  b1[4] = 2.058;
//
//  b2[0] = 0.23; b2[1] = 0.23; b2[2] = 0.23; b2[3] = 0.5; b2[4] = .5;
  b_mag[0] = -0.102;
  b_mag[1] = -0.102;
  b_mag[2] = -0.102;
  b_mag[3] = 1.06;
  b_mag[4] = 1.06;

  double A_ia=0.5, eta_ia=0.;
  double p_ia[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  p_ia[0] = A_ia; p_ia[1] = eta_ia;

  double mean_m[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double sigma_m[10]={0.005,0.005,0.005,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double bias_photoz_s[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  // double sigma_b_photoz_s[10]={0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08};
  double bias_photoz_l[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  // double sigma_b_photoz_l[10]={0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};

  // init_source_sample_mpp("./zdistris/nz_v0.16_smoothed.txt",4);
  // init_lens_sample_mpp("./zdistris/nz_y3_redmagic_v0.5.1_wide_gold_2.2.1_combined_hd3_hl2_z_samp.txt",5,b1,b2,0.0);

  init_source_sample_mpp("./zdistris/source_DESY6_3bins.nz",3);
  init_lens_sample_mpp("./zdistris/mcal_1101_lens.nz",5,b1,b2,0.0);


  // init_binning_mpp(20,2.5,250.);

  init_probes("3x2pt");
  
  // set_shear_priors_mpp(mean_m,sigma_m);
  //set_wlphotoz_priors_mpp(bias_photoz_s,sigma_b_photoz_s);
  //set_clphotoz_priors_mpp(bias_photoz_l,sigma_b_photoz_l);

  sprintf(like.MASK_FILE,"%s","none");
  printf("PATH TO MASK: %s\n",like.MASK_FILE);
  begin = clock();
  char datavfile[200];
  sprintf(datavfile,"datav/xi_Y6_3x2pt_3src");
//  sprintf(datavfile,"datav/xi_Y3_baseline+b2_MICE+bary_owls_AGN");
  compute_data_vector(datavfile,cosmology.Omega_m,NORM ,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,Omega_nuh2,cosmology.h0,0.0,0.0,cosmology.theta_s,
    b1, b_mag,
    bias_photoz_s, //source photo-z bias
    bias_photoz_l, //lens photo-z bias
    mean_m, //shear calibration
    p_ia); // IA
  
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
   printf("time spent %le\n",time_spent);
   printf("like.IA = %d\n", like.IA);

}

void test_likelihood_desy3_planck()
{
  clock_t begin, end;
  double time_spent;

  int i;
  sprintf(like.DATA_FILE, "datav/xi_Y3_test");
  sprintf(like.COV_FILE, "covs/cov_desy3xplanck");

  init_binning_fourier(20, 30., 3000.);

  sprintf(survey.name,"%s","DESY3xplanck");

  // init_bary(argv[2]);

  // init_priors(0.002,sigma_zphot_shear[sce],0.001,0.001,sigma_zphot_clustering[sce],0.001,0.001,3.0,1.2,3.8,2.0,16.0,5.0,0.8);


  double NORM;
  like.IA = 4; 
/* sample parameter values */
  // nuisance.oneplusz0_ia = 1.62;

/* here, do your time-consuming job */
  // init_cosmo_runmode("CLASS");
  init_cosmo_runmode("halofit");
//  init_bary("owls_AGN");
  cosmology.Omega_m   = 0.3;
  NORM    = 0.82355;
  // NORM  = 2.19e-09;
  cosmology.n_spec    = 0.97;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.048;
  cosmology.h0=0.69;
  double Omega_nuh2 = 0.0;//0.00083;//0.001743331232934258;
  cosmology.Omega_v=1.0-cosmology.Omega_m;
//  cosmology.theta_s  = -1.;
  double b1[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, b2[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},b_mag[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  b1[0] = 1.7;
  b1[1] = 1.7;
  b1[2] = 1.7;
  b1[3] = 2.0;
  b1[4] = 2.0;
//
//  b2[0] = 0.23; b2[1] = 0.23; b2[2] = 0.23; b2[3] = 0.5; b2[4] = .5;
  b_mag[0] = -0.19375;
  b_mag[1] = -0.6285407;
  b_mag[2] = -0.69319886;
  b_mag[3] = 1.17735723;
  b_mag[4] = 1.87509758;

  double A_ia=0., eta_ia=0.;
  double p_ia[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  p_ia[0] = A_ia; p_ia[1] = eta_ia;

  double mean_m[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double sigma_m[10]={0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03};
  double bias_photoz_s[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  // double sigma_b_photoz_s[10]={0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08};
  double bias_photoz_l[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  // double sigma_b_photoz_l[10]={0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};

  // init_source_sample_mpp("./zdistris/nz_v0.16_smoothed.txt",4);
  // init_lens_sample_mpp("./zdistris/nz_y3_redmagic_v0.5.1_wide_gold_2.2.1_combined_hd3_hl2_z_samp.txt",5,b1,b2,0.0);

  init_source_sample_mpp("./zdistris/mcal_1101_source.nz",4);
  init_lens_sample_mpp("./zdistris/mcal_1101_lens.nz",5,b1,b2,0.0);


  // init_binning_mpp(20,2.5,250.);

  init_probes("6x2pt");
  init_cmb("planck");
  
  // set_shear_priors_mpp(mean_m,sigma_m);
  //set_wlphotoz_priors_mpp(bias_photoz_s,sigma_b_photoz_s);
  //set_clphotoz_priors_mpp(bias_photoz_l,sigma_b_photoz_l);

  sprintf(like.MASK_FILE,"%s","none");
  printf("PATH TO MASK: %s\n",like.MASK_FILE);

//  sprintf(datavfile,"datav/xi_Y3_baseline+b2_MICE+bary_owls_AGN");
  double loglike;
  loglike = log_multi_like(cosmology.Omega_m,NORM ,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,Omega_nuh2,cosmology.h0,0.0,0.0,cosmology.theta_s,
    b1, b_mag,
    bias_photoz_s, //source photo-z bias
    bias_photoz_l, //lens photo-z bias
    mean_m, //shear calibration
    p_ia); // IA

   printf("loglike = %lg\n", loglike);

}

void test_likelihood_desy6_planck()
{
  clock_t begin, end;
  double time_spent;

  int i;
  sprintf(like.DATA_FILE, "datav/xi_Y6_6x2pt");
  sprintf(like.COV_FILE, "covs/cov_desy6xplanck");

  init_binning_fourier(20, 30., 3000.);
  init_scalecuts(20., 3000.);
  sprintf(survey.name,"%s","DESY6xplanck");

  // init_bary(argv[2]);

  // init_priors(0.002,sigma_zphot_shear[sce],0.001,0.001,sigma_zphot_clustering[sce],0.001,0.001,3.0,1.2,3.8,2.0,16.0,5.0,0.8);


  double NORM;
  like.IA = 4; 
/* sample parameter values */
  // nuisance.oneplusz0_ia = 1.62;

/* here, do your time-consuming job */
  // init_cosmo_runmode("CLASS");
  init_cosmo_runmode("halofit");
//  init_bary("owls_AGN");
  cosmology.Omega_m   = 0.3;
  NORM    = 0.82355;
  // NORM  = 2.19e-09;
  cosmology.n_spec    = 0.97;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.048;
  cosmology.h0=0.69;
  double Omega_nuh2 = 0.0;//0.00083;//0.001743331232934258;
  cosmology.Omega_v=1.0-cosmology.Omega_m;
//  cosmology.theta_s  = -1.;
//   double b1[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.} ,b_mag[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double b1[10]={1.12284, 1.59082, 1.45414, 1.61371, 1.75653, 2, 2, 2, 2, 2};

double b2[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//   b1[0] = 1.44;
//   b1[1] = 1.70;
//   b1[2] = 1.698;
//   b1[3] = 1.997;
//   b1[4] = 2.058;
// //
// //  b2[0] = 0.23; b2[1] = 0.23; b2[2] = 0.23; b2[3] = 0.5; b2[4] = .5;
//   b_mag[0] = -0.102;
//   b_mag[1] = -0.102;
//   b_mag[2] = -0.102;
//   b_mag[3] = 1.06;
//   b_mag[4] = 1.06;

  double A_ia=0.0, eta_ia=0.;
  double p_ia[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  p_ia[0] = A_ia; p_ia[1] = eta_ia;

  // double mean_m[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double sigma_m[10]={0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03};
  // double bias_photoz_s[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  // double sigma_b_photoz_s[10]={0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08};
  // double bias_photoz_l[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  // double sigma_b_photoz_l[10]={0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};

  // init_source_sample_mpp("./zdistris/nz_v0.16_smoothed.txt",4);
  // init_lens_sample_mpp("./zdistris/nz_y3_redmagic_v0.5.1_wide_gold_2.2.1_combined_hd3_hl2_z_samp.txt",5,b1,b2,0.0);

  init_source_sample_mpp("./zdistris/source_DESY6.nz",5);
  init_lens_sample_mpp("./zdistris/mcal_1101_lens.nz",5,b1,b2,0.0);


  // init_binning_mpp(20,2.5,250.);

  init_probes("6x2pt");
  init_cmb("planck");
  
  // set_shear_priors_mpp(mean_m,sigma_m);
  //set_wlphotoz_priors_mpp(bias_photoz_s,sigma_b_photoz_s);
  //set_clphotoz_priors_mpp(bias_photoz_l,sigma_b_photoz_l);

  sprintf(like.MASK_FILE,"%s","none");
  printf("PATH TO MASK: %s\n",like.MASK_FILE);

//  sprintf(datavfile,"datav/xi_Y3_baseline+b2_MICE+bary_owls_AGN");
  double loglike;
  // loglike = log_multi_like(cosmology.Omega_m,NORM ,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,Omega_nuh2,cosmology.h0,0.0,0.0,cosmology.theta_s,
  //   b1, b_mag,
  //   bias_photoz_s, //source photo-z bias
  //   bias_photoz_l, //lens photo-z bias
  //   mean_m, //shear calibration
  //   p_ia); // IA
  double b_mag[10] ={0.0763057, 0.0921815, 0.0315041, 0.114855, -0.0681918, 0, 0, 0, 0, 0};

 double bias_photoz_l[10]={0.000348045, 0.00580718, -0.00200309, 0.00210533, -0.0011209,0.0,0.0,0.0,0.0,0.0};
 double bias_photoz_s[10]={-0.00117484, -0.0010133, 0.000757794, -0.000348866, -0.0015606,0.0,0.0,0.0,0.0,0.0};
  double mean_m[10]={-0.00170711, 0.00136659, 0.00172648, -0.00149607, -0.0015889,0.0,0.0,0.0,0.0,0.0};
  p_ia[0] = 0.0168433; p_ia[1]= 0.0804398;
  loglike = log_multi_like(0.310246,0.815847,0.963305,cosmology.w0,cosmology.wa,0.0491755,Omega_nuh2,0.679244,0.0,0.0,cosmology.theta_s,
    b1, b_mag,
    bias_photoz_s, //source photo-z bias
    bias_photoz_l, //lens photo-z bias
    mean_m, //shear calibration
    p_ia); // IA
   printf("loglike = %lg\n", loglike);

}

int main(void)
{
  // test_Cls_desy3();
  test_Cls_desy6();
  // test_Cls_desy6_3src();
  // test_likelihood_desy6_planck();
  return 0;
}






