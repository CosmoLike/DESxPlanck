#include "like_fourier_3x2pt.c"

void test_Cls_Roman_DC1(int argc, char** argv)
{
  like.shearcalib=1;
  
  /********* parameter settings start *********/
  // set_cosmological_parameters_to_(argv[1],1);
  // set_survey_parameters_to_(argv[2],1);


  // ell binning
  double ell_min = 30, ell_max = 4000;
  int Nell = 15;

  // CMB band-power binning
  int l_min = 2, l_max = 2500;
  int Nbp = 0;
  char binmat_with_corr_file[500] = "./cmblensrec/plancksmica/pp_agr2_CMBmarged/binning_matrix_with_correction_table.txt";
  char ckk_offset_file[500] = "./cmblensrec/plancksmica/pp_agr2_CMBmarged/Ckk_bandpower_offset.txt";

  // CMB setting
  // NOTE: the scale-cuts and FWHM of Planck beam size is hard-coded
  // lmin/lmax_kappacmb = 40/2999, FWHM = 7 arcmin
  char cmbName[50] = "planck";
  sprintf(cmb.pathLensRecNoise, "./cmblensrec/plancksmica/cmb_lmax3000.txt");
  
  // scale-cuts
  double Rmin_bias = 20.0;
  double lmax_shear = 4000.0;
  double ggl_cut = 1.0e-20;

  // galaxy sample
  int ntomo_source = 8, ntomo_lens = 8;
  char source_nz[500] = "./zdistris/roman_50_8_0.03.nz";
  char lens_nz[500] = "./zdistris/roman_50_8_0.03.nz";

  // misc
  char runmode[50] = "halofit";
  char probes[50] = "3x2pt";
  int IA_model = 4; // 4 = NLA, power-law redshift evolution
  
  // data vector, mask, and covariance matrix
  char cov_file[500] = "./covs/cov_Roman_Fourier_3x2pt";
  //char cov_file[500] = "./yaml/Cl_Roman_3x2pt.mask";
  // char data_file[500] = "./datav/Cl_Roman_3x2pt.realvector";
  char data_file[500] = "./yaml/Cl_Roman_3x2pt.mask";
  char mask_file[500] = "./yaml/Cl_Roman_3x2pt.mask";
  char test_model_file[500] = "./datav/Cl_Roman_3x2pt_test.modelvector";
  // char baryon_pca_file[500] = "./datav/cosmic_shear_10sim.pca";
  char baryon_pca_file[500] = "./yaml/Cl_Roman_3x2pt.mask";

  // cosmological parameters
  input_cosmo_params_y3 ic = {
    .omega_m = 0.3156,
    .sigma_8 = 0.831,
    .n_s = 0.9645,
    .w0 = -1.0,
    .wa = 0.0,
    .omega_b = 0.0492,
    .omega_nuh2 = 0.0,
    .h0 = 0.6727,
    .MGSigma = 0.0,
    .MGmu = 0.0,
    .theta_s = 0.0104854,
  };

  // nuisance parameters
  double b2[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    input_nuisance_params_y3 in = {
    .bias = {1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1, 0.0, 0.0},
    //.b_mag = {1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    .b_mag = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    .lens_z_bias = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    .source_z_bias = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    .shear_m = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    //.p_ia = {0.6, -1.5, 1.0, 0.6, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0},
    .p_ia = {0.6, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    //.p_ia = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    .Q1 = 0.0, .Q2 = 0.0, .Q3 = 0.0,
  };

  /********* parameter setting end *********/

  clock_t begin, end;
  double time_spent;
  int i;
  begin = clock();

  // Initialization
  init_cosmo_runmode(runmode);
  init_source_sample_mpp(source_nz, ntomo_source);
  init_lens_sample_mpp(lens_nz, ntomo_lens, in.bias, b2, ggl_cut);
  init_binning_fourier(Nell, ell_min, ell_max);
  init_binning_bandpower(Nbp, l_min, l_max);
  init_binning_real(0, 0.0, 0.0);
  init_scalecuts(Rmin_bias, lmax_shear);
  init_probes_fourier(probes);
  init_cmb(cmbName);
  init_data_fourier(cov_file, mask_file, data_file);
  //init_data_bandpower(cov_file, mask_file, data_file, binmat_with_corr_file, 
  //  ckk_offset_file, baryon_pca_file);
  init_IA_mpp(IA_model);
  sprintf(survey.name,"%s","Roman_DC1");

  // calculate and write model vector
  printf("test Roman DC1\n\n");
  write_datavector_wrapper(test_model_file, ic, in);
  printf("model vector written to %s\n", test_model_file);

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("time spent %le\n",time_spent);

  /*for(int i=0;i<ntomo_source;i++){
    for(int j=0;j<ntomo_lens;j++){
      if(test_zoverlap(i,j)==0){printf("[test_zoverlap] ggl pair %d-%d excluded!\n", i, j);}
      if(test_zoverlap_c(i,j)==0){printf("[test_zoverlap_c] ggl pair %d-%d excluded!\n", i, j);}
    }
  }*/
    
    // Print Pdelta(k,a) for sanity check
/*
    FILE *F;
    F=fopen("Pdelta_table_cosmolike.txt", "w");
    if(F==NULL){
      printf("ERORR: Can not open file Pdelta_table.txt\nAborting...\n");
      exit(1);
    }
    for (int aa=0; aa<19; aa++){
      for (int kk=0; kk<41; kk++){
        double a = 0.1 + aa*0.05;
        double k = pow(10.0, -4+0.15*kk);
        fprintf(F, "%f %f %le\n", k, a, Pdelta(k*cosmology.coverH0, a));
      }
    }
    fclose(F);

    double Pk = Pdelta(0.1, 0.99);
    printf("P(0.1,0.99)=%le\n", Pk);
*/
}


int main(int argc, char**argv)
{
  test_Cls_Roman_DC1(argc, argv);
  return 0;
}






