double invcov_mask(int READ, int ci, int cj);
double data_read(int READ, int ci);
double bary_read(int READ, int PC, int cj);
void init_data_inv_bary(char *INV_FILE, char *DATA_FILE, char *BARY_FILE);
void init_priors(double M_Prior, double SigZ_source, double DeltaZ_source_Prior, double SigZ_source_Prior, double SigZ_lens, double DeltaZ_lens_Prior, double SigZ_lens_Prior, double A_IA_Prior, double betaIA_Prior, double etaIA_Prior, double etaZIA_Prior, double Q1_Prior, double Q2_Prior, double Q3_Prior);
void init_survey(char *surveyname, double nsource, double nlens, double area);
// void init_galaxies(char *SOURCE_ZFILE, char *LENS_ZFILE, char *lensphotoz, char *sourcephotoz, char *tomo_binning_source, char *tomo_binning_lens);
void init_cosmo_runmode(char *runmode);
void init_binning_fourier(int Ncl, double lmin, double lmax);
void init_scalecuts(double Rmin_bias, double lmax_shear);

void init_probes(char *probes);

void init_data_fourier(char *COV_FILE, char *MASK_FILE, char *DATA_FILE);

void init_sample_theta_s();

void init_source_sample_mpp(char *multihisto_file, int Ntomo);
void init_lens_sample_mpp(char *multihisto_file, int Ntomo, double *b1, double *b2, double ggl_cut);
void init_IA_mpp(int N);

void init_cmb(char * cmbName);
void set_cmb_planck();
void set_cmb_cmbs4();
void set_cmb_so_Y5();
void set_cmb_so_Y1();

int get_N_tomo_shear(void);
int get_N_tomo_clustering(void);
int get_N_ggl(void);
int get_N_ell(void);
int get_N_data_masked(void);

double mask(int ci);

void set_equal_tomo_bins(int Ntomo)
{
  int k,j;
  double frac, zi;
  tomo.shear_Nbin=Ntomo;
  tomo.shear_Npowerspectra=(int) (Ntomo*(Ntomo+1)/2);
  zdistr_histo_1(0.1, NULL);
  int zbins =2000;
  double da = (redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.shear_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+zdistr_histo_1(zi, NULL);
  }
  
  tomo.shear_zmin[0] = redshift.shear_zdistrpar_zmin;
  tomo.shear_zmax[tomo.shear_Nbin-1] = redshift.shear_zdistrpar_zmax;
  printf("\n");
  printf("Source Sample - Tomographic Bin limits:\n");
  for(k=0;k<tomo.shear_Nbin-1;k++){
    frac=(k+1.)/(1.*Ntomo)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.shear_zmax[k] = redshift.shear_zdistrpar_zmin+j*da;
    tomo.shear_zmin[k+1] = redshift.shear_zdistrpar_zmin+j*da;
    printf("min=%le max=%le\n",tomo.shear_zmin[k],tomo.shear_zmax[k]);
  }
  printf("min=%le max=%le\n",tomo.shear_zmin[tomo.shear_Nbin-1],tomo.shear_zmax[tomo.shear_Nbin-1]);
  printf("redshift.shear_zdistrpar_zmin=%le max=%le\n",redshift.shear_zdistrpar_zmin,redshift.shear_zdistrpar_zmax);
  free_double_vector(sum,0,zbins);
}

void init_source_sample_()
{
  int i;
  if(strcmp(survey.sourcephotoz,"none")==0) redshift.shear_photoz=0;
  if(strcmp(survey.sourcephotoz,"voigt")==0) redshift.shear_photoz=1;
  if(strcmp(survey.sourcephotoz,"voigt_out")==0) redshift.shear_photoz=2;
  if(strcmp(survey.sourcephotoz,"gaussian")==0) redshift.shear_photoz=3;
  if(strcmp(survey.sourcephotoz,"multihisto")==0) redshift.shear_photoz=4;
  //printf("Source Sample Redshift Errors set to %s: redshift.shear_photoz=%d\n",survey.sourcephotoz,redshift.shear_photoz);

  if (redshift.shear_photoz!=4)  set_equal_tomo_bins(tomo.shear_Nbin);
  for (i=0;i<tomo.shear_Nbin; i++)
  {
    //printf("zmean_source=%f\n",zmean_source(i));
    nuisance.bias_zphot_shear[i]=0.0;
  }
  tomo.shear_Npowerspectra = tomo.shear_Nbin*(tomo.shear_Nbin+1)/2;
}

void init_lens_sample_()
{
  int i,j,n;

  if(strcmp(survey.lensphotoz,"none")==0) redshift.clustering_photoz=0;
  if(strcmp(survey.lensphotoz,"voigt")==0) redshift.clustering_photoz=1;
  if(strcmp(survey.lensphotoz,"voigt_out")==0) redshift.clustering_photoz=2;
  if(strcmp(survey.lensphotoz,"gaussian")==0) redshift.clustering_photoz=3;
  if(strcmp(survey.lensphotoz,"multihisto")==0) redshift.clustering_photoz=4;
  //  printf("Lens Sample Redshift Errors set to %s: redshift.clustering_photoz=%d\n",survey.lensphotoz,redshift.clustering_photoz);
  tomo.clustering_Npowerspectra = tomo.clustering_Nbin;
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    // printf("zmean_lens=%f\n",zmean(i));
    nuisance.bias_zphot_clustering[i]=0.0;
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
    }
  }
  tomo.ggl_Npowerspectra = n;printf("ggl number: %d\n", n);

  //call test_kmax once to initialize look-up tables at reference cosmology
  test_kmax(1000.,1);
  // printf("end of lens sample init\n");
}

void init_sample_theta_s(){
  like.theta_s = 1;
}

int get_N_tomo_shear(void){
  return tomo.shear_Nbin;
}
int get_N_tomo_clustering(void){
  return tomo.clustering_Nbin;
}
int get_N_ggl(void){
  return tomo.ggl_Npowerspectra;
}
int get_N_ell(void){
  return like.Ncl;
}
int get_N_data_masked(void){
  int N = 0;
  for (int i = 0; i < like.Ndata; i++){ N += mask(i);}
  return N;
}

int count_rows(char* filename,const char delimiter){
  FILE *file = fopen (filename, "r" );
  char line [1000];
  if (file != NULL) {
    fgets(line,sizeof line,file);
    fclose(file);
  }
  else{
    printf("count_rows: file %s not found.\nEXIT\n",filename);
    exit(1);
  }
  int count = 1;
  char *p;

  p = line;
  while (*p != '\0')
  {
    if (*p == delimiter){
        while (*p == delimiter){p++;}
        count++;
    }
      p++;
    }
   return count;
}


// void init_data_inv_bary(char *INV_FILE, char *DATA_FILE, char *BARY_FILE)
// {
//   double init;
//   printf("\n");
//   printf("---------------------------------------\n");
//   printf("Initializing data vector and covariance\n");
//   printf("---------------------------------------\n");

//   sprintf(like.INV_FILE,"%s",INV_FILE);
//   printf("PATH TO INVCOV: %s\n",like.INV_FILE);
//   sprintf(like.DATA_FILE,"%s",DATA_FILE);
//   printf("PATH TO DATA: %s\n",like.DATA_FILE);
//   sprintf(like.BARY_FILE,"%s",BARY_FILE);
//   printf("PATH TO BARYONS: %s\n",like.BARY_FILE);
//   init=data_read(0,1);
//   init=bary_read(0,1,1);
//   init=invcov_read(0,1,1);
// }

void init_data_fourier(char *COV_FILE, char *MASK_FILE, char *DATA_FILE)
{
  double init;
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing data vector and covariance\n");
  printf("---------------------------------------\n");

  sprintf(like.MASK_FILE,"%s",MASK_FILE);
  printf("PATH TO MASK: %s\n",like.MASK_FILE);
  init=mask(1);

  sprintf(like.COV_FILE,"%s",COV_FILE);
  printf("PATH TO COV: %s\n",like.COV_FILE);
  init=invcov_mask(0,1,1);

  sprintf(like.DATA_FILE,"%s",DATA_FILE);
  printf("PATH TO DATA: %s\n",like.DATA_FILE);
  init=data_read(0,1);
}

double mask(int ci) // For fourier space
{
  int i,intspace;
  static double *mask =0;
  if(mask ==0){
    int N = 0;
    FILE *F;
    mask  = create_double_vector(0, like.Ndata-1); 
    double *maskc;
    maskc  = create_double_vector(0, like.Ndata-1); 
    if (strcmp(like.MASK_FILE,"none")==0){
      for (i=0;i<like.Ndata; i++){
        mask[i] = 1.0;
        maskc[i] =1.0;
      }
    }
    else{
      F=fopen(like.MASK_FILE,"r");
      if (!F){
        printf("init.c: invcov_mask: like.MASK_FILE = %s not found!\nEXIT!\n",like.MASK_FILE);
        exit(1);
      }
      for (i=0;i<like.Ndata; i++){
        fscanf(F,"%d %le\n",&intspace,&mask[i]);
        maskc[i] = mask[i];
        N += mask[i];
      }
     fclose(F);
     printf("%d bins within angular mask\n",N);
     printf("like.pos_pos = %d, like.shear_pos = %d,like.shear_shear = %d, like.ks = %d, like.gk = %d\n\n",like.pos_pos,like.shear_pos,like.shear_shear,like.ks,like.gk); 
    }
     int N3x2pt, N5x2pt, N6x2pt;
     N3x2pt = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);    
     N5x2pt = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra+tomo.shear_Nbin+tomo.clustering_Nbin);
     N6x2pt = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra+tomo.shear_Nbin+tomo.clustering_Nbin+1);
    //test whether Ndata assumes 3x2pt or 5x2pt format
    //if so, mask out probes excluded from the analysis
     if (N == N3x2pt || N== N5x2pt || N == N6x2pt){
      if(like.shear_shear==0){
        printf("masking out shear-shear bins\n");
       for (i = 0; i< like.Ncl*tomo.shear_Npowerspectra;i++){mask[i] = 0.;}
      }
      if(like.pos_pos==0){
        N = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra);
        printf("masking out clustering bins\n");
        for (i = N; i< N+like.Ncl*tomo.clustering_Npowerspectra;i++){mask[i] = 0.;}
      }
      if(like.shear_pos==0){
        N = like.Ncl*tomo.shear_Npowerspectra;
        printf("masking out ggl bins\n");
        for (i = N; i <N+like.Ncl*tomo.ggl_Npowerspectra; i++){mask[i] = 0.;}
      }
    }
    //test whether Ndata 5x2pt format
    //if so, mask out probes excluded from the analysis
    if (like.Ndata == N5x2pt){
      if(like.ks==0){
        printf("masking out shear x kappa bins\n");
        N = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
        for (i = N; i <N+like.Ncl*tomo.shear_Nbin; i++){mask[i] = 0.;}
      }
      if(like.gk==0){
        printf("masking out galaxies x kappa bins\n");
        N = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra+tomo.shear_Nbin);
        for (i = N; i < N+like.Ncl*tomo.clustering_Nbin; i++){mask[i] = 0.;}
      }
    }
    N = 0;
    for (i=0;i<like.Ndata; i++){
      //printf("mask(%d) = %.1f (was %.1f before probe cut)\n",i,mask[i],maskc[i]);
      N +=  mask[i];
    }
    printf("%d data points left after masking probes\n",N);
    if (N == 0){
      printf("init.c: mask: no data points left\nEXIT\n");
      exit(1);
    }
    printf("READ MASK FILE\n");
  }
  return mask[ci];
}

double invcov_mask(int READ, int ci, int cj)
{
  int i,j,intspace;
  static double **inv =0;

  if(READ==0 || inv == 0){
    inv   = create_double_matrix(0, like.Ndata-1, 0, like.Ndata-1);  

    gsl_matrix * cov   = gsl_matrix_calloc(like.Ndata, like.Ndata);      
    printf("Ndata: %d\n", like.Ndata);
    // gsl_matrix * inv_c   = gsl_matrix_calloc(like.Ndata, like.Ndata);
    gsl_matrix_set_zero(cov);
   
    double cov_G,cov_NG,doublespace,m;

    FILE *F;
    int n_rows =count_rows(like.COV_FILE,' ');
    F=fopen(like.COV_FILE,"r");

    if (!F){printf("init_LSSxCMB.c: invcov_mask: like.COV_FILE = %s not found!\nEXIT!\n",like.COV_FILE);exit(1);}
    switch (n_rows){
      case 3: while (fscanf(F,"%d %d %le\n", &i, &j, &cov_G) ==3) {
              m = 1.0;
              if (i < like.Ndata && j < like.Ndata){
              // apply mask to off-diagonal covariance elements
                if (i!=j){m = mask(i)*mask(j);}
              //  printf("%d %d (/%d) %e\n",i,j,like.Ndata,cov_G);
                gsl_matrix_set(cov,i,j,(cov_G)*m);
                gsl_matrix_set(cov,j,i,(cov_G)*m);
              }
            } break;
      case 4: while (fscanf(F,"%d %d %le %le\n", &i, &j, &cov_G, &cov_NG) ==4) {
              m = 1.0;
              if (i < like.Ndata && j < like.Ndata){
              // apply mask to off-diagonal covariance elements
                if (i!=j){m = mask(i)*mask(j);}
                //printf("%d %d (/%d) %e %e\n",i,j,like.Ndata,cov_G,cov_NG);
                gsl_matrix_set(cov,i,j,(cov_G+cov_NG)*m);
                gsl_matrix_set(cov,j,i,(cov_G+cov_NG)*m);
              }
            } break;
      case 10: while (fscanf(F,"%d %d %le %le %d %d %d %d %le %le\n", &i, &j, &doublespace, &doublespace,&intspace,&intspace,&intspace,&intspace,&cov_G,&cov_NG) ==10) {
              m = 1.0;
              if (i < like.Ndata && j < like.Ndata){
              // apply mask to off-diagonal covariance elements
                if (i!=j){m = mask(i)*mask(j);}
                //printf("%d %d (/%d) %e %e\n",i,j,like.Ndata,cov_G,cov_NG);
                gsl_matrix_set(cov,i,j,(cov_G+cov_NG)*m);
                gsl_matrix_set(cov,j,i,(cov_G+cov_NG)*m);
              }
            } break;
      default: printf("init_LSSxCMB.c:invcov_mask: covariance file %s has %d columns - unsupported format!\nEXIT\n",like.COV_FILE,n_rows);exit(1);
    }
    fclose(F);
    printf("READ COV_FILE\n");
    //SVD_inversion(cov,inv_c,like.Ndata);
    invert_matrix_colesky(cov);

    for (i=0;i<like.Ndata; i++){
      for (j=0;j<like.Ndata; j++){
        //apply mask again, to make sure numerical errors in matrix inversion don't cause problems...
        //also, set diagonal elements corresponding to datavector elements outside mask to zero, so that these elements don't contribute to chi2
        inv[i][j] =gsl_matrix_get(cov,i,j)*mask(i)*mask(j);
      }
    }
    gsl_matrix_free(cov);
    // gsl_matrix_free(inv_c);
    printf("FINISHED BUILDING INV COV\n");
  }
 return inv[ci][cj];
}

// double invcov_read(int READ, int ci, int cj)
// {
//   int i,j,intspace;
//   static double **inv =0;

//   if(READ==0 || inv == 0){
//     inv   = create_double_matrix(0, like.Ndata-1, 0, like.Ndata-1);      
//     FILE *F;
//     F=fopen(like.INV_FILE,"r");
//     for (i=0;i<like.Ndata; i++){
//       for (j=0;j<like.Ndata; j++){
//        fscanf(F,"%d %d %le\n",&intspace,&intspace,&inv[i][j]);  
//      }
//    }
//    fclose(F);
//    printf("FINISHED READING COVARIANCE\n");
//  }    
//  return inv[ci][cj];
// }


double data_read(int READ, int ci)
{
  int i,intspace;
  static double *data = 0;
  
  if(READ==0 || data ==0){
    data  = create_double_vector(0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.DATA_FILE,"r");
    for (i=0;i<like.Ndata; i++){  
      fscanf(F,"%d %le\n",&intspace,&data[i]);
    }
    fclose(F);
    printf("FINISHED READING DATA VECTOR\n");
  }    
  return data[ci];
}

double bary_read(int READ, int PC, int cj)
{
  int i,j,intspace, N_PC=6;
  static double **bary =0;

  if(READ==0 || bary == 0){
    bary   = create_double_matrix(0, N_PC-1, 0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.BARY_FILE,"r");
    for (i=0;i<like.Ndata; i++){
      fscanf(F,"%d %le %le %le %le %le %le\n",&intspace,&bary[0][i],&bary[1][i],&bary[2][i],&bary[3][i],&bary[4][i],&bary[5][i]);  
    }
    fclose(F);
    printf("FINISHED READING BARYON MATRIX\n");
  }    
  return bary[PC][cj];
}


void init_cosmo_runmode(char *runmode)
{
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Initializing Standard Runmode/Cosmology\n");
  printf("-------------------------------------------\n");
  set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP();
  sprintf(pdeltaparams.runmode,"%s",runmode);
  printf("pdeltaparams.runmode =%s\n",pdeltaparams.runmode);
}

void init_binning_fourier(int Ncl, double lmin, double lmax)
{
  printf("-------------------------------------------\n");
  printf("Initializing Binning\n");
  printf("-------------------------------------------\n");
  
  like.Ncl=Ncl;
  like.lmin= lmin; //std=20
  like.lmax= lmax; //15,000
  
  printf("number of ell bins Ncl: %d\n",like.Ncl);
  printf("minimum ell: %le\n",like.lmin);
  printf("maximum ell: %le\n",like.lmax);
  printf("init_binning_mpp complete\n");
}

void init_priors(double M_Prior, double SigZ_source, double DeltaZ_source_Prior, double SigZ_source_Prior, double SigZ_lens, double DeltaZ_lens_Prior, double SigZ_lens_Prior, double A_ia_Prior, double beta_ia_Prior, double eta_ia_Prior, double etaZ_ia_Prior, double Q1_Prior, double Q2_Prior, double Q3_Prior)
{
  int i;

  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing observational priors for marginalization\n");
  printf("---------------------------------------\n");
  
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = M_Prior;
  }
  like.shearcalib=1;
  
  for (i=0;i<tomo.shear_Nbin; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=SigZ_source; 
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i]; 
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = DeltaZ_source_Prior;
    prior.sigma_zphot_shear[i][1]= SigZ_source_Prior;
  }
  like.wlphotoz=1;

  for (i=0;i<tomo.clustering_Nbin; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=SigZ_lens; 
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    // rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = DeltaZ_lens_Prior;
    prior.sigma_zphot_clustering[i][1]= SigZ_lens_Prior;
  }
  like.clphotoz=1;

  prior.A_ia[0]=5.92; 
  prior.A_ia[1]=A_ia_Prior; 
  
  prior.beta_ia[0]=1.1; 
  prior.beta_ia[1]=beta_ia_Prior; 
  
  prior.eta_ia[0]=-0.47; 
  prior.eta_ia[1]=eta_ia_Prior; 
  
  prior.eta_ia_highz[0]=0.0; 
  prior.eta_ia_highz[1]=etaZ_ia_Prior; 
  like.IA=1;

  prior.bary_Q1[0]=0.0; 
  prior.bary_Q1[1]=Q1_Prior; 

  prior.bary_Q2[0]=0.0; 
  prior.bary_Q2[1]=Q2_Prior; 

  prior.bary_Q3[0]=0.0; 
  prior.bary_Q3[1]=Q3_Prior; 
  like.baryons=1;

  printf("\n");
  printf("---------------------------------------\n");
  printf("Shear Calibration Prior\n");
  printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[0][0],prior.shear_calibration_m[i][1]);

  printf("\n");
  printf("---------------------------------------\n");
  printf("Photo-z priors Weak Lensing\n");
  printf("Delta_z=%le, Sigma (Delta_z)=%le\n",prior.bias_zphot_shear[0][0],prior.bias_zphot_shear[0][1]);
  printf("Sigma_z=%le, Sigma (Sigma_z)=%le\n",prior.sigma_zphot_shear[0][0],prior.sigma_zphot_shear[0][1]); 

  printf("\n");
  printf("---------------------------------------\n");
  printf("Photo-z priors Clustering\n");
  printf("Delta_z=%le, Sigma (Delta_z)=%le\n",prior.bias_zphot_clustering[0][0],prior.bias_zphot_clustering[0][1]);
  printf("Sigma_z=%le, Sigma (Sigma_z)=%le\n",prior.sigma_zphot_clustering[0][0],prior.sigma_zphot_clustering[0][1]);

  printf("\n");
  printf("---------------------------------------\n");
  printf("IA Priors\n");
  printf("A_IA=%le, A_IA_Prior=%le\n",prior.A_ia[0],prior.A_ia[1]);
  printf("beta_ia=%le, betaIA_Prior=%le\n",prior.beta_ia[0],prior.beta_ia[1]);
  printf("eta_ia=%le, etaIA_Prior=%le\n",prior.eta_ia[0],prior.eta_ia[1]);
  printf("eta_ia_highz=%le, etaZIA_Prior=%le\n",prior.eta_ia_highz[0],prior.eta_ia_highz[1]);

  printf("\n");
  printf("---------------------------------------\n");
  printf("Baryon Priors\n");
  printf("Q1=%le, Sigma (Q1)=%le\n",prior.bary_Q1[0],prior.bary_Q1[1]);
  printf("Q2=%le, Sigma (Q2)=%le\n",prior.bary_Q2[0],prior.bary_Q2[1]);
  printf("Q3=%le, Sigma (Q3)=%le\n",prior.bary_Q3[0],prior.bary_Q3[1]);
}
 

void init_survey(char *surveyname, double nsource, double nlens, double area)
{
  printf("\n");
  printf("-------------------------------\n");
  printf("Initializing Survey Parameters\n");
  printf("-------------------------------\n");
  
  survey.area   = area;
  survey.n_gal   = nsource;
  survey.n_lens=nlens;
  survey.sigma_e   = 0.37;
  sprintf(survey.name,"%s",surveyname);
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=24.5;

  printf("Survey set to %s\n",survey.name);
  printf("Survey area: %le deg^2\n",survey.area);
  printf("Source Galaxy Density: %le galaxies/arcmin^2\n",survey.n_gal); 
}


void init_probes(char *probes)
{
  printf("\n");
  printf("------------------------------\n");
  printf("Initializing Probes\n");
  printf("------------------------------\n"); 
  printf("like.Ncl=%d\n",like.Ncl);
  printf("tomo.shear_Npowerspectra=%d\n",tomo.shear_Npowerspectra);
  printf("tomo.ggl_Npowerspectra=%d\n",tomo.ggl_Npowerspectra);
  printf("tomo.clustering_Npowerspectra=%d\n",tomo.clustering_Npowerspectra);

  sprintf(like.probes,"%s",probes);
  if(strcmp(probes,"shear_shear")==0){
    like.Ndata=like.Ncl*tomo.shear_Npowerspectra;
    like.shear_shear=1;
    printf("Shear-Shear computation initialized\n");
  }
  if(strcmp(probes,"pos_pos")==0){
    like.Ndata= like.Ncl*tomo.clustering_Npowerspectra;
    like.pos_pos=1;
    printf("Position-Position computation initialized\n");
  }
  if(strcmp(probes,"ggl_cl")==0){
    like.Ndata=like.Ncl*(tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  }
  if(strcmp(probes,"3x2pt")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  } 
  if(strcmp(probes,"6x2pt")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    like.kk = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("CMBkappa-Shear computation initialized\n");
    printf("CMBkappa-Position computation initialized\n");
    printf("CMBkappa-CMBkappa computation initialized\n");
  }
  if(strcmp(probes,"gg_gk_gs")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    printf("Position-Position computation initialized\n");
    printf("Position-Shear computation initialized\n");
    printf("Position-CMBkappa computation initialized\n");
  }
  if(strcmp(probes,"kk_ks_ss")==0) {
    like.Ndata = like.Ncl * (1+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.kk = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Shear-Shear computation initialized\n");
    printf("CMBkappa-Shear computation initialized\n");
    printf("CMBkappa-CMBkappa computation initialized\n");
  }
  printf("Total number of data points like.Ndata=%d\n",like.Ndata);
}



/************ CMB Settings ***********/
void init_cmb(char * cmbName) {
   printf("\n");
   printf("-----------------------------------\n");
   printf("Initializing CMB\n");
   printf("-----------------------------------\n");
   
   printf("CMB survey: %s\n", cmbName);
   if (strcmp(cmbName, "planck")==0)
      set_cmb_planck();
   if (strcmp(cmbName, "cmbs4")==0)
      set_cmb_cmbs4();
   if (strcmp(cmbName, "so_Y1")==0)
      set_cmb_so_Y1();
   if (strcmp(cmbName, "so_Y5")==0)
      set_cmb_so_Y5();
}

void set_cmb_planck() {
   sprintf(cmb.name, "planck");
   cmb.pathLensRecNoise = "./cmblensrec/plancksmica/cmb_lmax3000.txt";
   cmb.fsky = 0.673;
   cmb.fwhm = 7.0 * constants.arcmin;
   like.lmax_kappacmb = 2999.;
   printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
}

void set_cmb_cmbs4() {
   sprintf(cmb.name, "cmbs4");
   cmb.pathLensRecNoise = "./cmblensrec/cmbs4/cmb_lmax3000.txt";
   printf("ini_LSSxCMB.c: set_cmb_cmbs4(): need to specify cmb.fsky for cmbs4\n"); exit(1);
   like.lmax_kappacmb = 2999.;
   printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
}

void set_cmb_so_Y5() {
   sprintf(cmb.name, "so_Y5");
   // cmb.fwhm = 1.4 * (constants.pi/180.) / 60.;
   // cmb.sensitivity = 18.*(constants.pi/180.)/60.;
   printf("ini_LSSxCMB.c: set_cmb_so_Y5(): need to specify cmb.fsky for SO_Y5\n"); exit(1);
   cmb.pathLensRecNoise = "./cmblensrec/so/YEAR5_2colformat_nlkk_v3_1_0deproj0_SENS1_fsky0p4_it_lT30-3000_lP30-5000.dat";
   printf("ini_LSSxCMB.c: set_cmb_so_Y5(): need to specify cmb.fsky for SO_Y5\n"); exit(1);
   like.lmax_kappacmb = 2999.;
   printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
}
void set_cmb_so_Y1() {
   sprintf(cmb.name, "so_Y1");
   // cmb.fwhm = 1.4 * (constants.pi/180.) / 60.;
   // cmb.sensitivity = 18.*(constants.pi/180.)/60.;
   printf("ini_LSSxCMB.c: set_cmb_so_Y1(): need to specify cmb.fsky for SO_Y1\n"); exit(1);
   cmb.pathLensRecNoise = "./cmblensrec/so/YEAR1_nlkk_SOlike_y1_tt_SENS1_qe_fsky0p4_lT30-3000.dat";
   printf("ini_LSSxCMB.c: set_cmb_so_Y1(): need to specify cmb.fsky for SO_Y1\n"); exit(1);
   like.lmax_kappacmb = 2999.;
   printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
}

void init_ggl_tomo(){
  if (tomo.clustering_Nbin ==0){
    printf("WARNING! init_mpp.c: init_ggl_tomo called while tomo.clustering_Nbin =0\n");
  }
  if (tomo.shear_Nbin ==0){
    printf("WARNING! init_mpp.c: init_ggl_tomo called while tomo.shear_Nbin =0\n");
  }
  int n = 0;
  for (int i = 0; i < tomo.clustering_Nbin; i++){
    for(int j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      //printf("GGL combinations zl=%d zs=%d accept=%d; <z_l> = %.3f, <z_s> = %.3f\n",i,j,test_zoverlap(i,j), zmean(i),zmean_source(j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}

void init_lens_sample_mpp(char *multihisto_file, int Ntomo, double *b1, double *b2, double ggl_cut)
{
  sprintf(redshift.clustering_REDSHIFT_FILE,"%s",multihisto_file);
  redshift.clustering_photoz=4;
  tomo.clustering_Nbin = Ntomo;
  tomo.clustering_Npowerspectra = tomo.clustering_Nbin;
  survey.ggl_overlap_cut = ggl_cut;
  printf("Lens redshifts: multi-histo file %s, containing %d tomography bins\n",multihisto_file,tomo.clustering_Nbin);
  pf_photoz(0.1,0);
  for (int i=0;i<tomo.clustering_Nbin; i++)
  {
    gbias.b1_function = & b1_per_bin;
 //   tomo.n_lens[i]= n_lens[i];
    gbias.b[i] = b1[i];
    gbias.b2[i] = b2[i];
    nuisance.bias_zphot_clustering[i]=0.0;
 //   printf("bin %d: <z_l>=%.3f, b_1=%.3f, b_2=%.3f\n",i,zmean(i),gbias.b[i],gbias.b2[i]);
  }
  init_ggl_tomo();
  printf("init_lens_sample_mpp complete\n");
}


void init_source_sample_mpp(char *multihisto_file, int Ntomo)
{
  sprintf(redshift.shear_REDSHIFT_FILE,"%s",multihisto_file);
  redshift.shear_photoz=4;
  tomo.shear_Nbin = Ntomo;
  tomo.shear_Npowerspectra = tomo.shear_Nbin*(tomo.shear_Nbin+1)/2;
  printf("Source redshifts: multi-histo file %s, containing %d tomography bins\n",multihisto_file,tomo.shear_Nbin);
  for (int i=0;i<tomo.shear_Nbin; i++)
  {
    printf("bin %d: <z_s>=%f\n",i,zmean_source(i));
    //tomo.n_source[i]= n_source[i];
    nuisance.bias_zphot_shear[i]=0.0;
  }
  printf("init_source_sample_mpp complete\n");
}


void init_IA_mpp(int N)
{  
  if(N ==3){
    like.IA = N;
    printf("Set like.IA =3: NLA with per-z-bin amplitude\nSupply one nuisance.A_z[] parameter per source tomo bin\n");
  }
  else if(N ==4){
    like.IA = N;
    printf("Set like.IA =4; NLA with powerlaw z-evolution\nSupply nuisance.A_ia and nuisance.eta_ia\n");
  }
  else if(N ==5){
    like.IA = N;
    printf("Set like.IA =5; TATT with per-z-bin amplitude\nSupply nuisance.A_z[], nuisance.b_ta_z[], nuisance.A2_z[] parameters per source tomo bin\n");
  }
  else if(N ==6){
    like.IA = N;
    printf("Set like.IA =6; TATT with powerlaw z-evolution\nSupply nuisance.A_ia, nuisance.eta_ia, nuisance.b_ta_z[0], nuisance.A2_ia and nuisance.eta_ia_tt\n");
  }
  else{
    printf("like.IA = %d not supported in des_mpp\nEXIT\n", N);
    exit(1);
  }
}


void init_scalecuts(double Rmin_bias, double lmax_shear){
  like.Rmin_bias = Rmin_bias;
  like.lmax_shear = lmax_shear;
  printf("Initializing Rmin_bias = %lg, lmax_shear = %lg\n", Rmin_bias, lmax_shear);
}
