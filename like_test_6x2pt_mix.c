#include "like_mix_6x2pt.c"

void test_mix_desy1_planck()
{
	like.shearcalib=1;
	
	/********* parameter settings start *********/

	// theta binning
	double theta_min = 2.5, theta_max = 250.0;
	int Ntheta = 20;

	// CMB band-power binning
	int l_min = 2, l_max = 2500;
	int Nbp = 14;
	char binmat_with_corr_file[500] = "./cmblensrec/plancksmica/pp_agr2_CMBmarged/binning_matrix_with_correction_table.txt";
	char ckk_offset_file[500] = "./cmblensrec/plancksmica/pp_agr2_CMBmarged/Ckk_bandpower_offset.txt";

	// CMB setting
	// NOTE: the scale-cuts and FWHM of Planck beam size is hard-coded
	// lmin/lmax_kappacmb = 40/2999, FWHM = 7 arcmin
	char cmbName[50] = "planck";
	char cmb_lens_noise_file[500] = "./cmblensrec/plancksmica/cmb_lmax3000.txt";

	// scale-cuts
	double Rmin_bias = 20.0;
	double lmax_shear = 3000.0;
	double ggl_cut = 0.0;

	// galaxy sample
	int ntomo_source = 4, ntomo_lens = 5;
	char source_nz[500] = "./zdistris/mcal_1101_source.nz";
	char lens_nz[500] = "./zdistris/mcal_1101_lens.nz";

	// misc
	char runmode[50] = "halofit";
	char probes[50] = "6x2pt";
	int IA_model = 4; // 4 = NLA, power-law redshift evolution
	
	// data vector, mask, and covariance matrix
	char cov_file[500] = "./covs/cov_y1xplanck_mix6x2pt_pp_p18cosmo_agr2_withAnnulus_kkkkSimDR3";
	char data_file[500] = "./datav/xi_desy1xplanck_6x2pt_realdata_20_wA_ref_pp_agr2";
	char mask_file[500] = "./yaml/xi_desy1xplanck_6x2pt_realdata_pp_agr2_CMBmarged.mask";
	char test_model_file[500] = "./datav/xi_desy1xplanck_6x2pt_test_fullsky";
	char baryon_pca_file[500] = "./datav/cosmic_shear_9sim.pca";

	// cosmological parameters
	input_cosmo_params_y3 ic = {
		.omega_m = 0.3,
		//.sigma_8 = 0.8191,
		.A_s = 2.10e-9,
		.n_s = 0.96605,
		.w0 = -1.0,
		.wa = 0.0,
		.omega_b = 0.04,
		.omega_nuh2 = 0.06/93.14,// neutrino mass [0, 0, 0.06] eV
		.h0 = 0.6732,
		.MGSigma = 0.0,
		.MGmu = 0.0,
		.theta_s = 0.0104854,
	};

	// nuisance parameters
	double b2[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	input_nuisance_params_y3 in = {
		.bias = {1.72716, 1.65168, 1.61423, 1.92886, 2.11633, 
				0.0, 0.0, 0.0, 0.0, 0.0},
		.b_mag = {-0.19375, -0.6285407, -0.69319886, 1.17735723, 1.87509758,
				0.0, 0.0, 0.0, 0.0, 0.0},
		.lens_z_bias = {0.008, -0.005, 0.006, 0.0, 0.0, 
			0.0, 0.0, 0.0, 0.0, 0.0},
		.source_z_bias = {-0.001, -0.019, 0.009, -0.018,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		.shear_m = {0.012, 0.012, 0.012, 0.012, 
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		.p_ia = {0.606102, -1.51541, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
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
	init_binning_fourier(0, 0.0, 0.0);
	init_binning_bandpower(Nbp, l_min, l_max);
	init_binning_real(Ntheta, theta_min, theta_max);
	init_scalecuts(Rmin_bias, lmax_shear);
	init_probes(probes);
	init_cmb(cmbName, cmb_lens_noise_file);
	//init_data_fourier(cov_file, mask_file, data_file);
	init_data_bandpower(cov_file, mask_file, data_file, binmat_with_corr_file, 
		ckk_offset_file, baryon_pca_file);
	init_IA_mpp(IA_model);
	sprintf(survey.name,"%s","DESY1xplanck");

	// calculate and write model vector
	printf("test3\n\n");
	write_datavector_wrapper(test_model_file, ic, in);
	printf("model vector written to %s\n", test_model_file);

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spent %le\n",time_spent);
}

int main(void)
{
	test_mix_desy1_planck();
	return 0;
}

