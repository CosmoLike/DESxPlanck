opt_home := -std=c17 -Wno-missing-braces -Wno-missing-field-initializers \
-I/usr/local/include -L/usr/local/lib -lgsl -lfftw3 -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -L../cosmolike_core/class -lclass
opt_puma :=    -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib \
 -lgsl -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass -lfftw3

cfastpt_dir := ../cosmolike_core/cfastpt/
cfastpt := $(cfastpt_dir)cfastpt.c $(cfastpt_dir)utils.c $(cfastpt_dir)utils_complex.c

cfftlog_dir := ../cosmolike_core/cfftlog/
cfftlog := $(cfftlog_dir)cfftlog.c $(cfftlog_dir)utils.c $(cfftlog_dir)utils_complex.c

COV_BIN_SIMPLE := ../cosmolike_core/theory/covariances_binned_simple.c

# puma_lib_mix:
#	gcc -shared -o like_mix_6x2pt.so -fPIC like_mix_6x2pt.c -DSAMPLING $(opt_puma) $(cfftlog) $(cfastpt)

# puma_cov_real:
#	gcc compute_covariances_real_6x2pt.c -o ./compute_covariances_real_6x2pt $(opt_puma)


home: 
	make home_lib
	make home_datavs
	make home_cov_real

home_datavs:
	gcc like_test_6x2pt_mix.c -o ./like_test_6x2pt_mix $(opt_home)
	#gcc like_test_3x2pt_mix.c -o ./like_test_3x2pt_mix $(opt_home)

home_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_home)

home_cov_real:
	gcc compute_covariances_real_6x2pt.c -o ./compute_covariances_real_6x2pt $(opt_home)

home_cov_fourier_binned:
	gcc compute_covariances_fourier_6x2pt_binned.c -o ./compute_covariances_fourier_6x2pt_binned $(opt_home)

home_cov_planck:
	gcc compute_covariances_fourier_planck.c -o ./compute_covariances_fourier_planck $(opt_home)


home_lib:
	gcc -shared -o like_mix_6x2pt.so -fPIC like_mix_6x2pt.c -DSAMPLING $(opt_home) $(cfftlog) $(cfastpt)

################################################
hpc_mix:
	make hpc_lib_mix
	make hpc_datavs_mix
	make hpc_cov_mix


hpc_datavs_mix:
	gcc like_test_6x2pt_mix.c -o ./like_test_6x2pt_mix $(opt_puma) $(cfftlog) $(cfastpt)

hpc_lib_mix:
	gcc -shared -o like_mix_6x2pt.so -fPIC like_mix_6x2pt.c -DSAMPLING $(opt_puma) $(cfftlog) $(cfastpt)

hpc_cov_mix:
	gcc compute_covariances_real_6x2pt.c  -o ./compute_covariances_real_6x2pt $(opt_puma) $(cfftlog) $(cfastpt)

hpc_mix_clean:
	rm like_mix_6x2pt.so ./compute_covariances_real_6x2pt
