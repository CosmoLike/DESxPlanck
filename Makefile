opt_home := -std=c17 -Wno-missing-braces -Wno-missing-field-initializers \
-I/usr/local/include -L/usr/local/lib -lgsl -lfftw3 -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -L../cosmolike_core/class -lclass
opt_ocelote := -std=c17 -Wno-missing-braces -Wno-missing-field-initializers \
-I/cm/shared/uaapps/gsl/2.1/include -L/cm/shared/uaapps/gsl/2.1/lib \
-lfftw3 -lgsl -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
opt_puma := -std=c11 -Wno-missing-braces -Wno-missing-field-initializers \
-I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib \
-lfftw3 -lgsl -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
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
ocelote:
	make ocelote_lib
	make ocelote_datavs
	make ocelote_cov

ocelote_datavs:
	gcc like_fourier_6x2pt.c -o ./like_fourier_6x2pt $(opt_ocelote)

ocelote_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_ocelote)

ocelote_cov_planck:
	gcc compute_covariances_fourier_planck.c -o ./compute_covariances_fourier_planck $(opt_ocelote)

ocelote_lib:
	gcc -shared -o like_fourier_6x2pt.so -fPIC like_fourier_6x2pt.c -DSAMPLING $(opt_ocelote) $(cfftlog) $(cfastpt)


ocelote_cov_real:
	gcc compute_covariances_real_6x2pt.c -o ./compute_covariances_real_6x2pt $(opt_ocelote)

##########################################
puma_fourier:
	make puma_lib
	#make puma_datavs
	make puma_cov

puma_mix:
	make puma_lib_mix
	make puma_datavs
	make puma_cov_real


puma_datavs:
	gcc like_test_6x2pt_mix.c -o ./like_test_6x2pt_mix $(opt_puma)

puma_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_puma)

puma_cov_planck:
	gcc compute_covariances_fourier_planck.c -o ./compute_covariances_fourier_planck $(opt_puma)

puma_lib:
	gcc -shared -o like_fourier_6x2pt.so -fPIC like_fourier_6x2pt.c -DSAMPLING $(opt_puma) $(cfftlog) $(cfastpt)

puma_lib_mix:
	gcc -shared -o like_mix_6x2pt.so -fPIC like_mix_6x2pt.c -DSAMPLING $(opt_puma) $(cfftlog) $(cfastpt)

puma_cov_real:
	gcc compute_covariances_real_6x2pt.c  -o ./compute_covariances_real_6x2pt $(opt_puma)

puma_mix_clean:
	rm like_mix_6x2pt.so ./compute_covariances_real_6x2pt
