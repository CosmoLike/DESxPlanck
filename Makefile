opt_home := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/usr/local/include -L/usr/local/lib -lgsl -lfftw3 -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -L../cosmolike_core/class -lclass
opt_ocelote := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/cm/shared/uaapps/gsl/2.1/include -L/cm/shared/uaapps/gsl/2.1/lib \
-lfftw3 -lgsl -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
cfftlog_dir := ../cosmolike_core/cfftlog/
cfftlog := $(cfftlog_dir)cfftlog.c $(cfftlog_dir)utils.c $(cfftlog_dir)utils_complex.c


home: 
	make home_lib
	make home_datavs
	make home_cov

home_datavs:
	gcc like_test_6x2pt.c -o ./like_fourier_6x2pt $(opt_home)
	gcc like_test_3x2pt.c -o ./like_fourier_3x2pt $(opt_home)

home_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_home)

home_cov_real:
	gcc compute_covariances_real_6x2pt.c -o ./compute_covariances_real_6x2pt $(opt_home)

home_cov_fourier_binned:
	gcc compute_covariances_fourier_6x2pt_binned.c -o ./compute_covariances_fourier_6x2pt_binned $(opt_home)

home_cov_planck:
	gcc compute_covariances_fourier_planck.c -o ./compute_covariances_fourier_planck $(opt_home)


home_lib:
	gcc -shared -o like_fourier_6x2pt.so -fPIC like_fourier_6x2pt.c -DSAMPLING $(opt_home) $(cfftlog) $(cfastpt)


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

