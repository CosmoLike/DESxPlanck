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
	make home_shared
	make home_datav
	make home_cov

home_datav:
	gcc like_test_6x2pt.c -o ./like_fourier_6x2pt $(opt_home)

home_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_home)

home_shared:
	gcc -shared -o like_fourier.so -fPIC like_fourier.c $(opt_home)

ocelote:
	make ocelote_shared 
	make ocelote_datav
	make ocelote_cov

ocelote_datav:
	gcc like_fourier_6x2pt.c -o ./like_fourier_6x2pt $(opt_ocelote)

ocelote_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_ocelote)

ocelote_shared:
	gcc -shared -o like_fourier.so -fPIC like_fourier.c $(opt_ocelote)

