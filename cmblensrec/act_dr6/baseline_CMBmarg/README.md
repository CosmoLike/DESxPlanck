# Data description

The data files in this folder are extracted from the public [ACT DR6 CMB lensing likelihood dataset](https://lambda.gsfc.nasa.gov/product/act/actadv_dr6_lensing_lh_get.html)

- _CMBL reconstruction noise power spectrum_: `NLkk_ACT_DR6_lmax3000.txt`, noise property of the reconstructed CMBL convergence map (for cross-correlation usage). The NLkk in the original release is truncated at L=2100, we extend the public NLkk to L=3000 by measuring the CLkk bias between the reconstructed and the input fiducial CMBL convergence using the [400 Gaussian simulations](https://lambda.gsfc.nasa.gov/product/act/actadv_dr6_lensing_maps_info.html)

- _data vector_: `Ckk_bandpower_datavector.txt`, binned CLkk, i.e. the convergence power spectrum without the L prefactor.

- _data vector offset_: `Ckk_bandpower_offset.txt`, CLkk bias due to marginalization over primary CMB in the band-power likelihood. Should be added to the binned theoretical prediction. 

**Note:** unlike Planck, ACT does not has a bias in CLkk CMB-marginalization. We make up an empty file to follow the same convention

- _band power covariance_: `Ckk_bandpower_covariance.txt`, numerical covariance matrix of the CLkk band-power, including the extra noise due to primary-CMB marginalization. 

**Note**: Hartlap factor should be applied to it when getting the precision matrix during likelihood evaluation.

- _binning definition_: `bin_definition.txt`, define the [L_min, L_max] of each band power bin.

- _binning matrix_: `binning_matrix_table.txt`, binning matrix without CMB-marginalization, used for covariance matrix calculation between CLkk and other configuration-space probes. The L-range of the matrix is [8, 2048]. An extended version of it is `binning_matrix_table_extend.txt`, where the matrix is padded to L in [2,2500] with zeros.

- _binning matrix with correction_: `binning_matrix_with_correction_table.txt`, binning matrix + corrections due to primary-CMB marginalization. This is used when calculating CLkk model vector.

**Note**: ACT DR6 does differenciate between the original and CMB-marged binning matrix. They only include the effect in the inflated covariance matrix.
