Data description

The data files in this folder are extracted from the public Planck data release

> smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_agr2_CMBmarged.dataset

some keywords:
- pp: minimum-variance (MV) estimated
- agr2: aggressive scale cuts, 8 <= L <= 2048
- CMBmarged: marginalized over primary CMB

Extracted using the python scripts provided by Karim Benabed

The original data products are

- _data vector_: `Ckk_bandpower_datavector.txt`, binned CLkk, i.e. the convergence power spectrum without the L prefactor.

- _data vector offset_: `Ckk_bandpower_offset.txt`, CLkk bias due to marginalization over primary CMB in the band-power likelihood. Should be added to the binned theoretical prediction.

- _band power covariance_: `Ckk_bandpower_covariance.txt`, numerical covariance matrix of the CLkk band-power, including the extra noise due to primary-CMB marginalization. **Note**: Hartlap factor should be applied to it when getting the precision matrix during likelihood evaluation.

- _binning definition_: `bin_definition.txt`, define the [L_min, L_max] of each band power bin.

- _binning matrix_: `binning_matrix_table.txt`, binning matrix without CMB-marginalization, used for covariance matrix calculation between CLkk and other configuration-space probes. The L-range of the matrix is [8, 2048]. An extended version of it is `binning_matrix_table_extend.txt`, where the matrix is padded to L in [2,2500] with zeros.

- _binning matrix with correction_: `binning_matrix_with_correction_table.txt`, binning matrix + corrections due to primary-CMB marginalization. This is used when calculating CLkk model vector.
