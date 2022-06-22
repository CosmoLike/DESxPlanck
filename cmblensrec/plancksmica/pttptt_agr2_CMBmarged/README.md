Data description

The data files in this folder are extracted from the public Planck data release

> smicadx12_Dec5_ftl_mv2_ndclpttptt_p_teb_agr2_CMBmarged.dataset

some keywords:
- pttptt: Temperature-only estimated
- agr2: aggressive scale cuts, 8 <= L <= 2048
- CMBmarged: marginalized over primary CMB

Extracted using the python scripts provided by Karim Benabed

The original data products are
- data vector: `Ckk_bandpower_datavector.txt`, binned as $\frac{(L(L+1))^2}{2\pi}\times C^{\phi\phi}_L$, but we multiply the data by $\pi/2$ to make it as $C^{\kappa\kappa}_L$
- data vector offset: `Ckk_bandpower_offset.txt`, same binning post-processing as the data vector. When get theoretical.
- band power covariance: `Ckk_bandpower_covariance.txt`, same binning post-processing as the data vector (except apply the half-pi factor twice).
- binning definition: `bin_definition.txt`, define the [L_min, L_max] of each band power bin
- binning matrix: `binning_matrix_table.txt`:
- binning matrix with correction: `binning_matrix_with_correction_table.txt`: binning matrix but marginalized over primary CMB. 
