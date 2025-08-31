"""
Identifier:     csst_ifs_elfo/para.py
Name:           para.py
Description:    用拟合结果生成拟合初始猜测文件qsopar_i_j.fits
Author:         Hui Guo
Created:        2025-02-08
Modified-History:
    2025-02-08, created
"""

from astropy.io import fits
import numpy as np
from astropy.table import Table
import os


def create_qsopar(i: int, j: int, par_cube: np.ndarray):
    """
    Generate a fitting parameter FITS file for the given spectral coordinates.

    This function creates a FITS file containing the fitting parameters for the spectrum
    at the specified coordinates (i, j). The `par_cube` contains the fitting results
    from neighboring points, which are used to generate the parameter file.

    Parameters
    ----------
    i : int
        The row index of the spectral coordinates.
    j : int
        The column index of the spectral coordinates.
    par_cube : np.ndarray
        The fitting results from neighboring points, used to generate the fitting parameter file.

    Returns
    -------
    None
        Return None.

    Examples
    --------
    >>> create_qsopar(10, 15, par_cube_data)
    """

    path_ex = '.'  # os.path.join(pyqsofit.__path__[0], '..', 'example')

    # create a header
    hdr0 = fits.Header()
    hdr0['Author'] = 'Hengxiao Guo'
    primary_hdu = fits.PrimaryHDU(header=hdr0)
    """
    In this table, we par_cubeify the priors / initial conditions
    and boundaries for the line fitting parameters.

    """

    line_priors = np.rec.array([
        (np.e**par_cube['Ha_br1_1_centerwave'], 'Ha', 6400, 6800, 'Ha_br1', 1,
         par_cube['Ha_br1_1_scale'],
         0.0, 1e10, par_cube['Ha_br1_1_sigma'], 0.0017, 0.05, 0.015, 0, 0, 0,
         0.05, 1),
        (np.e**par_cube['Ha_br2_1_centerwave'], 'Ha', 6400, 6800, 'Ha_br2', 1,
         par_cube['Ha_br2_1_scale'],
         0.0, 1e10, par_cube['Ha_br2_1_sigma'], 0.0017, 0.05, 0.015, 0, 0, 0,
         0.05, 1),
        (np.e**par_cube['Ha_na_1_centerwave'], 'Ha', 6400, 6800, 'Ha_na', 1,
         par_cube['Ha_na_1_scale'],
         0.0, 1e10, par_cube['Ha_na_1_sigma'], 1.2e-4, 0.00169, 0.01, 1, 1, 0,
         0.002, 1),
        (np.e**par_cube['NII6549_1_centerwave'], 'Ha', 6400, 6800, 'NII6549',
         1, par_cube['NII6549_1_scale'],
         0.0, 1e10, par_cube['NII6549_1_sigma'], 1.2e-4, 0.00169, 5e-3, 1, 1,
         1, 0.001, 1),
        (np.e**par_cube['NII6585_1_centerwave'], 'Ha', 6400, 6800, 'NII6585',
         1, par_cube['NII6585_1_scale'],
         0.0, 1e10, par_cube['NII6585_1_sigma'], 1.2e-4, 0.00169, 5e-3, 1, 1,
         1, 0.003, 1),
        (np.e**par_cube['SII6718_1_centerwave'], 'Ha', 6400, 6800, 'SII6718',
         1, par_cube['SII6718_1_scale'],
         0.0, 1e10, par_cube['SII6718_1_sigma'], 1.2e-4, 0.00169, 5e-3, 1, 1,
         2, 0.001, 1),
        (np.e**par_cube['SII6732_1_centerwave'], 'Ha', 6400, 6800, 'SII6732',
         1, par_cube['SII6732_1_scale'],
         0.0, 1e10, par_cube['SII6732_1_sigma'], 1.2e-4, 0.00169, 5e-3, 1, 1,
         2, 0.001, 1),

        (np.e**par_cube['Hb_br1_1_centerwave'], 'Hb', 4640, 5100, 'Hb_br1', 1,
         par_cube['Hb_br1_1_scale'],
         0.0, 1e10, par_cube['Hb_br1_1_sigma'], 0.0017, 0.05, 0.01, 0, 0, 0,
         0.01, 1),
        (np.e**par_cube['Hb_br2_1_centerwave'], 'Hb', 4640, 5100, 'Hb_br2', 1,
         par_cube['Hb_br2_1_scale'],
         0.0, 1e10, par_cube['Hb_br2_1_sigma'], 0.0017, 0.05, 0.01, 0, 0, 0,
         0.01, 1),
        (np.e**par_cube['Hb_na_1_centerwave'], 'Hb', 4640, 5100, 'Hb_na', 1,
         par_cube['Hb_na_1_scale'],
         0.0, 1e10, par_cube['Hb_na_1_sigma'], 2.4e-4, 0.00169, 0.01, 1, 1, 0,
         0.002, 1),
        (np.e**par_cube['OIII4959c_1_centerwave'], 'Hb', 4640, 5100,
         'OIII4959c', 1, par_cube['OIII4959c_1_scale'],
         0.0, 1e10, par_cube['OIII4959c_1_sigma'], 2.4e-4, 0.00169, 0.01, 1,
         1, 0, 0.002, 1),
        (np.e**par_cube['OIII5007c_1_centerwave'], 'Hb', 4640, 5100,
         'OIII5007c', 1, par_cube['OIII5007c_1_scale'],
         0.0, 1e10, par_cube['OIII5007c_1_sigma'], 2.4e-4, 0.00169, 0.01, 1,
         1, 0, 0.004, 1),
        (np.e**par_cube['OIII4959w_1_centerwave'], 'Hb', 4640, 5100,
         'OIII4959w', 1, par_cube['OIII4959w_1_scale'],
         0.0, 1e10, par_cube['OIII4959w_1_sigma'], 2.4e-4, 0.004, 0.01, 2, 2,
         0, 0.001, 1),
        (np.e**par_cube['OIII5007w_1_centerwave'], 'Hb', 4640, 5100,
         'OIII5007w', 1, par_cube['OIII5007w_1_scale'],
         0.0, 1e10, par_cube['OIII5007w_1_sigma'], 2.4e-4, 0.004, 0.01, 2, 2,
         0, 0.002, 1),
    ],

        formats=('float32,    a20,  float32, float32,      a20,  int32,'
                 ' float32, float32, float32, float32, float32, float32,'
                 ' float32,   int32,  int32,  int32, float32, int32'),
        names=(' lambda, compname,   minwav,  maxwav, linename, ngauss,'
               '  inisca,  minsca,  maxsca,  inisig,  minsig,  maxsig,'
               '    voff,  vindex, windex, findex,  fvalue,  vary'))

    # Header
    hdr1 = fits.Header()
    hdr1['lambda'] = 'Vacuum Wavelength in Ang'
    hdr1['minwav'] = 'Lower complex fitting wavelength range'
    hdr1['maxwav'] = 'Upper complex fitting wavelength range'
    hdr1['ngauss'] = 'Number of Gaussians for the line'

    # Can be set to negative for absorption lines if you want
    hdr1['inisca'] = 'Initial guess of line scale [flux]'
    hdr1['minsca'] = 'Lower range of line scale [flux]'
    hdr1['maxsca'] = 'Upper range of line scale [flux]'

    hdr1['inisig'] = 'Initial guess of linesigma [lnlambda]'
    hdr1['minsig'] = 'Lower range of line sigma [lnlambda]'
    hdr1['maxsig'] = 'Upper range of line sigma [lnlambda]'

    hdr1['voff  '] = (
        'Limits on velocity offset from the central wavelength [lnlambda]')
    hdr1['vindex'] = (
        'Entries w/ same NONZERO vindex constrained to have same velocity')
    hdr1['windex'] = (
        'Entries w/ same NONZERO windex constrained to have same width')
    hdr1['findex'] = (
        'Entries w/ same NONZERO findex have constrained flux ratios')
    hdr1['fvalue'] = 'Relative scale factor for entries w/ same findex'

    hdr1['vary'] = (
        'Whether or not to vary the parameter (set to 0 to fix the line'
        ' parameter to initial values)')

    # Save line info
    hdu1 = fits.BinTableHDU(data=line_priors, header=hdr1, name='line_priors')
    """
    In this table, we par_cubeify the windows and priors / initial conditions
    and boundaries for the continuum fitting parameters.

    """

    conti_windows = np.rec.array([
        (1150., 1170.),
        (1275., 1290.),
        (1350., 1360.),
        (1445., 1465.),
        (1690., 1705.),
        (1770., 1810.),
        (1970., 2400.),
        (2480., 2675.),
        (2925., 3400.),
        (3775., 3832.),
        (4000., 4050.),
        (4200., 4230.),
        (4435., 4640.),
        (5100., 5535.),
        (6005., 6035.),
        (6110., 6250.),
        (6800., 7000.),
        (7160., 7180.),
        (7500., 7800.),
        # Continuum fitting windows (to avoid emission line, etc.)  [AA]
        (8050., 8150.),
    ],
        formats='float32,  float32',
        names='min,     max')

    hdu2 = fits.BinTableHDU(data=conti_windows, name='conti_windows')

    conti_priors = np.rec.array([
        # Normalization of the MgII Fe template [flux]
        ('Fe_uv_norm',  0.0,   0.0,   1e10,  1),
        # FWHM of the MgII Fe template [AA]
        ('Fe_uv_FWHM',  3000,  1200,  18000, 1),
        # Wavelength shift of the MgII Fe template [lnlambda]
        ('Fe_uv_shift', 0.0,   -0.01, 0.01,  1),
        # Normalization of the Hbeta/Halpha Fe template [flux]
        ('Fe_op_norm',  0.0,   0.0,   1e10,  1),
        # FWHM of the Hbeta/Halpha Fe template [AA]
        ('Fe_op_FWHM',  3000,  1200,  18000, 1),
        # Wavelength shift of the Hbeta/Halpha Fe template [lnlambda]
        ('Fe_op_shift', 0.0,   -0.01, 0.01,  1),
        # Normalization of the power-law (PL)
        # continuum f_lambda = (lambda/3000)^-alpha
        ('PL_norm',     1.0,   0.0,   1e10,  1),
        # Slope of the power-law (PL) continuum
        ('PL_slope',    -1.5,  -5.0,  3.0,   1),
        # Normalization of the Balmer continuum
        #  at < 3646 AA [flux] (Dietrich et al. 2002)
        ('Blamer_norm', 0.0,   0.0,   1e10,  1),
        # Te of the Balmer continuum at < 3646 AA [K?]
        ('Balmer_Te',   15000, 10000, 50000, 1),
        # Tau of the Balmer continuum at < 3646 AA
        ('Balmer_Tau',  0.5,   0.1,   2.0,   1),
        # 1st coefficient of the polynomial continuum
        ('conti_a_0',   0.0,   None,  None,  1),
        # 2nd coefficient of the polynomial continuum
        ('conti_a_1',   0.0,   None,  None,  1),
        # 3rd coefficient of the polynomial continuum
        ('conti_a_2',   0.0,   None,  None,  1),
        # Note: The min/max bounds on the conti_a_0
        #  coefficients are ignored by the code,
        # so they can be determined automatically for numerical stability.
    ],

        formats='a20,  float32, float32, float32, int32',
        names='parname, initial,   min,     max,     vary')

    hdr3 = fits.Header()
    hdr3['ini'] = 'Initial guess of line scale [flux]'
    hdr3['min'] = 'FWHM of the MgII Fe template'
    hdr3['max'] = 'Wavelength shift of the MgII Fe template'

    hdr3['vary'] = ('Whether or not to vary the parameter (set to 0 to fix'
                    ' the continuum parameter to initial values)')

    hdu3 = fits.BinTableHDU(
        data=conti_priors, header=hdr3, name='conti_priors')
    """
    In this table, we allow user to customized some key parameters in
      our result measurements.

    """

    measure_info = Table(
        [
            [[1350, 1450, 3000, 4200, 5100]],
            [[
                # [2240, 2650],
                [4435, 4685],
            ]]
        ],
        names=([
            'cont_loc',
            'Fe_flux_range'
        ]),
        dtype=([
            'float32',
            'float32'
        ])
    )
    hdr4 = fits.Header()
    hdr4['cont_loc'] = 'The wavelength of continuum luminosity in results'
    hdr4['Fe_flux_range'] = (
        'Fe emission wavelength range calculated in results')

    hdu4 = fits.BinTableHDU(
        data=measure_info, header=hdr4, name='measure_info')

    hdu_list = fits.HDUList([primary_hdu, hdu1, hdu2, hdu3, hdu4])
    hdu_list.writeto(os.path.join(
        path_ex, f'qsopar_{i}_{j}.fits'), overwrite=True)
