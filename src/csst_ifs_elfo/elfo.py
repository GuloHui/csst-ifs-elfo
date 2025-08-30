"""
Identifier:     csst_ifs_elfo/elfo.py
Name:           elfo.py
Description:    "IFU emission line fitting optimization"
Author:         Hui Guo
Created:        2025-02-08
Modified-History:
    2025-02-08, created
"""
import time
import multiprocessing as mp
from .pyqsofit.PyQSOFit_changedwaverange import QSOFit
from astropy.io import fits
import numpy as np
import os.path
import os
from csst_ifs_elfo.para import create_qsopar
# import shutil  # 用于文件复制操作


def safe_open_fits(file_path):
    """
    安全地打开FITS文件。

    Args:
        file_path (str): FITS文件路径

    Returns:
        fits.HDU 或 None: 如果文件存在返回FITS数据,否则返回None
    """
    if os.path.isfile(file_path):
        return fits.open(file_path)[1].data
    else:
        return None


# def median(i, j, path_from):
#     """
#     计算3x3网格中心波长的中位数值。

#     Args:
#         i (int): 行索引
#         j (int): 列索引
#         path_from (str): 源文件路径

#     Returns:
#         tuple: (速度值, 选中的参数立方体, 坐标)
#     """
#     par_cube_list = []
#     for row in range(3):
#         for col in range(3):
#             i_coord = i + row - 1
#             j_coord = j + col - 1
#             par_cube = safe_open_fits(
#                 f'./{path_from}/{i_coord}_{j_coord}.fits')
#             if (
#                 par_cube is not None
#                 and '1_line_status' in par_cube.dtype.names
#                 and '2_line_status' in par_cube.dtype.names
#             ):
#                 vel = (par_cube['Ha_na_1_centerwave'] if 'Ha_na_1_centerwave'
#                        in par_cube.dtype.names else np.inf)
#                 par_cube_list.append((vel, par_cube, (i_coord, j_coord)))

#     # 确保至少有3个有效点
#     if len(par_cube_list) < 3:
#         return None, None, None

#     par_cube_list.sort(key=lambda x: x[0])
#     median_index = len(par_cube_list) // 2
#     best_specs = par_cube_list[median_index]
#     vel, par_cube_choosen, coords = best_specs

#     create_qsopar(i, j, par_cube_choosen)
#     return vel, par_cube_choosen, coords


# def flat(i, j, path_from):
#     """
#     计算3x3网格中心波长的中位数值。

#     Args:
#         i (int): 行索引
#         j (int): 列索引
#         path_from (str): 源文件路径

#     Returns:
#         tuple: (速度值, 选中的参数立方体, 坐标)
#     """
#     par_cube_list = []
#     for row in range(3):
#         for col in range(3):
#             i_coord = i + row - 1
#             j_coord = j + col - 1
#             par_cube = safe_open_fits(
#                 f'./{path_from}/{i_coord}_{j_coord}.fits')
#             if (
#                 par_cube is not None
#                 and '1_line_status' in par_cube.dtype.names
#                 and '2_line_status' in par_cube.dtype.names
#             ):
#                 chi2 = (par_cube['2_line_red_chi2'] if '2_line_red_chi2'
#                         in par_cube.dtype.names else np.inf)
#                 par_cube_list.append((chi2, par_cube, (i_coord, j_coord)))

#     par_cube_list.sort(key=lambda x: x[0])
#     best_specs = par_cube_list[0]
#     chi2, par_cube_choosen, coords = best_specs

#     create_qsopar(i, j, par_cube_choosen)
#     return chi2, par_cube_choosen, coords

def downward(i: int, j: int, jmin: int, jmax: int, path_out: str):
    """
    Generate the fitting parameter file for the specified coordinates.

    This function selects the best-fitting result with the smallest chi-square value
    from the spectral fitting of adjacent row (i-1) and columns (jmin to jmax)
    in the IFU data and generates the fitting parameter file for the corresponding (i, j) spatial pixel spectrum.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    j : int
        The column index for the current spectrum.
    jmin : int
        The minimum `j` coordinate for the range of adjacent columns to consider.
    jmax : int
        The maximum `j` coordinate for the range of adjacent columns to consider.
    path_out : str
        The path to the output directory where the fitting parameter files will be saved.

    Returns
    -------
    chi2 : float
        The chi-square value of the best-fitting spectrum.
    par_cube_choosen : `np.recarray`
        The parameter cube corresponding to the best-fitting spectrum.
    coords : tuple
        The coordinates `(i_coord, j_coord)` of the best-fitting spectrum.

    Examples
    --------
    >>> chi2, par_cube, coords = downward(10, 20, 5, 25, './output/')
    """
    par_cube_list = []
    for j_coord in range(jmin, jmax):
        i_coord = i - 1
        par_cube = safe_open_fits(
            f'./{path_out}/{i_coord}_{j_coord}.fits')
        if (
            par_cube is not None
            and '1_line_status' in par_cube.dtype.names
            and '2_line_status' in par_cube.dtype.names
        ):
            chi2 = (par_cube['2_line_red_chi2'] if '2_line_red_chi2'
                    in par_cube.dtype.names else np.inf)
            par_cube_list.append((chi2, par_cube, (i_coord, j_coord)))

    par_cube_list.sort(key=lambda x: x[0])
    best_specs = par_cube_list[0]
    chi2, par_cube_choosen, coords = best_specs

    create_qsopar(i, j, par_cube_choosen)
    return chi2, par_cube_choosen, coords


def upward(i: int, j: int, jmin: int, jmax: int, path_out: str):
    """
    Generate the fitting parameter file for the specified coordinates.

    This function selects the best-fitting result with the smallest chi-square value
    from the spectral fitting of adjacent row (i+1) and columns (jmin to jmax)
    in the IFU data and generates the fitting parameter file for the corresponding (i, j) spatial pixel spectrum.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    j : int
        The column index for the current spectrum.
    jmin : int
        The minimum `j` coordinate for the range of adjacent columns to consider.
    jmax : int
        The maximum `j` coordinate for the range of adjacent columns to consider.
    path_out : str
        The path to the output directory where the fitting parameter files will be saved.

    Returns
    -------
    chi2 : float
        The chi-square value of the best-fitting spectrum.
    par_cube_choosen : `np.recarray`
        The parameter cube corresponding to the best-fitting spectrum.
    coords : tuple
        The coordinates `(i_coord, j_coord)` of the best-fitting spectrum.

    Examples
    --------
    >>> chi2, par_cube, coords = upward(10, 20, 5, 25, './output/')
    """
    par_cube_list = []
    for j_coord in range(jmin, jmax):
        i_coord = i + 1
        par_cube = safe_open_fits(
            f'./{path_out}/{i_coord}_{j_coord}.fits')
        if (
            par_cube is not None
            and '1_line_status' in par_cube.dtype.names
            and '2_line_status' in par_cube.dtype.names
        ):
            chi2 = (par_cube['2_line_red_chi2'] if '2_line_red_chi2'
                    in par_cube.dtype.names else np.inf)
            par_cube_list.append((chi2, par_cube, (i_coord, j_coord)))

    par_cube_list.sort(key=lambda x: x[0])
    best_specs = par_cube_list[0]
    chi2, par_cube_choosen, coords = best_specs

    create_qsopar(i, j, par_cube_choosen)
    return chi2, par_cube_choosen, coords


def leftward(i: int, j: int, imin: int, imax: int, path_out: str):
    """
    Generate the fitting parameter file for the specified coordinates.

    This function selects the best-fitting result with the smallest chi-square value
    from the spectral fitting of adjacent column (j+1) and rows (imin to imax)
    in the IFU data and generates the fitting parameter file for the corresponding (i, j) spatial pixel spectrum.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    j : int
        The column index for the current spectrum.
    imin : int
        The minimum `i` coordinate for the range of adjacent columns to consider.
    imax : int
        The maximum `i` coordinate for the range of adjacent columns to consider.
    path_out : str
        The path to the output directory where the fitting parameter files will be saved.

    Returns
    -------
    chi2 : float
        The chi-square value of the best-fitting spectrum.
    par_cube_choosen : `np.recarray`
        The parameter cube corresponding to the best-fitting spectrum.
    coords : tuple
        The coordinates `(i_coord, j_coord)` of the best-fitting spectrum.

    Examples
    --------
    >>> chi2, par_cube, coords = leftward(10, 20, 5, 25, './output/')
    """
    par_cube_list = []
    for i_coord in range(imin, imax):
        j_coord = j + 1

        par_cube = safe_open_fits(
            f'./{path_out}/{i_coord}_{j_coord}.fits')
        if (
            par_cube is not None
            and '1_line_status' in par_cube.dtype.names
            and '2_line_status' in par_cube.dtype.names
        ):
            chi2 = (par_cube['2_line_red_chi2'] if '2_line_red_chi2'
                    in par_cube.dtype.names else np.inf)
            par_cube_list.append((chi2, par_cube, (i_coord, j_coord)))

    par_cube_list.sort(key=lambda x: x[0])
    best_specs = par_cube_list[0]
    chi2, par_cube_choosen, coords = best_specs

    create_qsopar(i, j, par_cube_choosen)
    return chi2, par_cube_choosen, coords


def rightward(i: int, j: int, imin: int, imax: int, path_out: str):
    """
    Generate the fitting parameter file for the specified coordinates.

    This function selects the best-fitting result with the smallest chi-square value
    from the spectral fitting of adjacent column (j-1) and rows (imin to imax)
    in the IFU data and generates the fitting parameter file for the corresponding (i, j) spatial pixel spectrum.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    j : int
        The column index for the current spectrum.
    imin : int
        The minimum `i` coordinate for the range of adjacent columns to consider.
    imax : int
        The maximum `i` coordinate for the range of adjacent columns to consider.
    path_out : str
        The path to the output directory where the fitting parameter files will be saved.

    Returns
    -------
    chi2 : float
        The chi-square value of the best-fitting spectrum.
    par_cube_choosen : `np.recarray`
        The parameter cube corresponding to the best-fitting spectrum.
    coords : tuple
        The coordinates `(i_coord, j_coord)` of the best-fitting spectrum.

    Examples
    --------
    >>> chi2, par_cube, coords = rightward(10, 20, 5, 25, './output/')
    """
    par_cube_list = []
    for i_coord in range(imin, imax):
        j_coord = j - 1
        par_cube = safe_open_fits(
            f'./{path_out}/{i_coord}_{j_coord}.fits')
        if (
            par_cube is not None
            and '1_line_status' in par_cube.dtype.names
            and '2_line_status' in par_cube.dtype.names
        ):
            chi2 = (par_cube['2_line_red_chi2'] if '2_line_red_chi2'
                    in par_cube.dtype.names else np.inf)
            par_cube_list.append((chi2, par_cube, (i_coord, j_coord)))

    par_cube_list.sort(key=lambda x: x[0])
    best_specs = par_cube_list[0]
    chi2, par_cube_choosen, coords = best_specs

    create_qsopar(i, j, par_cube_choosen)
    return chi2, par_cube_choosen, coords


def process_spectrum(i: int, j: int, path_out: str, wave: np.ndarray, flux: np.ndarray, err: np.ndarray, z: float, ra: float,
                     dec: float, path_ex: str):
    """
    Fit the spectrum for the given (i, j) coordinates.

    This function performs spectral fitting for the given row (i) and column (j)
    using the fitting parameters defined in a given parameter file.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    j : int
        The column index for the current spectrum.
    path_out : str
        The directory path where the fitting results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux : np.ndarray
        The flux array containing the spectral data.
    err : np.ndarray
        The error array containing the uncertainty of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        This function does not return any value but saves the fitting results to the specified output path.

    Examples
    --------
        >>> process_spectrum(5, 3, './output/', wave_array, flux_array, err_array, 0.064, 343.524658, -17.58185, './path_to_fits/')
    """
    try:

        lam = wave

        print(f"拟合光谱的坐标：i = {i}, j = {j}", end='\r')

        q_mle = QSOFit(lam, flux, err, z=z, path=path_ex, ra=ra, dec=dec)
        q_mle.Fit(name=str(i) + '_' + str(j)+'', nsmooth=1,
                  and_mask=False, or_mask=False, reject_badpix=False,
                  deredden=True,
                  wave_range=None, wave_mask=None, decompose_host=False,
                  host_line_mask=True,
                  Fe_uv_op=True, poly=True, BC=False, initial_guess=None,
                  rej_abs_conti=False,
                  linefit=True, rej_abs_line=False, MC=False, MCMC=False,
                  param_file_name='qsopar_muse.fits',
                  save_result=True, save_fits_path=path_out, plot_fig=True,
                  save_fig=True,
                  verbose=False,
                  kwargs_plot={'save_fig_path': path_out, 'broad_fwhm': 1200},)
        result = 1

    except Exception as e:
        result = 0
        # try:
        #     os.remove(f'qsopar_{i}_{j}.fits')
        # except:
        #     pass
        print(f"Error processing i={i}, j={j}: {e}")
    return result


# def process_spectrum_median(i, j, path_out, path_from, wave, flux,
#                             err, z, ra, dec, path_ex):
#     """
#     用周围3x3网格的中位数来作为光谱拟合的初始值


#     Args:
#         i (int): 行索引
#         j (int): 列索引
#         path_out (str): 输出路径
#         path_from (str): 源文件路径
#     """
#     try:

#         lam = wave

#         print(f"拟合光谱的坐标：i = {i}, j = {j}", end='\r')
#         median(i, j, path_from)

#         q_mle = QSOFit(lam, flux, err, z=z, path=path_ex, ra=ra, dec=dec)
#         q_mle.Fit(name=str(i) + '_' + str(j)+'', nsmooth=1,
#                   and_mask=False, or_mask=False, reject_badpix=False,
#                   deredden=True,
#                   wave_range=None, wave_mask=None, decompose_host=False,
#                   host_line_mask=True,
#                   Fe_uv_op=True, poly=True, BC=False, initial_guess=None,
#                   rej_abs_conti=False,
#                   linefit=True, rej_abs_line=False, MC=False, MCMC=False,
#                   param_file_name=f'qsopar_{i}_{j}.fits',
#                   save_result=True, save_fits_path=path_out, plot_fig=True,
#                   save_fig=True,
#                   verbose=False, kwargs_plot={'save_fig_path': path_out,
#                                               'broad_fwhm': 1200},)
#         os.remove(f'qsopar_{i}_{j}.fits')

#     except Exception as e:
#         try:
#             os.remove(f'qsopar_{i}_{j}.fits')
#         except OSError:
#             pass
#         print(f"Error processing i={i}, j={j}: {e}")


# def process_spectrum_flat(i, j, path_out, path_from, wave, flux,
#                           err, z, ra, dec, path_ex):
#     """
#     用周围3x3网格的中位数来作为光谱拟合的初始值


#     Args:
#         i (int): 行索引
#         j (int): 列索引
#         path_out (str): 输出路径
#         path_from (str): 源文件路径
#     """
#     try:

#         lam = wave

#         print(f"拟合光谱的坐标：i = {i}, j = {j}", end='\r')
#         flat(i, j, path_from)

#         q_mle = QSOFit(lam, flux, err, z=z, path=path_ex, ra=ra, dec=dec)
#         q_mle.Fit(name=str(i) + '_' + str(j)+'', nsmooth=1,
#                   and_mask=False, or_mask=False, reject_badpix=False,
#                   deredden=True,
#                   wave_range=None, wave_mask=None, decompose_host=False,
#                   host_line_mask=True,
#                   Fe_uv_op=True, poly=True, BC=False, initial_guess=None,
#                   rej_abs_conti=False,
#                   linefit=True, rej_abs_line=False, MC=False, MCMC=False,
#                   param_file_name=f'qsopar_{i}_{j}.fits',
#                   save_result=True, save_fits_path=path_out, plot_fig=True,
#                   save_fig=True,
#                   verbose=False,
#                   kwargs_plot={'save_fig_path': path_out, 'broad_fwhm': 1200},)
#         os.remove(f'qsopar_{i}_{j}.fits')

#     except Exception as e:
#         try:
#             os.remove(f'qsopar_{i}_{j}.fits')
#         except OSError:
#             pass
#         print(f"Error processing i={i}, j={j}: {e}")


def process_spectrum_down(i: int, j: int, path_out: str, wave: np.ndarray, flux: np.ndarray, err: np.ndarray, z: float, ra: float,
                          dec: float, path_ex: str):
    """
    Fit the spectrum for the given (i, j) coordinates.

    This function performs spectral fitting for the given row (i) and column (j)
    using the best fitting results of the adjacent row (i-1) in the range (jmin-jmax) as initial values.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    j : int
        The column index for the current spectrum.
    path_out : str
        The directory path where the fitting results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux : np.ndarray
        The flux array containing the spectral data.
    err : np.ndarray
        The error array containing the uncertainty of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        This function does not return any value but saves the fitting results to the specified output path.

    Examples
    --------
        >>> process_spectrum_down(5, 3, './output/', wave_array, flux_array, err_array, 0.064, 343.524658, -17.58185, './path_to_fits/')
    """
    try:

        lam = wave

        print(f"拟合光谱的坐标：i = {i}, j = {j}", end='\r')
        downward(i, j, j-4, j+5, path_out)

        q_mle = QSOFit(lam, flux, err, z=z, path=path_ex, ra=ra, dec=dec)
        q_mle.Fit(name=str(i) + '_' + str(j)+'', nsmooth=1,
                  and_mask=False, or_mask=False, reject_badpix=False,
                  deredden=True,
                  wave_range=None, wave_mask=None, decompose_host=False,
                  host_line_mask=True,
                  Fe_uv_op=True, poly=True, BC=False, initial_guess=None,
                  rej_abs_conti=False,
                  linefit=True, rej_abs_line=False, MC=False, MCMC=False,
                  param_file_name=f'qsopar_{i}_{j}.fits',
                  save_result=True, save_fits_path=path_out, plot_fig=True,
                  save_fig=True,
                  verbose=False,
                  kwargs_plot={'save_fig_path': path_out, 'broad_fwhm': 1200},)
        os.remove(f'qsopar_{i}_{j}.fits')
        result = 1

    except Exception as e:
        result = 0
        try:
            os.remove(f'qsopar_{i}_{j}.fits')
        except OSError:
            pass
        print(f"Error processing i={i}, j={j}: {e}")

    return result


def process_spectrum_up(i: int, j: int, path_out: str, wave: np.ndarray, flux: np.ndarray, err: np.ndarray, z: float, ra: float,
                        dec: float, path_ex: str):
    """
    Fit the spectrum for the given (i, j) coordinates.

    This function performs spectral fitting for the given row (i) and column (j)
    using the best fitting results of the adjacent row (i+1) in the range (jmin-jmax) as initial values.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    j : int
        The column index for the current spectrum.
    path_out : str
        The directory path where the fitting results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux : np.ndarray
        The flux array containing the spectral data.
    err : np.ndarray
        The error array containing the uncertainty of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        This function does not return any value but saves the fitting results to the specified output path.

    Examples
    --------
        >>> process_spectrum_up(5, 3, './output/', wave_array, flux_array, err_array, 0.064, 343.524658, -17.58185, './path_to_fits/')
    """
    try:

        lam = wave

        print(f"拟合光谱的坐标：i = {i}, j = {j}", end='\r')
        upward(i, j, j-4, j+5, path_out)

        q_mle = QSOFit(lam, flux, err, z=z, path=path_ex, ra=ra, dec=dec)
        q_mle.Fit(name=str(i) + '_' + str(j)+'', nsmooth=1,
                  and_mask=False, or_mask=False, reject_badpix=False,
                  deredden=True,
                  wave_range=None, wave_mask=None, decompose_host=False,
                  host_line_mask=True,
                  Fe_uv_op=True, poly=True, BC=False, initial_guess=None,
                  rej_abs_conti=False,
                  linefit=True, rej_abs_line=False, MC=False, MCMC=False,
                  param_file_name=f'qsopar_{i}_{j}.fits',
                  save_result=True, save_fits_path=path_out, plot_fig=True,
                  save_fig=True,
                  verbose=False,
                  kwargs_plot={'save_fig_path': path_out, 'broad_fwhm': 1200},)
        os.remove(f'qsopar_{i}_{j}.fits')
        result = 1

    except Exception as e:
        result = 0
        try:
            os.remove(f'qsopar_{i}_{j}.fits')
        except OSError:
            pass
        print(f"Error processing i={i}, j={j}: {e}")
    return result


def process_spectrum_left(i: int, j: int, path_out: str, wave: np.ndarray, flux: np.ndarray, err: np.ndarray, z: float, ra: float,
                          dec: float, path_ex: str):
    """
    Fit the spectrum for the given (i, j) coordinates.

    This function performs spectral fitting for the given row (i) and column (j)
    using the best fitting results of the adjacent column (j+1) in the range (imin-imax) as initial values.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    j : int
        The column index for the current spectrum.
    path_out : str
        The directory path where the fitting results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux : np.ndarray
        The flux array containing the spectral data.
    err : np.ndarray
        The error array containing the uncertainty of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        This function does not return any value but saves the fitting results to the specified output path.

    Examples
    --------
        >>> process_spectrum_left(5, 3, './output/', wave_array, flux_array, err_array, 0.064, 343.524658, -17.58185, './path_to_fits/')
    """
    try:

        lam = wave
        print(f"拟合光谱的坐标：i = {i}, j = {j}", end='\r')
        leftward(i, j, i-4, i+5, path_out)

        q_mle = QSOFit(lam, flux, err, z=z, path=path_ex, ra=ra, dec=dec)
        q_mle.Fit(name=str(i) + '_' + str(j)+'', nsmooth=1,
                  and_mask=False, or_mask=False, reject_badpix=False,
                  deredden=True,
                  wave_range=None, wave_mask=None, decompose_host=False,
                  host_line_mask=True,
                  Fe_uv_op=True, poly=True, BC=False, initial_guess=None,
                  rej_abs_conti=False,
                  linefit=True, rej_abs_line=False, MC=False, MCMC=False,
                  param_file_name=f'qsopar_{i}_{j}.fits',
                  save_result=True, save_fits_path=path_out, plot_fig=True,
                  save_fig=True,
                  verbose=False,
                  kwargs_plot={'save_fig_path': path_out, 'broad_fwhm': 1200},)
        os.remove(f'qsopar_{i}_{j}.fits')
        result = 1

    except Exception as e:
        result = 0
        try:
            os.remove(f'qsopar_{i}_{j}.fits')
        except OSError:
            pass
        print(f"Error processing i={i}, j={j}: {e}")
    return result


def process_spectrum_right(i: int, j: int, path_out: str, wave: np.ndarray, flux: np.ndarray, err: np.ndarray, z: float, ra: float,
                           dec: float, path_ex: str):
    """
    Fit the spectrum for the given (i, j) coordinates.

    This function performs spectral fitting for the given row (i) and column (j)
    using the best fitting results of the adjacent column (j-1) in the range (imin-imax) as initial values.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    j : int
        The column index for the current spectrum.
    path_out : str
        The directory path where the fitting results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux : np.ndarray
        The flux array containing the spectral data.
    err : np.ndarray
        The error array containing the uncertainty of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        This function does not return any value but saves the fitting results to the specified output path.

    Examples
    --------
        >>> process_spectrum_right(5, 3, './output/', wave_array, flux_array, err_array, 0.064, 343.524658, -17.58185, './path_to_fits/')
    """
    try:

        lam = wave

        print(f"拟合光谱的坐标：i = {i}, j = {j}", end='\r')
        rightward(i, j, i-4, i+5, path_out)

        q_mle = QSOFit(lam, flux, err, z=z, path=path_ex, ra=ra, dec=dec)
        q_mle.Fit(name=str(i) + '_' + str(j)+'', nsmooth=1,
                  and_mask=False, or_mask=False, reject_badpix=False,
                  deredden=True,
                  wave_range=None, wave_mask=None, decompose_host=False,
                  host_line_mask=True,
                  Fe_uv_op=True, poly=True, BC=False, initial_guess=None,
                  rej_abs_conti=False,
                  linefit=True, rej_abs_line=False, MC=False, MCMC=False,
                  param_file_name=f'qsopar_{i}_{j}.fits',
                  save_result=True, save_fits_path=path_out, plot_fig=True,
                  save_fig=True,
                  verbose=False,
                  kwargs_plot={'save_fig_path': path_out, 'broad_fwhm': 1200},)
        os.remove(f'qsopar_{i}_{j}.fits')
        result = 1

    except Exception as e:
        result = 0
        try:
            os.remove(f'qsopar_{i}_{j}.fits')
        except OSError:
            pass
        print(f"Error processing i={i}, j={j}: {e}")
    return result


def process_i(i: int, path_out: str, wave: np.ndarray, flux_cube: np.ndarray, var_cube: np.ndarray, z: float, ra: float, dec: float, path_ex: str):
    """
    Fit the spectrum for a given row and generate the fitting results.

    This function fits the spectrum of row `i` by using the fitting parameters
    defined in a given parameter file.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    path_out : str
        The directory path where the output results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux_cube : np.ndarray
        The flux cube containing the spectral data.
    var_cube : np.ndarray
        The variance cube of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        The function does not return any value but saves the results to the specified output path.

    Examples
    --------
    >>> process_i(5, './output/', wave_array, flux_cube, var_cube, 0.5, 150.7, -1.5, './path_to_fits/')
    """
    start_time_i = time.time()  # Record the start time of the outer loop
    print(f"开始i = {i}行 ")
    pool = mp.Pool(processes=5)  # Create a pool of processes

    # Use starmap to pass multiple arguments to the function
    pool.starmap(process_spectrum, [
        (i, j, path_out, wave, flux_cube[:, i, j], np.sqrt(var_cube[:, i, j]),
         z, ra, dec, path_ex)
        for j in range(0, flux_cube.shape[2])
    ])

    pool.close()

    pool.join()  # Wait for the worker processes to exit

    end_time_i = time.time()  # Record the end time of the outer loop
    # Calculate the elapsed time for the outer loop
    elapsed_time_i = end_time_i - start_time_i
    # Convert elapsed time to minutes and seconds
    minutes, seconds = divmod(elapsed_time_i, 60)
    print(f"i = {i} 的拟合时间 = {int(minutes)}m{seconds:.2f}s")
    return elapsed_time_i


def process_j(j: int, path_out: str, wave: np.ndarray, flux_cube: np.ndarray, var_cube: np.ndarray, z: float, ra: float, dec: float, path_ex: str):
    """
    Fit the spectrum for a given column and generate the fitting results.

    This function fits the spectrum of column `j` by using the fitting parameters
    defined in a given parameter file.

    Parameters
    ----------
    j : int
        The column index for the current spectrum.
    path_out : str
        The directory path where the output results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux_cube : np.ndarray
        The flux cube containing the spectral data.
    var_cube : np.ndarray
        The variance cube of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        The function does not return any value but saves the results to the specified output path.

    Examples
    --------
    >>> process_j(5, './output/', wave_array, flux_cube, var_cube, 0.5, 150.7, -1.5, './path_to_fits/')
    """
    start_time_j = time.time()  # Record the start time of the outer loop
    print(f"开始 j = {j} 列")
    pool = mp.Pool(processes=5)  # Create a pool of processes

    # Use starmap to pass multiple arguments to the function
    pool.starmap(process_spectrum, [
        (i, j, path_out, wave, flux_cube[:, i, j], np.sqrt(var_cube[:, i, j]),
         z, ra, dec, path_ex)
        for i in range(0, flux_cube.shape[1])
    ])
    pool.close()
    pool.join()  # Wait for the worker processes to exit

    end_time_j = time.time()  # Record the end time of the outer loop
    # Calculate the elapsed time for the outer loop
    elapsed_time_j = end_time_j - start_time_j
    # Convert elapsed time to minutes and seconds
    minutes, seconds = divmod(elapsed_time_j, 60)
    print(f"j = {j} 的拟合时间 = {int(minutes)}m{seconds:.2f}s")
    return elapsed_time_j


# def process_median(i, path_out, path_from, wave, flux_cube, var_cube, z, ra,
#                    dec, path_ex):
#     """
#     用周围3x3网格的中位数来作为光谱拟合的初始值
#     并行处理一行数据

#     Args:
#         i (int): 行索引
#         path_out (str): 输出路径
#         path_from (str): 源文件路径
#     """
#     start_time_i = time.time()  # Record the start time of the outer loop
#     print(f"开始i = {i}行 ")
#     pool = mp.Pool(processes=6)  # Create a pool of processes

#     pool.starmap(process_spectrum_median, [
#             (i, j, path_out, path_from, wave,
#              flux_cube[:, i, j], np.sqrt(var_cube[:, i, j]), z, ra, dec,
#              path_ex)
#             for j in range(0, flux_cube.shape[2])
#         ])

#     pool.close()
#     pool.join()  # Wait for the worker processes to exit

#     end_time_i = time.time()  # Record the end time of the outer loop
#     # Calculate the elapsed time for the outer loop
#     elapsed_time_i = end_time_i - start_time_i
#     # Convert elapsed time to minutes and seconds
#     minutes, seconds = divmod(elapsed_time_i, 60)
#     print(f"i = {i} 的拟合时间 = {int(minutes)}m{seconds:.2f}s")


# def process_flat(i, path_out, path_from, wave, flux_cube, var_cube, z, ra,
#                  dec, path_ex):
#     """
#     用周围3x3网格的中位数来作为光谱拟合的初始值
#     并行处理一行数据

#     Args:
#         i (int): 行索引
#         path_out (str): 输出路径
#         path_from (str): 源文件路径
#     """
#     start_time_i = time.time()  # Record the start time of the outer loop
#     print(f"开始i = {i}行 ")
#     pool = mp.Pool(processes=6)  # Create a pool of processes

#     pool.starmap(process_spectrum_flat, [
#             (i, j, path_out, path_from, wave,
#              flux_cube[:, i, j], np.sqrt(var_cube[:, i, j]), z, ra, dec,
#              path_ex)
#             for j in range(0, flux_cube.shape[2])
#         ])
#     pool.close()
#     pool.join()  # Wait for the worker processes to exit

#     end_time_i = time.time()  # Record the end time of the outer loop
#     # Calculate the elapsed time for the outer loop
#     elapsed_time_i = end_time_i - start_time_i
#     # Convert elapsed time to minutes and seconds
#     minutes, seconds = divmod(elapsed_time_i, 60)
#     print(f"i = {i} 的拟合时间 = {int(minutes)}m{seconds:.2f}s")


def process_right(j: int, path_out: str, wave: np.ndarray, flux_cube: np.ndarray, var_cube: np.ndarray, z: float, ra: float, dec: float, path_ex: str):
    """
    Fit the spectrum for a given column and generate the fitting results.

    This function fits the spectrum of column `j` using the fitting results
    of the neighboring columns in the specified direction.

    Parameters
    ----------
    j : int
        The column index for the current spectrum.
    path_out : str
        The directory path where the output results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux_cube : np.ndarray
        The flux cube containing the spectral data.
    var_cube : np.ndarray
        The variance cube of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        The function does not return any value but saves the results to the specified output path.

    Examples
    --------
    >>> process_right(5, './output/', wave_array, flux_cube, var_cube, 0.5, 150.7, -1.5, './path_to_fits/')
    """

    start_time_j = time.time()  # Record the start time of the outer loop

    print(f"开始j = {j}行 ")
    pool = mp.Pool(processes=6)  # Create a pool of processes

    pool.starmap(process_spectrum_right, [
            (i, j, path_out, wave, flux_cube[:, i, j],
             np.sqrt(var_cube[:, i, j]), z, ra, dec, path_ex)
            for i in range(0, flux_cube.shape[1])
        ])
    pool.close()
    pool.join()  # Wait for the worker processes to exit

    end_time_j = time.time()  # Record the end time of the outer loop
    # Calculate the elapsed time for the outer loop
    elapsed_time_j = end_time_j - start_time_j
    # Convert elapsed time to minutes and seconds
    minutes, seconds = divmod(elapsed_time_j, 60)
    print(f"j = {j} 的拟合时间 = {int(minutes)}m{seconds:.2f}s")
    return elapsed_time_j


def process_left(j: int, path_out: str, wave: np.ndarray, flux_cube: np.ndarray, var_cube: np.ndarray, z: float, ra: float, dec: float, path_ex: str):
    """
    Fit the spectrum for a given column and generate the fitting results.

    This function fits the spectrum of column `j` using the fitting
    results of the neighboring columns in the specified direction.

    Parameters
    ----------
    j : int
        The column index for the current spectrum.
    path_out : str
        The directory path where the output results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux_cube : np.ndarray
        The flux cube containing the spectral data.
    var_cube : np.ndarray
        The variance cube of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        The function does not return any value but saves the results to the specified output path.

    Examples
    --------
    >>> process_left(5, './output/', wave_array, flux_cube, var_cube, 0.5, 150.7, -1.5, './path_to_fits/')
    """

    start_time_j = time.time()  # Record the start time of the outer loop

    print(f"开始j = {j}行 ")
    pool = mp.Pool(processes=6)  # Create a pool of processes

    pool.starmap(process_spectrum_left, [
            (i, j, path_out, wave, flux_cube[:, i, j],
             np.sqrt(var_cube[:, i, j]), z, ra, dec, path_ex)
            for i in range(0, flux_cube.shape[1])
        ])
    pool.close()
    pool.join()  # Wait for the worker processes to exit

    end_time_j = time.time()  # Record the end time of the outer loop
    # Calculate the elapsed time for the outer loop
    elapsed_time_j = end_time_j - start_time_j
    # Convert elapsed time to minutes and seconds
    minutes, seconds = divmod(elapsed_time_j, 60)
    print(f"j = {j} 的拟合时间 = {int(minutes)}m{seconds:.2f}s")
    return elapsed_time_j


def process_down(i: int, path_out: str, wave: np.ndarray, flux_cube: np.ndarray, var_cube: np.ndarray, z: float, ra: float, dec: float, path_ex: str):
    """
    Fit the spectrum for a given row and generate the fitting results.

    This function fits the spectrum of row `i` using the fitting results
    of the neighboring rows in the specified direction.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    path_out : str
        The directory path where the output results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux_cube : np.ndarray
        The flux cube containing the spectral data.
    var_cube : np.ndarray
        The variance cube of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        The function does not return any value but saves the results to the specified output path.

    Examples
    --------
    >>> process_down(5, './output/', wave_array, flux_cube, var_cube, 0.5, 150.7, -1.5, './path_to_fits/')
    """

    start_time_i = time.time()  # Record the start time of the outer loop
    print(f"开始i = {i}行 ")
    pool = mp.Pool(processes=6)  # Create a pool of processes

    pool.starmap(process_spectrum_down, [
            (i, j, path_out, wave, flux_cube[:, i, j],
             np.sqrt(var_cube[:, i, j]), z, ra, dec, path_ex)
            for j in range(0, flux_cube.shape[2])
        ])
    pool.close()
    pool.join()  # Wait for the worker processes to exit

    end_time_i = time.time()  # Record the end time of the outer loop
    # Calculate the elapsed time for the outer loop
    elapsed_time_i = end_time_i - start_time_i
    # Convert elapsed time to minutes and seconds
    minutes, seconds = divmod(elapsed_time_i, 60)
    print(f"i = {i} 的拟合时间 = {int(minutes)}m{seconds:.2f}s")
    return elapsed_time_i


def process_up(i: int, path_out: str, wave: np.ndarray, flux_cube: np.ndarray, var_cube: np.ndarray, z: float, ra: float, dec: float, path_ex: str):
    """
    Fit the spectrum for a given row and generate the fitting results.

    This function fits the spectrum of row `i` using the fitting results
    of the neighboring rows in the specified direction.

    Parameters
    ----------
    i : int
        The row index for the current spectrum.
    path_out : str
        The directory path where the output results will be stored.
    wave : np.ndarray
        The wavelength array of the spectral data.
    flux_cube : np.ndarray
        The flux cube containing the spectral data.
    var_cube : np.ndarray
        The variance cube of the spectral data.
    z : float
        The redshift value to apply to the spectral data.
    ra : float
        The right ascension coordinate of the target object.
    dec : float
        The declination coordinate of the target object.
    path_ex : str
        The path to the external FITS file.

    Returns
    -------
    None
        The function does not return any value but saves the results to the specified output path.

    Examples
    --------
    >>> process_up(5, './output/', wave_array, flux_cube, var_cube, 0.5, 150.7, -1.5, './path_to_fits/')
    """

    start_time_i = time.time()  # Record the start time of the outer loop
    print(f"开始i = {i}行 ")
    pool = mp.Pool(processes=6)  # Create a pool of processes

    pool.starmap(process_spectrum_up, [
            (i, j, path_out, wave, flux_cube[:, i, j],
             np.sqrt(var_cube[:, i, j]), z, ra, dec, path_ex)
            for j in range(0, flux_cube.shape[2])
        ])
    pool.close()
    pool.join()  # Wait for the worker processes to exit

    end_time_i = time.time()  # Record the end time of the outer loop
    # Calculate the elapsed time for the outer loop
    elapsed_time_i = end_time_i - start_time_i
    # Convert elapsed time to minutes and seconds
    minutes, seconds = divmod(elapsed_time_i, 60)
    print(f"i = {i} 的拟合时间 = {int(minutes)}m{seconds:.2f}s")
    return elapsed_time_i


def process_i_refit(i_1: int, i_2: int, fits_file: str, z: float, scale_factor: int = 1, flux_cube_path: str = None, var_cube_path: str = None, format: str = 'muse'):
    """
    Fit all spectra using the fitting results of adjacent rows.

    This function applies the fitting results of neighboring rows to refit all spectra.
    The spectra in the range from row `i_1` to row `i_2` are fitted twice.

    Parameters
    ----------
    i_1 : int
        The starting row index for the first fitting.
    i_2 : int
        The starting row index for the second fitting.
    fits_file : str
        The path to the input fits file.
    z : float
        The redshift of the target object.
    scale_factor : int, optional
        The scale factor for rebinning the IFU spectra. The default is 1.
    flux_cube_path : str, optional
        The file path to the flux cube of the IFU data. The default is None.
    var_cube_path : str, optional
        The file path to the var cube of the IFU data. The default is None.
    format : str, optional
        The format of the input fits file. The default is 'muse'.

    Returns
    -------
    str
        The path to the output directory where the results is saved.

    Examples
    --------
    >>> process_i_refit(80, 90, 2, 'mmt_1.fits', 0.02, 2, 'reduced2_flux.npy', 'reduced2_var.npy, 'muse')
    """
    # 检查文件是否存在
    if not os.path.exists(fits_file):
        raise FileNotFoundError(f"Error: FITS file '{fits_file}' does not exist.")

    if scale_factor != 1:
        if not os.path.exists(flux_cube_path):
            raise FileNotFoundError(f"Error: Flux cube file '{flux_cube_path}' does not exist.")

        if not os.path.exists(var_cube_path):
            raise FileNotFoundError(f"Error: Variance cube file '{var_cube_path}' does not exist.")
    # 读取数据
    if format == 'muse':
        hdu = fits.open(fits_file)
        hdu.info()

        ra = hdu[0].header['RA']
        dec = hdu[0].header['DEC']
        header = hdu[1].header
        nx, ny, nz = header['NAXIS1'], header['NAXIS2'], header['NAXIS3']
        print(f"IFU cube的shape:nx = {nx}, ny = {ny}, nz = {nz}")

        # 计算波长数组
        wave = np.array(
            header['CRVAL3'] + header['CD3_3'] * np.arange(header['NAXIS3'])
        )
        print('波长：', wave[0], wave[-1], wave.shape)
        print(f'flux的单位:unit = {header["BUNIT"]}')
        flux_in = hdu[1].data
        var_in = hdu[2].data
    if format == 'csst':
        hdu = fits.open(fits_file)
        hdu.info()

        ra = hdu[0].header['RA_cen']
        dec = hdu[0].header['DEC_cen']
        header = hdu[1].header
        nx, ny, nz = header['NAXIS1'], header['NAXIS2'], header['NAXIS3']
        print(f"IFU cube的shape:nx = {nx}, ny = {ny}, nz = {nz}")

        # 计算波长数组
        wave = np.array(hdu[0].header['CRVAL3'] + hdu[0].header['PC3_3']*np.arange(hdu[1].header['NAXIS3']))
        print('波长：', wave[0], wave[-1], wave.shape)
        print(f'flux的单位:unit = {hdu[0].header["F_UNIT"]}')
        flux_in = hdu[1].data
        var_in = 1/hdu[2].data

    path_ex = '.'

    if scale_factor == 1:
        flux_cube = flux_in
        var_cube = var_in
    else:
        flux_cube = np.load(flux_cube_path)
        var_cube = np.load(var_cube_path)
    print(f'flux_cube.shape = {flux_cube.shape}')
    print(f'var_cube.shape = {var_cube.shape}')

    path_out = f'./reduced{scale_factor}_i{i_1}_{i_2}'
    # 检查文件夹是否存在，不存在则创建
    os.makedirs(path_out, exist_ok=True)

    if i_1 < i_2:
        process_i(i_1, path_out, wave, flux_cube, var_cube, z, ra, dec,
                  path_ex)
        for i in range(i_1+1, flux_cube.shape[1]):
            process_down(i, path_out, wave, flux_cube, var_cube, z, ra, dec,
                         path_ex)
        for i in range(i_2, -1, -1):
            process_up(i, path_out, wave, flux_cube, var_cube, z, ra, dec,
                       path_ex)

    else:
        process_i(i_1, path_out, wave, flux_cube, var_cube, z, ra, dec,
                  path_ex)
        for i in range(i_1-1, -1, -1):
            process_up(i, path_out, wave, flux_cube, var_cube, z, ra, dec,
                       path_ex)

        for i in range(i_2, flux_cube.shape[1]):
            process_down(i, path_out, wave, flux_cube, var_cube, z, ra, dec,
                         path_ex)
    return path_out


def process_j_refit(j_1: int, j_2: int, fits_file: str, z: float, scale_factor: int = 1, flux_cube_path: str = None, var_cube_path: str = None, format: str = 'muse'):
    """
    Fit all spectra using the fitting results of adjacent columns.

    This function applies the fitting results of neighboring columns to refit all spectra.
    The spectra in the range from column `j_1` to column `j_2` are fitted twice.

    Parameters
    ----------
    j_1 : int
        The starting column index for the first fitting.
    j_2 : int
        The starting column index for the second fitting.
    fits_file : str
        The path to the input fits file.
    z : float
        The redshift of the target object.
    scale_factor : int, optional
        The scale factor for rebinning the IFU spectra. The default is 1.
    flux_cube_path : str, optional
        The file path to the flux cube of the IFU data. The default is None.
    var_cube_path : str, optional
        The file path to the var cube of the IFU data. The default is None.
    format : str, optional
        The format of the input fits file. The default is 'muse'.

    Returns
    -------
    str
        The path to the output directory where the results is saved.

    Examples
    --------
    >>> process_i_refit(80, 90, 2, 'mmt_1.fits', 0.02, 2, 'reduced2_flux.npy', 'reduced2_var.npy, 'csst')
    """
    # 检查文件是否存在
    if not os.path.exists(fits_file):
        raise FileNotFoundError(f"Error: FITS file '{fits_file}' does not exist.")

    if scale_factor != 1:
        if not os.path.exists(flux_cube_path):
            raise FileNotFoundError(f"Error: Flux cube file '{flux_cube_path}' does not exist.")

        if not os.path.exists(var_cube_path):
            raise FileNotFoundError(f"Error: Variance cube file '{var_cube_path}' does not exist.")
    # 读取数据
    if format == 'muse':
        hdu = fits.open(fits_file)
        hdu.info()

        ra = hdu[0].header['RA']
        dec = hdu[0].header['DEC']
        header = hdu[1].header
        nx, ny, nz = header['NAXIS1'], header['NAXIS2'], header['NAXIS3']
        print(f"IFU cube的shape:nx = {nx}, ny = {ny}, nz = {nz}")

        # 计算波长数组
        wave = np.array(
            header['CRVAL3'] + header['CD3_3'] * np.arange(header['NAXIS3'])
        )
        print('波长：', wave[0], wave[-1], wave.shape)
        print(f'flux的单位:unit = {header["BUNIT"]}')
        flux_in = hdu[1].data
        var_in = hdu[2].data
    if format == 'csst':
        hdu = fits.open(fits_file)
        hdu.info()

        ra = hdu[0].header['RA_cen']
        dec = hdu[0].header['DEC_cen']
        header = hdu[1].header
        nx, ny, nz = header['NAXIS1'], header['NAXIS2'], header['NAXIS3']
        print(f"IFU cube的shape:nx = {nx}, ny = {ny}, nz = {nz}")

        # 计算波长数组
        wave = np.array(hdu[0].header['CRVAL3'] + hdu[0].header['PC3_3']*np.arange(hdu[1].header['NAXIS3']))
        print('波长：', wave[0], wave[-1], wave.shape)
        print(f'flux的单位:unit = {hdu[0].header["F_UNIT"]}')
        flux_in = hdu[1].data
        var_in = 1/hdu[2].data

    path_ex = '.'

    if scale_factor == 1:
        flux_cube = flux_in
        var_cube = var_in
    else:
        flux_cube = np.load(flux_cube_path)
        var_cube = np.load(var_cube_path)
    print(f'flux_cube.shape = {flux_cube.shape}')
    print(f'var_cube.shape = {var_cube.shape}')

    path_out = f'./reduced{scale_factor}_j{j_1}_{j_2}'
    # 检查文件夹是否存在，不存在则创建
    os.makedirs(path_out, exist_ok=True)

    if j_1 < j_2:
        process_j(j_1, path_out, wave, flux_cube, var_cube, z, ra, dec,
                  path_ex)
        for j in range(j_1+1, flux_cube.shape[2]):
            process_right(j, path_out, wave, flux_cube, var_cube, z, ra, dec,
                          path_ex)

        for j in range(j_2, -1, -1):
            process_left(j, path_out, wave, flux_cube, var_cube, z, ra, dec,
                         path_ex)

    else:
        process_j(j_1, path_out, wave, flux_cube, var_cube, z, ra, dec,
                  path_ex)
        for j in range(j_1-1, -1, -1):
            process_left(j, path_out, wave, flux_cube, var_cube, z, ra, dec,
                         path_ex)
        for j in range(j_2, flux_cube.shape[2]):
            process_right(j, path_out, wave, flux_cube, var_cube, z, ra, dec,
                          path_ex)
    return path_out


# def fill_nan(path_out=f'./reduced{scale_factor}_selected_filled/',
#              path_from=f'./reduced{scale_factor}_selected/'):
#     """
#     处理并填充缺失值（NaN）的数据。

#     从源目录读取fits文件，复制文件到目标目录，
#     并处理含有NaN值的光谱数据。
#     """
#     if os.path.exists(path_out):
#         shutil.rmtree(path_out)

#     os.makedirs(path_out, exist_ok=True)
#     ha_vel = np.full((flux_cube.shape[1], flux_cube.shape[2]), np.nan)

#     for i in range(0, flux_cube.shape[1]):
#         for j in range(0, flux_cube.shape[2]):

#             try:
#                 ha_vel[i, j] = fits.open(
#                     path_from+str(i)+'_'+str(j)+'.fits'
#                                          )[1].data['Ha_na_1_centerwave']
#                 file_name_fits = f'{i}_{j}.fits'
#                 file_name_pdf = f'{i}_{j}.pdf'

#                 # 源文件路径
#                 source_path_fits = os.path.join(path_from, file_name_fits)
#                 source_path_pdf = os.path.join(path_from, file_name_pdf)

#                 # 目标文件路径
#                 destination_path_fits = os.path.join(
#                     path_out, file_name_fits)
#                 destination_path_pdf = os.path.join(path_out, file_name_pdf)

#                 # 尝试复制 .fits 文件
#                 try:
#                     shutil.copy(source_path_fits, destination_path_fits)

#                     print(f"复制文件 {file_name_fits}  成功！")

#                 except FileNotFoundError:
#                     print(f"文件 {file_name_fits} 中未找到，跳过该文件。")

#                 print(f"当前循环进度：i = {i}, j = {j}", end='\r')
#                 # 尝试复制 .pdf 文件
#                 try:
#                     shutil.copy(source_path_pdf, destination_path_pdf)
#                     print(f"复制文件 {file_name_pdf} 成功！")
#                 except FileNotFoundError:
#                     print(f"文件 {file_name_pdf} 中未找到，跳过该文件。")
#             except (IOError, KeyError) as e:  # 明确指定要捕获的异常类型
#                 print(f"处理文件 {i},{j} 时出错: {e}")  # 打印错误信息以便调试
#                 continue

#     nan_positions = np.argwhere(np.isnan(ha_vel))

#     pool = mp.Pool(processes=4)
#     args = [(i, j, path_out, path_from) for i, j in nan_positions]

#     pool.starmap(process_spectrum_median, args)

#     pool.close()
#     pool.join()

#     print("\n处理完成!")
