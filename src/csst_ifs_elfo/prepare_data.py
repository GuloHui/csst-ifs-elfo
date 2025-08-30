"""
Identifier:     csst_ifs_elfo/prepare_data.py
Name:           prepare_data
Description:    合并IFU光谱
Author:         Hui Guo
Created:        2025-02-08
Modified-History:
    2025-02-08, created
"""
from astropy.io import fits
import numpy as np


def combine_pixels(ifu_data: np.ndarray, scale_factor: int):
    """
    Rebin the IFU data by combining pixels based on the given scale factor.

    This function takes the original IFU data and re-bins it by combining the pixels
    in blocks of size `scale_factor x scale_factor`. The resulting data has reduced
    height and width, and the pixel values are the sum of the corresponding original pixels.

    Parameters
    ----------
    ifu_data : np.ndarray
        The original IFU data with shape (spectra, height, width).
    scale_factor : int
        The factor by which to reduce the resolution of the IFU data by combining `scale_factor x scale_factor` blocks.

    Returns
    -------
    np.ndarray
        The re-binned IFU data with reduced height and width.

    Examples
    --------
    >>> new_ifu_data = combine_pixels(ifu_data, 2)
    """
    spectra, height, width = ifu_data.shape
    new_height = height // scale_factor
    new_width = width // scale_factor
    print(f"new_height = {new_height}, new_width = {new_width}")
    new_ifu_data = np.full((spectra, new_height, new_width), np.nan)
    for i in range(new_height):
        for j in range(new_width):
            new_ifu_data[:, i, j] = np.sum(
                ifu_data[:, i*scale_factor:(i+1)*scale_factor,
                         j*scale_factor:(j+1)*scale_factor], axis=(1, 2))
    return new_ifu_data


# if __name__ == '__main__':
    # scale_factor = 2
    # hdu = fits.open('mr2251.fits')
    # hdu.info()

    # ra = hdu[0].header['RA']
    # dec = hdu[0].header['DEC']
    # header = hdu[1].header
    # nx, ny, nz = header['NAXIS1'], header['NAXIS2'], header['NAXIS3']
    # print(f"IFU cube的shape:nx = {nx}, ny = {ny}, nz = {nz}")

    # # 计算波长数组
    # wave = np.array(header['CRVAL3'] + header['CD3_3']*np.arange(header['NAXIS3']))
    # print('波长：', wave[0], wave[-1], wave.shape)
    # print(f'flux的单位:unit = {header["BUNIT"]}')

    # flux = hdu[1].data
    # var = hdu[2].data
    # scale_factor = 2
    # # 计算合并后的 flux 和 var
    # reduced_flux = combine_pixels(flux, scale_factor)
    # reduced_var = combine_pixels(var, scale_factor)

    # # 保存合并后的 flux 和 var，使用 f-string 创建文件名
    # np.save(f'reduced{scale_factor}_flux.npy', reduced_flux)
    # np.save(f'reduced{scale_factor}_var.npy', reduced_var)
