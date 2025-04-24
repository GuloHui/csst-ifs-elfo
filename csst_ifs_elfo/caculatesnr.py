"""
Identifier:     csst_ifs_elfo/caculatesnr.py
Name:           caculatesnr.py
Description:    用ha的拟合结果计算信噪比
Author:         Hui Guo
Created:        2025-02-08
Modified-History:
    2025-02-08, created
"""
# 用ha的拟合结果计算信噪比
# 提取Hα相关的所有成分参数

import os
import numpy as np
from PyAstronomy import pyasl
from astropy.io import fits
import multiprocessing as mp
from functools import partial
import time
from sfdmap2 import sfdmap


def _Manygauss(xval, pp):
    """

    xval: wavelength array in AA

    pp: Paramaters array [ngauss, 3]
        scale: line amplitude
        wave: central ln wavelength in AA
        sigma: width in km/s
    """

    return np.sum(
        pp[:, 0] * np.exp(
            -(xval[:, np.newaxis] - pp[:, 1])**2 / (2*pp[:, 2]**2)
        ), axis=1)


def caculate_sn(pp, wave, err):

    pp = np.array(pp).astype(float)
    ngauss = len(pp)//3

    pp_shaped = pp.reshape([len(pp)//3, 3])

    if ngauss == 0:
        snr = 0
    else:
        cen = pp_shaped[:, 1]
        sig = pp_shaped[:, 2]

        # print cen,sig,area
        left = np.min(cen - 3*sig)
        right = np.max(cen + 3*sig)
        disp = 1e-4*np.log(10)
        npix = int((right - left)/disp)

        xx = np.linspace(left, right, npix)
        yy = _Manygauss(xx, pp_shaped)

        ypeak = yy.max()
        ypeak_ind = np.argmax(yy)
        peak_wave = np.exp(xx[ypeak_ind])

        # 找到最接近峰值波长的两个观测点
        peak_idx = np.argmin(np.abs(wave - peak_wave))

        # 确保我们有左右两个点进行插值
        if wave[peak_idx] > peak_wave and peak_idx > 0:
            idx_left = peak_idx - 1
            idx_right = peak_idx
        elif wave[peak_idx] < peak_wave and peak_idx < len(wave) - 1:
            idx_left = peak_idx
            idx_right = peak_idx + 1
        else:
            # 如果在边界，就只用最近的点
            noise = err[peak_idx]
            return ypeak/noise
        # plot_multigauss_fit(wave, line_flux, err,pp_shaped)
        # 线性插值计算误差
        wave_left = wave[idx_left]
        wave_right = wave[idx_right]
        err_left = err[idx_left]
        err_right = err[idx_right]

        # 线性插值公式
        noise = err_left + (peak_wave - wave_left) * (err_right - err_left) / \
            (wave_right - wave_left)

        snr = ypeak/noise

        # print(peak_wave)
        # print(peak_idx)
        # print(ypeak)
        # print(noise)

        return snr


def process_pixel(i, j, a, b, flux_cube, var_cube, wave, header, z, ra, dec,
                  scale_factor, dustmap_path):
    print(f"Processing pixel ({i}, {j})")
    """
    处理单个像素点的函数，用于并行化。
    """
    try:
        flux = flux_cube[:, i, j]
        err = np.sqrt(var_cube[:, i, j])
        m = sfdmap.SFDMap(dustmap_path)
        zero_flux = np.where(flux == 0, True, False)

        flux[zero_flux] = 1e-10
        flux_unred = pyasl.unred(wave, flux, m.ebv(ra, dec))
        err_unred = err*flux_unred/flux
        flux_unred[zero_flux] = 0
        del flux, err
        flux = flux_unred
        err = err_unred
        wave = wave/(1+z)
        flux = flux*(1+z)
        err = err*(1+z)
        data = fits.open(f'./reduced2_selected/{i}_{j}.fits')[1].data
        params = []

        # 提取各成分参数
        components = [
            'Ha_br1', 'Ha_br2', 'Ha_na',
            'NII6549', 'NII6585', 'SII6718', 'SII6732'
        ]
        for comp in components:
            scale = data[f'{comp}_1_scale'][0]
            center = data[f'{comp}_1_centerwave'][0]
            sigma = data[f'{comp}_1_sigma'][0]
            params.extend([scale, center, sigma])

        params = np.array(params)

        # 计算信噪比
        return caculate_sn(params, wave, err)

    except Exception as e:
        print(f"Error processing pixel ({i}, {j}): {e}")

        return np.nan


def parallel_process_snr(a, b, scale_factor, flux_cube, var_cube, wave, header,
                         z, ra, dec, dustmap_path):
    """
    并行处理所有像素点，返回信噪比矩阵。
    """
    hasnr = np.full((a, b), np.nan)

    # 创建进度条

    # 定义并行池
    with mp.Pool(processes=5) as pool:
        # 使用偏函数简化参数传递
        func = partial(
            process_pixel, a=a, b=b, flux_cube=flux_cube, var_cube=var_cube,
            wave=wave, header=header, z=z, ra=ra, dec=dec,
            scale_factor=scale_factor, dustmap_path=dustmap_path
        )

        # 创建任务列表
        tasks = [(i, j) for i in range(0, a) for j in range(0, b)]
        print(len(tasks))  # 检查是否有 33630 个任务
        print(tasks[:5])   # 查看前几个任务的格式是否正确
        # aaa = func(50,50)     # 测试单个任务的格式是否正确
        # print(aaa)

        # 并行执行
        # results = pool.starmap(func, tasks)
        results = []
        for res in pool.starmap(func, tasks):
            results.append(res)

    # 填充结果矩阵
    for idx, (i, j) in enumerate(tasks):
        hasnr[i, j] = results[idx]

    return hasnr


if __name__ == '__main__':
    start = time.time()

    scale_factor = 2
    z = 0.064
    hdu = fits.open('mr2251.fits')
    header = hdu[1].header
    a = hdu[1].data.shape[1] // scale_factor
    b = hdu[1].data.shape[2] // scale_factor
    wave = np.array(
        header['CRVAL3'] + header['CD3_3'] * np.arange(header['NAXIS3'])
    )
    flux_cube = np.load(f'reduced{scale_factor}_flux.npy')
    var_cube = np.load(f'reduced{scale_factor}_var.npy')
    dustmap_dir = './pyqsofit'
    dustmap_path = os.path.join(dustmap_dir, 'sfddata')

    ra = hdu[0].header['RA']
    dec = hdu[0].header['DEC']

    # 并行计算
    hasnr = parallel_process_snr(
        a, b, scale_factor, flux_cube, var_cube, wave, header, z, ra, dec,
        dustmap_path
    )

    # 保存结果
    np.save('./hasnr_selected.npy', hasnr)
    end = time.time()
    elapsed_time = end - start

    elapsed_time_minutes = elapsed_time / 60

    print(f"函数运行时间: {elapsed_time_minutes:.2f} 分钟")
