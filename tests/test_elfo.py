"""
Identifier:     tests/test_elfo.py
Name:           test_elfo.py
Description:    test for elfo
Author:         Hui Guo
Created:        2025-02-08
Modified-History:
    2025-02-08, created
"""
import unittest
import numpy as np
import os
import shutil
from csst_ifs_elfo import create_qsopar
from csst_ifs_elfo import combine_pixels
from csst_ifs_elfo import process_i_refit
from csst_ifs_elfo import process_j_refit
from csst_ifs_elfo import process_down
from csst_ifs_elfo import process_spectrum_down
from csst_ifs_elfo import downward
from csst_ifs_elfo import process_up
from csst_ifs_elfo import process_spectrum_up
from csst_ifs_elfo import upward
from csst_ifs_elfo import process_left
from csst_ifs_elfo import process_spectrum_left
from csst_ifs_elfo import leftward
from csst_ifs_elfo import process_right
from csst_ifs_elfo import process_spectrum_right
from csst_ifs_elfo import rightward
from csst_ifs_elfo import process_i
from csst_ifs_elfo import process_j
from csst_ifs_elfo import process_spectrum


class Test_elfo(unittest.TestCase):

    def setUp(self):
        # 准备测试数据
        dir_parcube = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/parcube.npy")
        self.parcube = np.load(dir_parcube)  # 加载 parcube.npy 数据
        self.i_1 = 20
        self.i_2 = 30
        self.j_1 = 30
        self.j_2 = 40
        self.test_dir = f'./reduced4_i{self.i_1}_{self.i_2}'
        self.test_dir_2 = f'./reduced4_j{self.j_1}_{self.j_2}'
        # 将测试目录设为类属性，这样类方法才能访问
        Test_elfo.test_dir = self.test_dir
        Test_elfo.test_dir_2 = self.test_dir_2

    def test_combine_pixels(self):
        # 模拟IFU数据
        spectra, height, width = 3, 6, 12  # 示例数据大小
        ifu_data = np.random.random((spectra, height, width))  # 模拟flux数据
        scale_factor = 3

        # 使用combine_pixels函数进行合并
        reduced_data = combine_pixels(ifu_data, scale_factor)

        # 保存合并后的数据
        np.save(f'reduced{scale_factor}_flux.npy', reduced_data)

        # 检查文件是否生成
        filename = f'reduced{scale_factor}_flux.npy'
        print(f"Checking if {filename} was created...")
        self.assertTrue(os.path.exists(filename), f"File {filename} was not created!")

    def test_create_qsopar(self):
        # 测试文件是否生成
        filename = 'qsopar_0_0.fits'

        # 在执行之前确保文件不存在
        if os.path.exists(filename):
            os.remove(filename)

        # 调用函数
        create_qsopar(000, 000, self.parcube)

        # 检查文件是否创建
        self.assertTrue(os.path.exists(filename), f"File {filename} was not created!")

    def test_process_i_refit(self):
        """
        测试 process_i_refit 是否创建文件夹，并且文件夹内文件数大于500。
        """
        # 调用 process_i_refit 函数
        dir_mr2251 = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/mr2251.fits")
        dir_csst = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/CSST_IFS_CUBE_104636.8+631325_R8202_LIN_60671.64611.fits")
        flux_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_flux.npy")
        var_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_var.npy")

        path_out = process_i_refit(i_1=self.i_1, i_2=self.i_2, fits_file=dir_mr2251, z=0.064, scale_factor=4, flux_cube_path=flux_cube_path, var_cube_path=var_cube_path, format='muse')

        # 检查文件夹是否存在
        self.assertTrue(os.path.exists(path_out), f"Directory {path_out} was not created!")

        # 检查文件夹内的文件数量
        files = os.listdir(path_out)
        file_count = len(files)

        # 判断文件数是否大于500
        self.assertGreater(file_count, 500, f"Expected more than 500 files, but got {file_count}.")

        path_out = process_i_refit(i_1=31, i_2=21, fits_file=dir_csst, z=0.003373, format='csst')
        # 检查文件夹是否存在
        self.assertTrue(os.path.exists(path_out), f"Directory {path_out} was not created!")

        # 检查文件夹内的文件数量
        files = os.listdir(path_out)
        file_count = len(files)

        # 判断文件数是否大于500
        self.assertGreater(file_count, 500, f"Expected more than 500 files, but got {file_count}.")

    def test_process_j_refit(self):
        """
        测试 process_j_refit 是否创建文件夹，并且文件夹内文件数大于500。
        """
        # 调用 process_j_refit 函数
        dir_mr2251 = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/mr2251.fits")
        dir_csst = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/CSST_IFS_CUBE_104636.8+631325_R8202_LIN_60671.64611.fits")
        flux_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_flux.npy")
        var_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_var.npy")

        path_out = process_j_refit(j_1=self.j_1, j_2=self.j_2, fits_file=dir_mr2251, z=0.064, scale_factor=4, flux_cube_path=flux_cube_path, var_cube_path=var_cube_path, format='muse')

        # 检查文件夹是否存在
        self.assertTrue(os.path.exists(path_out), f"Directory {path_out} was not created!")

        # 检查文件夹内的文件数量
        files = os.listdir(path_out)
        file_count = len(files)

        # 判断文件数是否大于500
        self.assertGreater(file_count, 500, f"Expected more than 500 files, but got {file_count}.")

        path_out = process_j_refit(j_1=31, j_2=21, fits_file=dir_csst, z=0.003373, format='csst')
        # 检查文件夹是否存在
        self.assertTrue(os.path.exists(path_out), f"Directory {path_out} was not created!")

        # 检查文件夹内的文件数量
        files = os.listdir(path_out)
        file_count = len(files)

        # 判断文件数是否大于500
        self.assertGreater(file_count, 500, f"Expected more than 500 files, but got {file_count}.")

    def test_process_refit_missing_files(self):
        """测试缺失 FITS、flux_cube 和 var_cube 时是否抛出 FileNotFoundError"""

        # 1. 测试 FITS 文件不存在
        with self.assertRaises(FileNotFoundError) as cm:
            process_i_refit(80, 90, 'missing.fits', 0.02)
        self.assertIn("Error: FITS file 'missing.fits' does not exist.", str(cm.exception))

        # 创建一个假 FITS 文件
        fake_fits = 'fake.fits'
        with open(fake_fits, 'w') as f:
            f.write('')

        # 2. 测试 flux_cube 缺失
        with self.assertRaises(FileNotFoundError) as cm:
            process_i_refit(80, 90, fake_fits, 0.02, scale_factor=2, flux_cube_path='missing_flux.npy', var_cube_path='var.npy')
        self.assertIn("Error: Flux cube file 'missing_flux.npy' does not exist.", str(cm.exception))

        # 创建一个假 flux.npy
        fake_flux = 'fakeflux.npy'
        with open(fake_flux, 'wb') as f:
            f.write(b'')

        # 3. 测试 var_cube 缺失
        with self.assertRaises(FileNotFoundError) as cm:
            process_i_refit(80, 90, fake_fits, 0.02, scale_factor=2, flux_cube_path=fake_flux, var_cube_path='missing_var.npy')
        self.assertIn("Error: Variance cube file 'missing_var.npy' does not exist.", str(cm.exception))

        # 清理测试文件
        os.remove(fake_fits)
        os.remove(fake_flux)

        # 1. 测试 FITS 文件不存在
        with self.assertRaises(FileNotFoundError) as cm:
            process_j_refit(80, 90, 'missing.fits', 0.02)
        self.assertIn("Error: FITS file 'missing.fits' does not exist.", str(cm.exception))

        # 创建一个假 FITS 文件
        fake_fits = 'fake.fits'
        with open(fake_fits, 'w') as f:
            f.write('')

        # 2. 测试 flux_cube 缺失
        with self.assertRaises(FileNotFoundError) as cm:
            process_j_refit(80, 90, fake_fits, 0.02, scale_factor=2, flux_cube_path='missing_flux.npy', var_cube_path='var.npy')
        self.assertIn("Error: Flux cube file 'missing_flux.npy' does not exist.", str(cm.exception))

        # 创建一个假 flux.npy
        fake_flux = 'fakeflux.npy'
        with open(fake_flux, 'wb') as f:
            f.write(b'')

        # 3. 测试 var_cube 缺失
        with self.assertRaises(FileNotFoundError) as cm:
            process_j_refit(80, 90, fake_fits, 0.02, scale_factor=2, flux_cube_path=fake_flux, var_cube_path='missing_var.npy')
        self.assertIn("Error: Variance cube file 'missing_var.npy' does not exist.", str(cm.exception))

        # 清理测试文件
        os.remove(fake_fits)
        os.remove(fake_flux)

    def test_z_process_i(self):
        flux_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_flux.npy")
        var_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_var.npy")

        flux_cube = np.load(flux_cube_path)
        var_cube = np.load(var_cube_path)
        seconds = process_i(self.i_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube, var_cube, 0.064, 343.524658, -17.58185, '.')
        self.assertGreater(seconds, 20, f"Expected processing time to be greater than 20 seconds, but got {seconds} s.")

    def test_z_process_j(self):
        flux_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_flux.npy")
        var_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_var.npy")

        flux_cube = np.load(flux_cube_path)
        var_cube = np.load(var_cube_path)
        seconds = process_j(self.j_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube, var_cube, 0.064, 343.524658, -17.58185, '.')
        self.assertGreater(seconds, 20, f"Expected processing time to be greater than 20 seconds, but got {seconds} s.")

    def test_z_process_spectrum(self):
        flux_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_flux.npy")
        var_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_var.npy")

        flux_cube = np.load(flux_cube_path)
        var_cube = np.load(var_cube_path)
        result = process_spectrum(self.i_1, self.j_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube[:, self.i_1, self.j_1], var_cube[:, self.i_1, self.j_1], 0.064, 343.524658, -17.58185, '.')
        self.assertEqual(result, 1, f"Expected 1, but got {result}")

    def test_z_process_row(self):
        flux_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_flux.npy")
        var_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_var.npy")

        flux_cube = np.load(flux_cube_path)
        var_cube = np.load(var_cube_path)
        seconds = process_down(self.i_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube, var_cube, 0.064, 343.524658, -17.58185, '.')
        self.assertGreater(seconds, 20, f"Expected processing time to be greater than 20 seconds, but got {seconds} s.")
        seconds = process_up(self.i_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube, var_cube, 0.064, 343.524658, -17.58185, '.')
        self.assertGreater(seconds, 20, f"Expected processing time to be greater than 20 seconds, but got {seconds} s.")

    def test_z_process_spectrum_row(self):
        flux_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_flux.npy")
        var_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_var.npy")

        flux_cube = np.load(flux_cube_path)
        var_cube = np.load(var_cube_path)
        result = process_spectrum_down(self.i_1, self.j_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube[:, self.i_1, self.j_1], var_cube[:, self.i_1, self.j_1], 0.064, 343.524658, -17.58185, '.')
        self.assertEqual(result, 1, f"Expected 1, but got {result}")
        result = process_spectrum_up(self.i_1, self.j_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube[:, self.i_1, self.j_1], var_cube[:, self.i_1, self.j_1], 0.064, 343.524658, -17.58185, '.')
        self.assertEqual(result, 1, f"Expected 1, but got {result}")

    def test_z_row(self):
        chi2, par_cube_choosen, coords = downward(self.i_1+1, self.j_1, self.j_1, self.j_1+1, self.test_dir)
        self.assertIsInstance(chi2[0], float)
        self.assertIsInstance(par_cube_choosen, np.ndarray)
        self.assertEqual(tuple(coords), (self.i_1, self.j_1), f"Expected coords to be ({self.i_1}, {self.j_1}), but got {coords}")
        chi2, par_cube_choosen, coords = upward(self.i_1-1, self.j_1, self.j_1, self.j_1+1, self.test_dir)
        self.assertIsInstance(chi2[0], float)
        self.assertIsInstance(par_cube_choosen, np.ndarray)
        self.assertEqual(tuple(coords), (self.i_1, self.j_1), f"Expected coords to be ({self.i_1}, {self.j_1}), but got {coords}")

    def test_z_process_column(self):
        flux_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_flux.npy")
        var_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_var.npy")

        flux_cube = np.load(flux_cube_path)
        var_cube = np.load(var_cube_path)
        seconds = process_left(self.j_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube, var_cube, 0.064, 343.524658, -17.58185, '.')
        self.assertGreater(seconds, 20, f"Expected processing time to be greater than 20 seconds, but got {seconds} s.")
        seconds = process_right(self.j_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube, var_cube, 0.064, 343.524658, -17.58185, '.')
        self.assertGreater(seconds, 20, f"Expected processing time to be greater than 20 seconds, but got {seconds} s.")

    def test_z_process_spectrum_column(self):
        flux_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_flux.npy")
        var_cube_path = os.path.join(os.environ["UNIT_TEST_DATA_ROOT"], "csst_ifs/csst_ifs_elfo/reduced4_var.npy")

        flux_cube = np.load(flux_cube_path)
        var_cube = np.load(var_cube_path)
        result = process_spectrum_left(self.i_1, self.j_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube[:, self.i_1, self.j_1], var_cube[:, self.i_1, self.j_1], 0.064, 343.524658, -17.58185, '.')
        self.assertEqual(result, 1, f"Expected 1, but got {result}")
        result = process_spectrum_right(self.i_1, self.j_1, self.test_dir, np.array(4750.17578125 + 1.25*np.arange(3681)), flux_cube[:, self.i_1, self.j_1], var_cube[:, self.i_1, self.j_1], 0.064, 343.524658, -17.58185, '.')
        self.assertEqual(result, 1, f"Expected 1, but got {result}")

    def test_z_column(self):
        chi2, par_cube_choosen, coords = leftward(self.i_1, self.j_1-1, self.i_1, self.i_1+1, self.test_dir)
        self.assertIsInstance(chi2[0], float)
        self.assertIsInstance(par_cube_choosen, np.ndarray)
        self.assertEqual(tuple(coords), (self.i_1, self.j_1), f"Expected coords to be ({self.i_1}, {self.j_1}), but got {coords}")
        chi2, par_cube_choosen, coords = rightward(self.i_1, self.j_1+1, self.i_1, self.i_1+1, self.test_dir)
        self.assertIsInstance(chi2[0], float)
        self.assertIsInstance(par_cube_choosen, np.ndarray)
        self.assertEqual(tuple(coords), (self.i_1, self.j_1), f"Expected coords to be ({self.i_1}, {self.j_1}), but got {coords}")

    @classmethod
    def tearDownClass(cls):
        # 清理生成的文件
        files_to_remove = [
            'reduced3_flux.npy',  # 假设scale_factor是2
            'qsopar_0_0.fits',
        ]

        for filename in files_to_remove:
            if os.path.exists(filename):
                os.remove(filename)
        if os.path.exists(cls.test_dir):
            shutil.rmtree(cls.test_dir)
        if os.path.exists(cls.test_dir_2):
            shutil.rmtree(cls.test_dir_2)


if __name__ == "__main__":
    unittest.main()
