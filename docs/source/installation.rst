.. include:: include/links.rst

.. highlight:: python
   :linenothreshold: 3

Installation
============

Quick Install
-------------

**csst-ifs-elfo** can be installed as follows::

  git clone https://csst-tb.bao.ac.cn/code/csst-pipeline/ifs/csst_ifs_elfo.git

How to Fit IFS Data
====================

Prepare the IFS data
----------------------
The IFS simulated observation data from CSST can be downloaded from the following link: `Download <http://202.127.29.3/~shen/NGC3359/>`_.

Import Necessary Functions
--------------------------

Import the  functions from **csst-ifs-elfo** to perform the fitting:

.. note:: 
  Both functions share the same fitting logic, fitting the spectra column by column.  
  The only difference lies in the direction of progression  ``process_i_refit`` proceeds along rows,  
  while ``process_j_refit`` proceeds along columns.
   
.. code-block:: python

  from csst_ifs_elfo import process_i_refit
  from csst_ifs_elfo import process_j_refit

The fitting of the IFS data is performed using the two imported functions.  
In the following, we take ``process_i_refit`` as an example, which fits all spectra using the fitting results of adjacent rows.  
For a detailed description of the fitting strategy, see the :ref:`workflow` section below.



**Parameters**

- ``i_1``: The starting row index for the first fitting. Type: ``int``.
- ``i_2``: The starting row index for the second fitting. Type: ``int``.
- ``fits_file``: The path to the input FITS file. Type: ``str``.
- ``z``: Redshift of the target object. Type: ``float``.
- ``scale_factor`` *(optional)*: The scale factor for rebinning the IFU spectra. Default is ``1``. Type: ``int``.
- ``flux_cube_path`` *(optional)*: The file path to the flux cube of the IFU data. Default is ``None``. Type: ``str``.
- ``var_cube_path`` *(optional)*: The file path to the variance cube of the IFU data. Default is ``None``. Type: ``str``.
- ``format`` *(optional)*: The format of the input FITS file. Default is ``'muse'``. Type: ``str``.

**Returns**

- ``str``: The path to the output directory where the results are saved.

**Usage Example**

.. code-block:: python

    # Using required parameters only
    process_i_refit(80, 90, 'mmt_1.fits', 0.02, 'csst')
    # Using optional parameters for custom flux/var cubes and scale factor
    process_i_refit(80, 90, 'mmt_1.fits', 0.02, 2, 'reduced2_flux.npy', 'reduced2_var.npy', 'csst')






