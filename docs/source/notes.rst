.. include:: include/links.rst

.. highlight:: python
   :linenothreshold: 3


Fit All Spectra with Default Parameters
-----------------------------------------

In ELFO, the ``process_i`` and ``process_j`` functions use the default parameter file 
``qsopar_muse.fits`` to fit an entire row of spectra. You can refer to the example file on the `PyQSOFit`_ website 
for guidance on how to create this file. To fit all spectra, you may refer to the code below.


.. code-block:: python

    for i in range(npixel):
        process_i(i, path_out, wave, flux_cube, var_cube, 0.05668, ra, dec, path_ex)
        
        # or

    for j in range(npixel):
        process_j(j, path_out, wave, flux_cube, var_cube, 0.05668, ra, dec, path_ex)



The effect of the code is the same as using a simple `for` loop with `PyQSOFit`_ to fit each pixel. However, 
the function here uses multiprocessing to perform the fitting, which will be more efficient.