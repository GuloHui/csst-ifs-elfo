.. include:: include/links.rst




.. highlight:: python
   :linenothreshold: 3


.. _workflow:




ELFO workflow
==============





This section summarizes the fitting strategy employed in the ``process_i_refit`` and ``process_j_refit`` functions of the ELFO package.


User-defined Parameters
-----------------------

- **Fitting start rows**: ``x0``, ``x1`` (correspond to ``i_1`` and ``i_2`` in ``process_i_refit``)
- **Fitting start columns**: ``y0``, ``y1`` (correspond to ``j_1`` and ``j_2`` in ``process_i_refit``)
- **Neighborhood range**: ``k`` (number of nearby pixels to consider, default is 9.  
  This value is hard-coded for simplicity and not exposed as a function parameter.  
  If you need to modify it, feel free to contact me.)

Fitting Flow
------------

``process_i_refit``
************************

**Fitting procedure when x0 < x1**

``x0`` ---> edge fitting


- Row ``x0`` is fitted with default initial guesses.
- When fitting specific pixels in row ``x0 + 1``, the initial guess for the emission-line fit is taken from the result with the smallest :math:`\chi^2_{\mathrm{red}}` value among the emission-line fits of the ``k`` neighboring pixels in row ``x0``.
- The process repeats until the edge of the data is reached.

``x1``----> ``x = 0`` fitting

- Restarting from row ``x1``
- When fitting specific pixels in row ``x1``, the initial guess for the emission-line fit is taken from the result with the smallest :math:`\chi^2_{\mathrm{red}}` value among the emission-line fits of the ``k`` neighboring pixels in row ``x1 + 1``.
- The process continues by fitting row ``x1 - 1``, and this repeats until reaching row ``0``.



**Fitting procedure when x0 > x1**



``x0`` -----> ``x = 0`` fitting


- Row ``x0`` is fitted with default initial guesses.
- When fitting specific pixels in row ``x0 - 1``, the initial guess for the emission-line fit is taken from the result with the smallest :math:`\chi^2_{\mathrm{red}}` value among the emission-line fits of the ``k`` neighboring pixels in row ``x0``.
- The process repeats until reaching row ``0``.

``x1``----> edge fitting

- Restarting from row ``x1``
- When fitting specific pixels in row ``x1``, the initial guess for the emission-line fit is taken from the result with the smallest :math:`\chi^2_{\mathrm{red}}` value among the emission-line fits of the ``k`` neighboring pixels in row ``x1 - 1``.
- The process continues by fitting row ``x1 + 1``, and this repeats until the edge of the data is reached.


A simple summary of ``process_i_refit`` 's fitting is that it first fits half of the spectra from bottom to top, and then fits the other half from top to bottom. 
The following diagram illustrates the fitting sequence when ``x0`` < ``x1`` — it might help make things clearer.



.. figure:: figures/flow.png
    :height: 200px
    :align: center

    ELFO workflow



``process_j_refit``
************************



**Left-to-right and right-to-left fitting** 

- The process works similarly on the columns.
- Each pixel uses the best nearby fit as its initial guess.

Goal
----

Ensure every pixel’s spectrum is fitted more robustly by leveraging information from spatially neighboring fits.




