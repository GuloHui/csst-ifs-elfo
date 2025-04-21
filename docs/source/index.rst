.. csst-ifs-elfo documentation master file, created by
   sphinx-quickstart on Fri Apr 18 15:23:34 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. title:: ELFO documentation



.. image:: figures/LOGO.png
    :width: 600px
    :align: center
    :alt: ELFO logo

| 

ELFO Documentation
=================================================


This document describes the principles behind ELFO and explains how to install and use ELFO.

**ELFO** (Emission Line Fitting Optimization) is a Python package designed to improve emission‐line 
fitting in integral‐field spectroscopy (IFS) data taken with the Chinese Space Station Telescope (CSST‑IFS).



----------


**ELFO** uses `PyQSOFit`_ to fit each spectrum individually and leverages the spatial correlation inherent in 
IFU data—stemming from the continuity of physical 
processes across neighboring spaxels—to improve emission‑line fitting.
`PyQSOFit`_ employs the Levenberg–Marquardt least‐squares algorithm, and choosing appropriate initial guesses 
for the model parameters is critical. **ELFO** processes spectra in a user‑defined spatial sequence, 
using the fitted results of neighboring spaxels to set the initial parameters for each fit. 

We also provide a selection algorithm to pick the best results from different fitting orders. 
Currently **ELFO** is used to improve Hα emission‐line fits in quasar spectra, but its general framework makes 
it easy to apply to other quasar lines or to emission‐line fitting in any other IFS data.


|

.. toctree::
   :maxdepth: 2
   :caption: Usage

   installation
   


.. toctree::
   :maxdepth: 2
   :caption: Methods

   workflow

.. toctree::
   :maxdepth: 2
   :caption: API/Code Reference

   api/elfo.rst
   api/para.rst
   api/prepare_data.rst

|

----

:Author: Hui Guo
:Institute: university of science and technology of china
:Contact: guohui@mail.ustc.edu.cn
:Last updated: 2025-04


----


|

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. include:: include/links.rst