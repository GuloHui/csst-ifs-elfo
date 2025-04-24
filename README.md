# csst-ifs-elfo

A program that uses the fitting results of adjacent points to improve the IFU spectral emission line fitting.

<div align="center">
  <img src="docs/source/figures/logo.jpg" alt="img" width="600"/>
</div>



**ELFO** (Emission Line Fitting Optimization) is a Python package designed to improve emission‐line 
fitting in integral‐field spectroscopy (IFS) data taken with the Chinese Space Station Telescope (CSST‑IFS).



-------------

**ELFO** uses [PyQSOFit] to fit each spectrum individually and leverages the spatial correlation inherent in 
IFS data—stemming from the continuity of physical 
processes across neighboring spaxels—to improve emission‑line fitting.
[PyQSOFit] employs the Levenberg–Marquardt least‐squares algorithm, and choosing appropriate initial guesses 
for the model parameters is critical. **ELFO** processes spectra in a user‑defined spatial sequence, 
using the fitted results of neighboring spaxels to set the initial parameters for each fit. 

We also provide a selection algorithm to pick the best results from different fitting orders. 
Currently **ELFO** is used to improve Hα emission‐line fits in quasar spectra, but its general framework makes 
it easy to apply to other quasar lines or to emission‐line fitting in any other IFS data.











**See our full documentation at [csst-ifs-elfo.readthedocs.io](https://csst-ifs-elfo.readthedocs.io).**








----------
[PyQSOFit]:  https://github.com/legolason/PyQSOFit