# csst-ifs-elfo

<div align="center">
  <img src="docs/source/figures/logo.png" alt="img" width="600"/>
</div>

**ELFO** (Emission Line Fitting Optimization) is a Python package for emission line fitting optimization in integral field spectroscopy data, based on [PyQSOFit][PyQSOFit]. It is designed for integral‐field spectroscopy (IFS) data taken with the Chinese Space Station Telescope (CSST‑IFS).


**ELFO** uses the results of neighboring spectra to determine initial guesses and selects the most spatially smooth solutions from multiple fitting attempts. The method has already been validated to CSST’s simulated IFS data and with slight modifications, it can be applied to other IFS data and different emission lines.

## Install

**csst-ifs-elfo** can be installed as follows:

``git clone https://github.com/GuloHui/csst-ifs-elfo``<br>

``cd csst-ifs-elfo``<br>

``python -m pip install .``<br>

## Usage

A brief example of how to use **ELFO** is shown in [example][example] notebook.

**See our full documentation at [csst-ifs-elfo.readthedocs.io](https://csst-ifs-elfo.readthedocs.io).**

A detailed description of the method and its applications has been presented in

[PyQSOFit]: https://github.com/legolason/PyQSOFit
[example]: https://github.com/GuloHui/csst-ifs-elfo/blob/main/example/example.ipynb
