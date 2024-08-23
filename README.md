# Joint Decorrelation for Phase-locked Artifact Attenuation

![GitHub](https://img.shields.io/github/license/atr-iine/joint_decorrelation_for_phase-locked_artifact_attenuation)
[![DOI](https://img.shields.io/badge/doi-10.1162/imag_a_00272-informational)](https://doi.org/10.1162/imag_a_00272)


This repository provides a reference implementations of the joint decorrelation (JD) [^1] method to attenuate phase-locked artifacts in neuroimaging data, as utilized in the accompanying research paper [^2].


[^1]: A. de Cheveigné and L. C. Parra, *Joint decorrelation, a versatile tool for multichannel data analysis,* **NeuroImage**, vol. 98, pp. 487–505, 2014, doi: [10.1016/j.neuroimage.2014.05.068](https://doi.org/10.1016/j.neuroimage.2014.05.068).

[^2]: T. Kuroda, R. J. Kobler, T. Ogawa, M. Tsutsumi, T. Kishi, and M. Kawanabe, *Test-retest reliability of EEG microstate metrics for evaluating noise reductions in simultaneous EEG-fMRI,* **Imaging Neuroscience**, 2024, doi: [10.1162/imag_a_00272](https://doi.org/10.1162/imag_a_00272).


## File List
The following files are provided in this repository:

- `demo.ipynb`: Jupyter notebook that demonstrates the Python reference implementation to use joint decorrelation to attenuate eye blink and ECG artifacts in EEG data.
- `jointdecorrelation.py`: Python module that implements the joint decorrelation algorithm.

- `jointdecorrelation.m`: Matlab function that implements the joint decorrelation algorithm.

## Usage

### Local Installation
If you want to run it locally on your machine, Python3 and Jupyter are needed.
The present code was developed and tested with the following packages:
```
- Python >= 3.8
- numpy
- scipy
- jupyter
- pyriemann
```

Make sure you have [Python3](https://www.python.org/downloads/) installed on
your computer.
You can then install the required packages (including Jupyter) by running
```bash
pip install -r requirements.txt
```

Finally, you can start a Jupyter session with
```bash
jupyter notebook
```

and open and run the `demo.ipynb` notebook.


## Acknowledgements
This work was supported by Innovative Science and Technology Initiative for
Security Grant Number JPJ004596, ATLA, Japan

This README file is based on the [template](https://github.com/klb2/reproducible-paper-python-template) by @klb2.

## License and Referencing
This program is licensed under the MIT license. If you in any way use this
code for research that results in publications, please cite the original article introducing joint decorrelation and our original article listed above.

You can use the following BibTeX entry
```bibtex

@article{de_cheveigne_joint_2014,
	title = {Joint decorrelation, a versatile tool for multichannel data analysis},
	volume = {98},
	doi = {10.1016/j.neuroimage.2014.05.068},
	journal = {NeuroImage},
	author = {de Cheveigné, Alain and Parra, Lucas C.},
	year = {2014},
	pages = {487--505}
}

@article{kuroda_test-retest_2024,
	title = {Test-retest reliability of {EEG} microstate metrics for evaluating noise reductions in simultaneous {EEG}-{fMRI}},
	doi = {10.1162/imag_a_00272},
	journal = {Imaging Neuroscience},
	author = {Kuroda, Toshikazu and Kobler, Reinmar J. and Ogawa, Takeshi and Tsutsumi, Mizuki and Kishi, Tomohiko and Kawanabe, Motoaki},
	year = {2024}
}
```