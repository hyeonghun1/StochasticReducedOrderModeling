# StochasticReducedOrderModeling
This repository contains a MATLAB (R2024b) implementation of the Operator Inference (OpInf) framework for stochastic differential equations (SDEs).

The implementation is based on and extends the methodology introduced in:

[1] M. A. Freitag, J. M. Nicolaus, M. Redmann
 	[Learning Stochastic Reduced Models from Data: A Nonintrusive Approach](https://epubs.siam.org/doi/full/10.1137/24M1679756)<details><summary>BibTex</summary><pre>
@article{freitag2025learning,
  title={Learning stochastic reduced models from data: A nonintrusive approach},
  author={Freitag, MA and Nicolaus, JM and Redmann, M},
  journal={SIAM Journal on Scientific Computing},
  volume={47},
  number={5},
  pages={A2851--A2880},
  year={2025},
  publisher={SIAM}
}</pre></details>

# Scope of This Repository
Compared to the original stochastic OpInf implementation in [1], this repository extends the framework to experimental data. In particular, it applies stochastic operator inference to capillary wave turbulence data acquired using a custom ultra–high-speed (UHS) digital holographic microscope (DHM) [2].

The goal is to learn reduced-order stochastic models that capture both:
the dominant deterministic dynamics, and the stochastic effects inherent in experimental measurements.

Experimental Data Source

The experimental measurements are obtained via off-axis digital holography microscopy, as described in:

[2] Y. Emery, T. Colomb, E. Cuche
 [Metrology applications using off-axis digital holography microscopy]
 (https://doi.org/10.1088/2515-7647/ac0957)
Journal of Physics: Photonics, 3(3):034016, 2021.
@article{emery2021metrology,
  title={Metrology applications using off-axis digital holography microscopy},
  author={Emery, Yves and Colomb, Tristan and Cuche, Etienne},
  journal={Journal of Physics: Photonics},
  volume={3},
  number={3},
  pages={034016},
  year={2021},
  publisher={IOP Publishing}
}

Requirements

MATLAB R2024b (or newer recommended)

Toolboxes as required by the OpInf implementation (see code comments)
