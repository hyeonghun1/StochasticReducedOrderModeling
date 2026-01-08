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
Compared to the original stochastic OpInf implementation in [1], this repository extends the framework to experimental data. In particular, it applies stochastic operator inference to capillary wave turbulence data acquired using a custom ultra–high-speed (UHS) digital holographic microscope (DHM), as described in:

[2] Y. Emery, T. Colomb, E. Cuche
 [Metrology applications using off-axis digital holography microscopy](https://doi.org/10.1088/2515-7647/ac0957)<details><summary>BibTex</summary><pre>
@article{emery2021metrology,
  title={Metrology applications using off-axis digital holography microscopy},
  author={Emery, Yves and Colomb, Tristan and Cuche, Etienne},
  journal={Journal of Physics: Photonics},
  volume={3},
  number={3},
  pages={034016},
  year={2021},
  publisher={IOP Publishing}
}</pre></details>

The goal is to learn reduced-order stochastic models that capture both:
the dominant deterministic dynamics, and the stochastic effects inherent in experimental measurements.

Experimental Data Source

The experimental measurements are obtained via DHM, as described in:
[3] J. Orosco, W. Connacher, J. Friend
[Identification of weakly to strongly-turbulent three-wave processes in a micro-scale system](https://www.sciencedirect.com/science/article/pii/S0960077923005167)<details><summary>BibTex</summary><pre>
@article{orosco2023identification,
  title={Identification of weakly to strongly-turbulent three-wave processes in a micro-scale system},
  author={Orosco, Jeremy and Connacher, William and Friend, James},
  journal={Chaos, Solitons \& Fractals},
  volume={172},
  pages={113615},
  year={2023},
  publisher={Elsevier}
}</pre></details>


Requirements

MATLAB R2024b (or newer recommended)
Toolboxes as required by the OpInf implementation (see code comments)
