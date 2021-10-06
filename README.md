| **Authors**  | **Project** | **Build Status**              | **Latest Version** | **License** |
|:------------:|:-----------:|:-----------------------------:|:------------------:|:-----------:|
|   N. Curti   |   **SysDyn**    | **Linux/MacOs** : [![Travis](https://travis-ci.com/Nico-Curti/SysDyn.svg?branch=master)](https://travis-ci.com/Nico-Curti/SysDyn) <br/> **Windows** : [![appveyor](https://ci.appveyor.com/api/projects/status/iwcl52c5ngx92w94?svg=true)](https://ci.appveyor.com/project/Nico-Curti/sysdyn)  | ![version](https://img.shields.io/badge/PyPI-v1.0.0-orange.svg?style=plastic) | [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://github.com/Nico-Curti/SysDyn/blob/master/LICENSE)

<a href="https://github.com/physycom">
  <div class="image">
    <img src="https://cdn.rawgit.com/physycom/templates/697b327d/logo_unibo.png" width="90" height="90">
  </div>
</a>

# Systems Dynamics Functions and Examples

The project is developed as tutorial for the courses *Complex System*, *Numerical Methods* and *Physical Methods of Biology* at the University of Bologna.

Several methods about **Systems Dynamics** research are shown from the simplest integration schemes, to a series of more advanced topics.
For the major part of codes you can find both the **Python** and **C++** versions: the simplest algorithms are written in **Python**, while the most computational expensive ones are written (only) in **C++**.

## Documentation

* [Metropolis-Hastings Algorithm](https://github.com/Nico-Curti/SysDyn/blob/master/docs/source/examples/Metropolis_algorithm.ipynb): Simple example of Monte Carlo algorithm for the computation of $\pi$.
  Starting from a naive implementation, we show a series of incremental refinements of the Python code, for the computational performances improvement.
  At the end of the example you can find also a (faster) C++ implementation of the algorithm based on multi-threading computation (ref. [here](https://github.com/Nico-Curti/SysDyn/blob/master/cpp/metropolis.cpp)).

* [Rate Equations](https://github.com/Nico-Curti/SysDyn/blob/master/docs/source/examples/lecture1_rate_equations.ipynb): Introduction about mathematical background of kinetic rate equations, starting from the zero-order kinetic to the more advanced conversion systems.
  For each example, the Python code is provided with a description of the fundamental steps.
  At the end of each sub-section, we show the corresponding C++ implementation (a putative implementation!) (ref. [zero-order](https://github.com/Nico-Curti/SysDyn/blob/master/cpp/zero_order_kinetic.cpp) & [first-order](https://github.com/Nico-Curti/SysDyn/blob/master/cpp/first_order_kinetic.cpp))

* [Michaelis-Menten](https://github.com/Nico-Curti/SysDyn/blob/master/docs/source/examples/lecture2_michaelis_menten.ipynb): Description of the Michaelis-Menten equations system.
  In this example we introduce the Runge-Kutta integration scheme, providing a possible implementation of it in the C++ version of the code (ref. [here](https://github.com/Nico-Curti/SysDyn/blob/master/cpp/michaelis_menten_rk4.cpp)).
  The RK method will be used in next examples, showing also a Python implementation of it.

* [Lotka Volterra](https://github.com/Nico-Curti/SysDyn/blob/master/docs/source/examples/lecture3_lotka_volterra.ipynb): Description and implementation of the Lotka Volterra model.
  In this example we introduce the basis of the Symbolic Computation for the stability analysis of the Lotka-Volterra system.
  For this example no C++ implementation is provided, since a simple edit of the Michaelis-Menten code can lead to the Lotka-Volterra integration.

* [Chemical Master Equation](https://github.com/Nico-Curti/SysDyn/blob/master/docs/source/examples/lecture5_CME.ipynb): Description of the mathematical background of Chemical Master Equation.
  In this example we introduce the Gillespie Algorithm for the integration of CME models, **TODO**

* [Brusselator](https://github.com/Nico-Curti/SysDyn/blob/master/docs/source/examples/advanced_application.ipynb): Description of the Brusselator kinetic system.
  In this example we combine together all the techniques that we have seen up to now, performing a complete analysis of the system for both the deterministic and stochastic version of the model.
  At the end of the example, the extension to 2D (with the integration of diffusion coefficients) is showed, leading to the creation of Turing patterns.

## Installation

### C++ version

To build the C++ scripts:

```bash
git clone https://github.com/Nico-Curti/SysDyn.git
cd SysDyn
mkdir build
cd build
cmake ..
cmake --build . --target install --config Release
```

**NOTE**: make sure to have a c++ compiler which supports the minimum standard required!
If some troubles occur, you can follow the instruction at [intrphysycom](https://github.com/physycom/intrphysycom) page to configure your machine; for issues related to softwares installation, you can use the scripts in [ShUt](https:://github.com//Nico-Curti/shut) if you are looking for *no root users* solutions.

### Python version

To use the python scripts install the prerequisites:

```bash
python -m pip install -r prerequisites.txt
```

## Contributions

Any contribution is more than welcome.
Just fill an issue or a pull request and I will check ASAP!

## Authors

* <img src="https://avatars0.githubusercontent.com/u/24650975?s=400&v=4" width="25px"> **Nico Curti** [git](https://github.com/Nico-Curti), [unibo](https://www.unibo.it/sitoweb/nico.curti2)

## License

The `SysDyn` package is licensed under the [GPLv3](https://github.com/Nico-Curti/SysDyn/blob/master/LICENSE) License.

## Citation

If you have found `SysDyn` helpful in your research, please consider citing the project repository

```BibTeX
@misc{SysDyn,
  author = {Curti, Nico},
  title = {SysDyn - System Dynamics functions and examples},
  year = {2021},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/Nico-Curti/SysDyn}},
}
```
