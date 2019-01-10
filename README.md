| **Authors**  | **Project** | **Build Status**              | **Latest Version** | **License** |
|:------------:|:-----------:|:-----------------------------:|:------------------:|:-----------:|
|   N. Curti   |   SysDyn    | **Linux/MacOs** : [![Travis](https://travis-ci.com/Nico-Curti/SysDyn.svg?branch=master)](https://travis-ci.com/Nico-Curti/SysDyn) <br/> **Windows** : [![appveyor](https://ci.appveyor.com/api/projects/status/iwcl52c5ngx92w94?svg=true)](https://ci.appveyor.com/project/Nico-Curti/sysdyn)  | ![version](https://img.shields.io/badge/PyPI-v1.0.0-orange.svg?style=plastic) | [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://github.com/Nico-Curti/SysDyn/blob/master/LICENSE)

<a href="https://github.com/physycom">
<div class="image">
<img src="https://cdn.rawgit.com/physycom/templates/697b327d/logo_unibo.png" width="90" height="90">
</div>
</a>

# Systems Dynamics Functions and Examples

The project is developed as tutorial for the courses *Complex System* and *Numerical Methods* at the University of Bologna.

Different methods about **Systems Dynamics** research are shown from the simplest ones to more advanced topics. The simplest code are written in **Python** language while the most (computationally) expensive ones are written in **C++**.

## Installation

### C++ version

To build the C++ scripts:

```bash
git clone https://github.com/Nico-Curti/SysDyn.git
cd SysDyn
mkdir build
cd build
cmake ..
cmake --build .
```

**NOTE**: make sure that your c++ compiler supports the standard required! If some troubles occur follow the instruction at [intrphysycom](https://github.com/physycom/intrphysycom) page to configure your machine or use the scripts in [ShUt](https:://github.com//Nico-Curti/shut) if you are looking for *no root users* solutions.

### Python version

To use the python scripts install the prerequisites:

```bash
pip install -r prerequisites.txt
```

## Contributions

Any contribution is more than welcome. Just fill an issue or a pull request and I will check ASAP!

## Authors

* **Nico Curti** [git](https://github.com/Nico-Curti), [unibo](https://www.unibo.it/sitoweb/nico.curti2)

## License

This project is released under GPL license.
