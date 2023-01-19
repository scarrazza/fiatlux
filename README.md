# fiatlux
API for LUXqed methodology in global PDF fits.

[![DOI](https://zenodo.org/badge/113965474.svg)](https://zenodo.org/badge/latestdoi/113965474)

## Project summary and aim

The aim of `libfiatlux` is to provide a blackbox tool which computed the photon PDF at a given Q value
using the LUX approach by Manohar, Nason, Salam and Zanderighi in [arXiv:1607.04266](https://arxiv.org/abs/1607.04266) and [arXiv:1708.01256](https://arxiv.org/abs/1708.01256). The output of this repository is a C++ library
which can be imported and shared to other programs.

The library implements following features:
- Computes LUX photon by subdividing in elastic, inelastic and msbar components
- Allow variations of parameters to estimate uncertainties
- Generic interface to F2, FL and alpha QED: you can plug APFEL or any other evolution code.

### Release and Tag policy

The library is tagged and released when a major and stable status is achieved.

### Testing

Manual testes are available in the `examples` folder.

## Installation

### Python library

```Shell
pip install .
```

### C++ library

`libfialux` depends on the following libraries:

- pkg-config
- yaml-cpp

optinally to build the examples:
- lhapdf
- apfel

please ensure to have the dependencies correctly installed and in your PATH before installing libfiatlux.

#### Configurations

Possible configurations:

```Shell
cmake .

```
or (recommended):

```Shell
mkdir build
cd build
cmake ..

```
You can control the optional flags with `ccmake` or from cmd line, the most relevant flags are:

```Shell
CMAKE_INSTALL_PREFIX
ENABLE_EXAMPLES
```

On the command line, options are controlled appending a `-D` flag. For
example:

```
cmake .. -DENABLE_EXAMPLES=on
```

## Documentation

### Code documentation

The code is documented with Doxygen (folder doc/), if you find methods or classes not fully documented open a issue request.

## Citation policy

If you decide to use this code please cite the following papers:

- The NNPDF3.1QED paper which is the fundamental motivation for this library [arXiv:1712.07053](https://arxiv.org/abs/1712.07053)
- The original LUX paper [arXiv:1607.04266](https://arxiv.org/abs/1607.04266)
- The long/complete version of LUX [arXiv:1708.01256](https://arxiv.org/abs/1708.01256)
- The GD11-P fit code from: The HERMES Collaboration [A. Airapetian et al.], JHEP 05 (2011) 126.
- The CLAS parametrization used in hep-ph/0301204 (CLAS) and described in hep-ph/9901360.
