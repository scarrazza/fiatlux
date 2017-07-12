# fiatlux
API for LUXqed methodology in global PDF fits.

## Project summary and aim

The aim of `libfiatlux` is to provide a blackbox tool which computed the photon PDF at a given Q value
using the LUX approach by Salam, Nason et al. The output of this repository is a C++ library
which can be imported and shared to other native C++ programs like fiatlux in nnpdfcpp. 

The library implements following features:
- Computes LUX photon by subdiving in elastic, inelastic and msbar components
- Allow variations of parameters to estimate uncertainties
- Generic interface to F2, FL and alpha QED: you can plug APFEL, HOPPET or any other code.

### Release and Tag policy

The library is tagged and released when a major and stable status is achieved. 
Tags and releases do not necessarily follows the NNPDF releases.

### Code development policy/rules

Currently the project is maintained by 1 developer.

### Code style

The code uses C++11.

### Continuous integration (CI)

We implement CI at CERN's gitlab. 

### Testing

Manual testes are available in the `examples` folder.

## Installation

`libfialux` depends on the following libraries:

- pkg-config
- gsl
- yaml-cpp

optinally to build the examples:
- lhapdf
- apfel
- hoppet QED

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
