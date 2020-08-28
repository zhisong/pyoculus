# pyoculus
A Python version of Oculus - The eye into the chaos: a comprehensive magnetic field diagnostic package for non-integrable, toroidal magnetic fields (and more general 1 1/2-D or 2D Hamiltonian system). Oculus is the Latin word for 'eye'.

## Usage

To use the package, simply import it in Python:
```python
import pyoculus
```

Examples can be found in the `examples` subfolder.

## Documentation

The documentation of pyoculus is managed by [Doxygen](https://www.doxygen.nl/index.html).
To generate the most up-to-date version of the documentation, please run

```
make doxygen
```

## SPEC magnetic field and Pressure Jump Hamiltonian (PJH) 

### Pre-requisite
The `py_spec` package that handles SPEC output.

### Compilation
To use the pyoculus on SPEC magnetic field and PJH, please compile the Fortran modules for SPEC by
```
make
```
One may need to change `F90` and `F90FLAGS` in `Makefile` to your local Fortran compiler.

The Makefile uses `f2py` from the `numpy` package to create an interface with Python.
Documentation for `f2py` can be found [here](https://numpy.org/doc/stable/f2py/).

## Link to the original Oculus package:

Github: https://github.com/SRHudson/Oculus

Documentation: https://w3.pppl.gov/~shudson/Oculus/oculus.pdf

