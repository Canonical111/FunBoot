# FunBoot

FunBoot is a functional conformal bootstrap solver leveraging the power of outer approximation. This project aims to help with computations related to conformal bootstrap, a technique that utilizes consistency conditions of Conformal Field Theories (CFTs).

## Installation

Ensure you're using Julia version 1.9.1 or above. You can add the FunBoot package via:

```julia
using Pkg; Pkg.add(url="https://github.com/Canonical111/FunBoot")
```

## Performance Boost

For enhanced performance, it is recommended to load the package `AppleAccelerate` on MacOS with Apple Silicon or `mkl` on Intel CPUs.

To use these, run:

```julia
using AppleAccelerate
```

or

```julia
using MKL
```

## Testing

An example test, reproducing the 2D Ising model, is provided within `test/runtest.jl`. It loads the h5 file `test/2dIsing.h5`.

To execute the test file, run the following:

```julia
using Pkg; Pkg.test("FunBoot")
```

## Upcoming Features

We aim to supplement FunBoot with a Mathematica package that will generate the h5 functional table. This will be added in future updates, enhancing the overall functionality of the FunBoot project.
```
