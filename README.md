# FunBoot
A functional conformal bootstrap solver using outer approximation


To install, make sure you are using Julia 1.9.1.

using Pkg; Pkg.add(url="https://github.com/Canonical111/FunBoot")

For better performance, load the package AppleAccelerate on MacOS with Apple Silicon or mkl on Intel CPU.

"using AppleAccelerate" or "using mkl"
