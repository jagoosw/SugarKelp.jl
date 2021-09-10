# Kelp.jl

Implementation of the [Broch and Slagstad, 2012 model of the growth and composition of _Saccharina latissima_](https://link.springer.com/article/10.1007/s10811-011-9695-y).

The main way to solve a single frond is `Kelp.solvekelp` and grids can be solved by `Kelp.solvegrid`.

Changes from the stated parameter values in the paper are detailed in [changes.pdf](https://github.com/jagoosw/Kelp/blob/main/changes.pdf).

The package is not yet registered so to use, download this repository and then install the dependencies by executing (from this directory):
```
>julia
julia> import Pkg
julia> ] activate .
julia> instantiate
```

## Examples
In the examples folder are is currently one example that attempts to reproduce the origional papers results. Here is the output:

![Figure 3 equivalent.](img/paper_comparison.png)

I will replace this with a better example and examples of solving a grid at some point.

## Documentation
The code currently has (hopefully) accurate docstrings which should explain how everything works. If you would like to view this locally like a readthedocs page execute:
```
>cd doc
>julia make.jl
```
You will then find documentation at doc/build/index.html
Please don't push this to the repository, it will eventually go in another branch.