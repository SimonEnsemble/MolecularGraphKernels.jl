# MolecularGraphKernels.jl

Graph kernels compare two graphs and return a measure of their similarity:

$\phi(G_1, G_2) = k$

This is useful, for example, when one wishes to train a model on graph-structured inputs of varying size, topology, and node/edge metadata (as one would do in molecule-based machine learning).
In practice, this frequently entails computing the product graph of $G_1$ and $G_2$ and passing this new graph to a function which operates on its adjacency matrix:

$G_{1x2} = \langle G_1, G_2\rangle$

$\psi(Adj(G_{1x2})) = k$

`MolecularGraphKernels.jl` provides an interface for calculating product graphs, adjacency matrices thereof, and kernel functions between graph representations of molecules.

## Installation

To install `MolecularGraphKernels.jl` simply add it via the Julia package manager:

```julia
import Pkg
Pkg.add("MolecularGraphKernels")
```

## Collaboration

`MolecularGraphKernels.jl` is in active development.
To report bugs, request features, share ideas, or ask questions, please [open an issue on GitHub](https://github.com/SimonEnsemble/MolecularGraphKernels.jl/issues).
Pull requests are welcome!
