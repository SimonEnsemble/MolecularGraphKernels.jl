```@meta
DocTestSetup = quote
    using MolecularGraph, MolecularGraphKernels
end
```

# Molecules as Graphs

For convenience, `MolecularGraphKernels` can use either [`MetaGraph`](https://juliagraphs.org/MetaGraphs.jl/dev/) or [`MolecularGraph.GraphMol`](https://mojaie.github.io/MolecularGraph.jl/dev/) objects as inputs representing molecular structures.
A function to convert `GraphMol` objects to `MetaGraph` objects is exported; this conversion will be handled automatically if not done ahead of time, but should be done manually for applications where an input will be reused many times (for efficiency).

```jldoctest; output=false
graph_mol = smilestomol("c1ccccc1")
g = MetaGraph(graph_mol)
# output
{6, 6} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)
```
