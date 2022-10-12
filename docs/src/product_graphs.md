```@meta
DocTestSetup = quote
    using InteractiveUtils, MolecularGraphKernels
    import MolecularGraphKernels.AbstractProductGraph
end
```

# Product Graphs

Product graphs are contained within the [`ProductGraph`](@ref) type.
The [`ProductGraph`](@ref) type is parameterized by subtypes of `AbstractProductGraph`.

```jldoctest
subtypes(AbstractProductGraph)

# output

2-element Vector{Any}:
 Direct
 Modular
```

While a [`ProductGraph`](@ref) contains the full representation of the product graph between two molecules, the adjacency matrix may be all that is required.
As such, it may be more efficient to work with the adjacency matrix (see [`product_graph_adjacency_matrix`](@ref)).

To create a new [`ProductGraph`](@ref), define a new `AbstractProductGraph` subtype and two functions, `record_adjacency!` and `product_graph_edge_label`, parameterized by the new type.

```@docs
ProductGraph
product_graph_adjacency_matrix
```
