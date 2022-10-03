```@meta
DocTestSetup = quote
    using InteractiveUtils, MolecularGraphKernels
    import MolecularGraphKernels.AbstractProductGraph
end
```

# Product Graphs

Product graphs are contained within the [`ProductGraph`](@ref) type.
The [`ProductGraph`](@ref) type is parameterized by subtypes of [`AbstractProductGraph`](@ref).

```jldoctest
subtypes(AbstractProductGraph)
# output
2-element Vector{Any}:
 Direct
 Factor
```

The adjacency matrix of a product graph is stored in the [`ProductGraphMatrix`](@ref) type.
[`ProductGraphMatrix`](@ref) is also parameterized by subtypes of [`AbstractProductGraph`](@ref).

While a [`ProductGraph`](@ref) contains the full representation of the product graph between two molecules, the adjacency matrix may be all that is required.
As such, it may be more efficient to work with the [`ProductGraphMatrix`](@ref).

To create a new [`ProductGraph`](@ref) or [`ProductGraphMatrix`](@ref), define a new [`AbstractProductGraph`](@ref) subtype and two functions, `mark_A!` and `product_graph_edge_label`, parameterized by the new type.

```@docs
```
