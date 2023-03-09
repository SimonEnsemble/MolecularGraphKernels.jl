var documenterSearchIndex = {"docs":
[{"location":"graph_kernels/#Graph-Kernels","page":"Graph Kernels","title":"Graph Kernels","text":"","category":"section"},{"location":"molecular_graphs/","page":"Molecular Graphs","title":"Molecular Graphs","text":"DocTestSetup = quote\n    using MolecularGraph, MolecularGraphKernels\nend","category":"page"},{"location":"molecular_graphs/#Molecules-as-Graphs","page":"Molecular Graphs","title":"Molecules as Graphs","text":"","category":"section"},{"location":"molecular_graphs/","page":"Molecular Graphs","title":"Molecular Graphs","text":"For convenience, MolecularGraphKernels can use either MetaGraph or MolecularGraph.GraphMol objects as inputs representing molecular structures. A function to convert GraphMol objects to MetaGraph objects is exported; this conversion will be handled automatically if not done ahead of time, but should be done manually for applications where an input will be reused many times (for efficiency).","category":"page"},{"location":"molecular_graphs/","page":"Molecular Graphs","title":"Molecular Graphs","text":"graph_mol = smilestomol(\"c1ccccc1\")\ng = MetaGraph(graph_mol)","category":"page"},{"location":"product_graphs/","page":"Product Graphs","title":"Product Graphs","text":"DocTestSetup = quote\n    using InteractiveUtils, MolecularGraphKernels\n    import MolecularGraphKernels.AbstractProductGraph\nend","category":"page"},{"location":"product_graphs/#Product-Graphs","page":"Product Graphs","title":"Product Graphs","text":"","category":"section"},{"location":"product_graphs/","page":"Product Graphs","title":"Product Graphs","text":"Product graphs are contained within the ProductGraph type. The ProductGraph type is parameterized by subtypes of AbstractProductGraph.","category":"page"},{"location":"product_graphs/","page":"Product Graphs","title":"Product Graphs","text":"subtypes(AbstractProductGraph)\n\n# output\n\n3-element Vector{Any}:\n Direct\n Modular\n Weighted","category":"page"},{"location":"product_graphs/","page":"Product Graphs","title":"Product Graphs","text":"While a ProductGraph contains the full representation of the product graph between two molecules, the adjacency matrix may be all that is required. As such, it may be more efficient to work with the adjacency matrix (see GraphMatrix).","category":"page"},{"location":"product_graphs/","page":"Product Graphs","title":"Product Graphs","text":"To create a new ProductGraph, define a new AbstractProductGraph subtype and two functions, record_adjacency! and product_graph_edge_label, parameterized by the new type.","category":"page"},{"location":"product_graphs/","page":"Product Graphs","title":"Product Graphs","text":"ProductGraph\nGraphMatrix","category":"page"},{"location":"product_graphs/#MolecularGraphKernels.ProductGraph","page":"Product Graphs","title":"MolecularGraphKernels.ProductGraph","text":"type-parameterized struct for product graphs\n\n\n\n\n\n","category":"type"},{"location":"visualization/#Vizualization","page":"Visualization","title":"Vizualization","text":"","category":"section"},{"location":"visualization/","page":"Visualization","title":"Visualization","text":"viz_graph","category":"page"},{"location":"visualization/#MolecularGraphKernels.viz_graph","page":"Visualization","title":"MolecularGraphKernels.viz_graph","text":"Visualize a molecular or product graph\n\n\n\n\n\n","category":"function"},{"location":"#MolecularGraphKernels.jl","page":"MolecularGraphKernels","title":"MolecularGraphKernels.jl","text":"","category":"section"},{"location":"","page":"MolecularGraphKernels","title":"MolecularGraphKernels","text":"Graph kernels compare two graphs and return a measure of their similarity:","category":"page"},{"location":"","page":"MolecularGraphKernels","title":"MolecularGraphKernels","text":"phi(G_1 G_2) = k","category":"page"},{"location":"","page":"MolecularGraphKernels","title":"MolecularGraphKernels","text":"This is useful, for example, when one wishes to train a model on graph-structured inputs of varying size, topology, and node/edge metadata (as one would do in molecule-based machine learning). In practice, this frequently entails computing the product graph of G_1 and G_2 and passing this new graph to a function which operates on its adjacency matrix:","category":"page"},{"location":"","page":"MolecularGraphKernels","title":"MolecularGraphKernels","text":"G_1x2 = langle G_1 G_2rangle","category":"page"},{"location":"","page":"MolecularGraphKernels","title":"MolecularGraphKernels","text":"psi(Adj(G_1x2)) = k","category":"page"},{"location":"","page":"MolecularGraphKernels","title":"MolecularGraphKernels","text":"MolecularGraphKernels.jl provides an interface for calculating product graphs, adjacency matrices thereof, and kernel functions between graph representations of molecules.","category":"page"},{"location":"#Installation","page":"MolecularGraphKernels","title":"Installation","text":"","category":"section"},{"location":"","page":"MolecularGraphKernels","title":"MolecularGraphKernels","text":"To install MolecularGraphKernels.jl simply add it via the Julia package manager:","category":"page"},{"location":"","page":"MolecularGraphKernels","title":"MolecularGraphKernels","text":"import Pkg\nPkg.add(\"MolecularGraphKernels\")","category":"page"},{"location":"#Collaboration","page":"MolecularGraphKernels","title":"Collaboration","text":"","category":"section"},{"location":"","page":"MolecularGraphKernels","title":"MolecularGraphKernels","text":"MolecularGraphKernels.jl is in active development. To report bugs, request features, share ideas, or ask questions, please open an issue on GitHub. Pull requests are welcome!","category":"page"}]
}
