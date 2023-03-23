module Test_visualization

using GraphMakie, Graphs, MetaGraphs, MolecularGraph, MolecularGraphKernels, Test
import MolecularGraphKernels: VizGraphKwargs

function test_vis(graph, graph_name, set_name; kwargs...)
    @testset "$set_name" begin
        rm("$graph_name.pdf"; force=true)
        @test !isnothing(viz_graph(graph; savename=graph_name, kwargs...))
        vis = isfile("$graph_name.pdf")
        @test vis
        if vis
            @info "$set_name test visualization in $graph_name.pdf"
        end
    end
end

@testset verbose = true "VizGraphKwargs" begin
    @test_throws Exception VizGraphKwargs()
    @test !isnothing(VizGraphKwargs(smilestomol("C")))
    @test !isnothing(VizGraphKwargs(smilestomol("C"); savename="foo"))
    @test !isnothing(VizGraphKwargs(MetaGraph(1)))
end

@testset verbose = true "Graph Types" begin
    test_vis(smilestomol("c1ccccc1"), "benzene", "MetaGraph"; layout=MolecularGraphKernels.Molecular)
    test_vis(
        ProductGraph{Direct}(smilestomol("NC=O"), smilestomol("C(NC=O)NC=O")),
        "dpg",
        "Direct Product Graph"
    )
    test_vis(
        ProductGraph{Modular}(smilestomol("NC=O"), smilestomol("C(NC=O)NC=O")),
        "mpg",
        "Modular Product Graph"
    )
end

@testset verbose = true "Graph Plot Styles" begin
    mol = smilestomol("C(NC=O)NC=O")
    styles = Dict(
        "spring" => GraphMakie.Spring,  
        "spectral" => GraphMakie.Spectral, 
        "stress" => GraphMakie.Stress,
        "molecular" => MolecularGraphKernels.Molecular
    )
    for (k, v) in styles
        test_vis(mol, k, k; layout=v)
    end
end

@testset verbose = true "Alpha Mask" begin
    mpg = ProductGraph{Modular}(smilestomol("C1C(C=O)C1"), smilestomol("C1C(=O)C1"))
    @testset "edges" begin
        edge_mask = [get_prop(mpg, e, :label) for e in edges(mpg)]
        test_vis(mpg, "edge_mask", "edge mask"; edge_alpha_mask=edge_mask)
    end
    @testset "nodes" begin
        node_mask = ones(nv(mpg))
        node_mask[[1, 3, 5, 7]] .= [0]
        test_vis(mpg, "node_mask", "node mask"; node_alpha_mask=node_mask)
    end
end

@testset verbose = true "Examples" begin
    g₁ = MetaGraph(smilestomol("COP(=O)(OC)OC(Br)C(Cl)(Cl)Br"))
    g₂ = MetaGraph(smilestomol("COP(N)(=O)SC"))
    dpg = ProductGraph{Direct}(g₁, g₂)
    test_vis(dpg, "bee-tox", "example")
end

end
