using BenchmarkTools, MolecularGraph, MolecularGraphKernels

g₁, g₂ = MetaGraph.(
    smilestomol.(
        [
            "Oc1c(c(O)cc(c1)CCCCC)[C@@H]2\\C=C(/CC[C@H]2\\C(=C)C)C" # cannabidiol
            "Oc1cc(cc(O)c1C/C=C(\\C)CC\\C=C(/C)C)CCCCC" # cannabigerol
        ]
    )
)

@btime ccsi(g₁, g₂)
