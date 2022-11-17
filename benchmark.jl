using BenchmarkTools, IOCapture, MolecularGraph
import Pkg
IOCapture.capture() do
    Pkg.activate(".")
end
using MolecularGraphKernels

g₁, g₂ = MetaGraph.(
    smilestomol.(
        [
            "Oc1c(c(O)cc(c1)CCCCC)[C@@H]2\\C=C(/CC[C@H]2\\C(=C)C)C" # cannabidiol
            "Oc1cc(cc(O)c1C/C=C(\\C)CC\\C=C(/C)C)CCCCC" # cannabigerol
        ]
    )
)

capture = IOCapture.capture() do 
    return eval(quote @btime ccsi(g₁, g₂) end)
end

k = capture.value
t = capture.output

@info "CCSI on Two Graphs: $t"
if k ≠ 6205
    @error "k = $k"
end

capture = IOCapture.capture() do
    return eval(quote @btime ProductGraph{Modular}(g₁, g₂) end)
end

mpg = capture.value
t = capture.output

@info "Product Graph: $t"

capture = IOCapture.capture() do 
    return eval(quote @btime ccsi(mpg) end)
end

k = capture.value
t = capture.output

@info "Kernel: $t"
if k ≠ 6205
    @error "k = $k"
end
