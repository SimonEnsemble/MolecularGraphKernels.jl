module Test_maccs

using MolecularGraphKernels, Test
import MolecularGraphKernels: maccs_fp

@testset verbose = true "MACCS Fingerprint" begin
    if !Sys.iswindows()
        smiles = "CCOP(=S)(OCC)Oc1cc(C)nc(C(C)C)n1"

        # result from RDKit
        rdkit_on_bits = [
            29,
            38,
            48,
            65,
            67,
            73,
            74,
            77,
            80,
            81,
            88,
            92,
            97,
            98,
            102,
            106,
            109,
            110,
            112,
            113,
            114,
            117,
            120,
            121,
            124,
            126,
            127,
            128,
            130,
            137,
            138,
            141,
            142,
            143,
            146,
            148,
            149,
            153,
            155,
            156,
            157,
            159,
            160,
            161,
            162,
            163,
            164,
            165
        ]

        fingerprint = maccs_fp(smiles)

        @test rdkit_on_bits == findall(fingerprint)
    else
        @warn "MACCS fingerprinting not supported on Windows! Skipping test."
    end
end

end
