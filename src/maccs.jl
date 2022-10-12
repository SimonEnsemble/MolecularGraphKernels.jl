"""
returns the MACCS fingerprint from a SMILES string

!!! warning
    
    Not supported in Windows (RDKitMinimalLib incompatibility).
"""
function maccs_fp(smiles::String)::BitVector
    mol = get_mol(smiles)
    match_counts = zeros(length(maccs_queries))
    match_counts[query_notnothing] .=
        [length(get_substruct_matches(mol, maccs_queries[i])) for i in query_notnothing_idx]
    return match_counts .> maccs_counts
end
