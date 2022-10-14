"""
From RDKit
https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/MACCSkeys.py
"""
maccs_queries = [
    (raw"?", 0),
    #(raw"[#104,#105,#106,#107,#106,#109,#110,#111,#112]",0),  # atomic num >103 Not complete
    (raw"[#104]", 0),  # limit the above def"n since the RDKit only accepts up to #104
    (raw"[#32,#33,#34,#50,#51,#52,#82,#83,#84]", 0),  # Group IVa,Va,VIa Rows 4-6
    (raw"[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]", 0),  # actinide
    (raw"[Sc,Ti,Y,Zr,Hf]", 0),  # Group IIIB,IVB (Sc...)
    (raw"[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]", 0),  # Lanthanide
    (raw"[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]", 0),  # Group VB,VIB,VIIB
    (raw"[!#6;!#1]1~*~*~*~1", 0),  # QAAA@1
    (raw"[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]", 0),  # Group VIII (Fe...)
    (raw"[Be,Mg,Ca,Sr,Ba,Ra]", 0),  # Group IIa (Alkaline earth)
    (raw"*1~*~*~*~1", 0),  # 4M Ring
    (raw"[Cu,Zn,Ag,Cd,Au,Hg]", 0),  # Group IB,IIB (Cu..)
    (raw"[#8]~[#7](~[#6])~[#6]", 0),  # ON(C)C
    (raw"[#16]-[#16]", 0),  # S-S
    (raw"[#8]~[#6](~[#8])~[#8]", 0),  # OC(O)O
    (raw"[!#6;!#1]1~*~*~1", 0),  # QAA@1
    (raw"[#6]#[#6]", 0),  #CTC
    (raw"[#5,#13,#31,#49,#81]", 0),  # Group IIIA (B...)
    (raw"*1~*~*~*~*~*~*~1", 0),  # 7M Ring
    (raw"[#14]", 0),  #Si
    (raw"[#6]=[#6](~[!#6;!#1])~[!#6;!#1]", 0),  # C=C(Q)Q
    (raw"*1~*~*~1", 0),  # 3M Ring
    (raw"[#7]~[#6](~[#8])~[#8]", 0),  # NC(O)O
    (raw"[#7]-[#8]", 0),  # N-O
    (raw"[#7]~[#6](~[#7])~[#7]", 0),  # NC(N)N
    (raw"[#6]=;@[#6](@*)@*", 0),  # C$=C($A)$A
    (raw"[I]", 0),  # I
    (raw"[!#6;!#1]~[CH2]~[!#6;!#1]", 0),  # QCH2Q
    (raw"[#15]", 0),  # P
    (raw"[#6]~[!#6;!#1](~[#6])(~[#6])~*", 0),  # CQ(C)(C)A
    (raw"[!#6;!#1]~[F,Cl,Br,I]", 0),  # QX
    (raw"[#6]~[#16]~[#7]", 0),  # CSN
    (raw"[#7]~[#16]", 0),  # NS
    (raw"[CH2]=*", 0),  # CH2=A
    (raw"[Li,Na,K,Rb,Cs,Fr]", 0),  # Group IA (Alkali Metal)
    (raw"[#16R]", 0),  # S Heterocycle
    (raw"[#7]~[#6](~[#8])~[#7]", 0),  # NC(O)N
    (raw"[#7]~[#6](~[#6])~[#7]", 0),  # NC(C)N
    (raw"[#8]~[#16](~[#8])~[#8]", 0),  # OS(O)O
    (raw"[#16]-[#8]", 0),  # S-O
    (raw"[#6]#[#7]", 0),  # CTN
    (raw"F", 0),  # F
    (raw"[!#6;!#1;!H0]~*~[!#6;!#1;!H0]", 0),  # QHAQH
    (raw"[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]", 0),  # OTHER
    (raw"[#6]=[#6]~[#7]", 0),  # C=CN
    (raw"Br", 0),  # BR
    (raw"[#16]~*~[#7]", 0),  # SAN
    (raw"[#8]~[!#6;!#1](~[#8])(~[#8])", 0),  # OQ(O)O
    (raw"[!+0]", 0),  # CHARGE
    (raw"[#6]=[#6](~[#6])~[#6]", 0),  # C=C(C)C
    (raw"[#6]~[#16]~[#8]", 0),  # CSO
    (raw"[#7]~[#7]", 0),  # NN
    (raw"[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]", 0),  # QHAAAQH
    (raw"[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]", 0),  # QHAAQH
    (raw"[#8]~[#16]~[#8]", 0),  #OSO
    (raw"[#8]~[#7](~[#8])~[#6]", 0),  # ON(O)C
    (raw"[#8R]", 0),  # O Heterocycle
    (raw"[!#6;!#1]~[#16]~[!#6;!#1]", 0),  # QSQ
    (raw"[#16]!:*:*", 0),  # Snot%A%A
    (raw"[#16]=[#8]", 0),  # S=O
    (raw"*~[#16](~*)~*", 0),  # AS(A)A
    (raw"*@*!@*@*", 0),  # A$!A$A
    (raw"[#7]=[#8]", 0),  # N=O
    (raw"*@*!@[#16]", 0),  # A$A!S
    (raw"c:n", 0),  # C%N
    (raw"[#6]~[#6](~[#6])(~[#6])~*", 0),  # CC(C)(C)A
    (raw"[!#6;!#1]~[#16]", 0),  # QS
    (raw"[!#6;!#1;!H0]~[!#6;!#1;!H0]", 0),  # QHQH (&...) SPEC Incomplete
    (raw"[!#6;!#1]~[!#6;!#1;!H0]", 0),  # QQH
    (raw"[!#6;!#1]~[#7]~[!#6;!#1]", 0),  # QNQ
    (raw"[#7]~[#8]", 0),  # NO
    (raw"[#8]~*~*~[#8]", 0),  # OAAO
    (raw"[#16]=*", 0),  # S=A
    (raw"[CH3]~*~[CH3]", 0),  # CH3ACH3
    (raw"*!@[#7]@*", 0),  # A!N$A
    (raw"[#6]=[#6](~*)~*", 0),  # C=C(A)A
    (raw"[#7]~*~[#7]", 0),  # NAN
    (raw"[#6]=[#7]", 0),  # C=N
    (raw"[#7]~*~*~[#7]", 0),  # NAAN
    (raw"[#7]~*~*~*~[#7]", 0),  # NAAAN
    (raw"[#16]~*(~*)~*", 0),  # SA(A)A
    (raw"*~[CH2]~[!#6;!#1;!H0]", 0),  # ACH2QH
    (raw"[!#6;!#1]1~*~*~*~*~1", 0),  # QAAAA@1
    (raw"[NH2]", 0),  #NH2
    (raw"[#6]~[#7](~[#6])~[#6]", 0),  # CN(C)C
    (raw"[C;H2,H3][!#6;!#1][C;H2,H3]", 0),  # CH2QCH2
    (raw"[F,Cl,Br,I]!@*@*", 0),  # X!A$A
    (raw"[#16]", 0),  # S
    (raw"[#8]~*~*~*~[#8]", 0),  # OAAAO
    (
        raw"[$([!#6;!#1;!H0]~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]",
        0
    ),  # QHAACH2A
    (
        raw"[$([!#6;!#1;!H0]~*~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]",
        0
    ),  # QHAAACH2A
    (raw"[#8]~[#6](~[#7])~[#6]", 0),  # OC(N)C
    (raw"[!#6;!#1]~[CH3]", 0),  # QCH3
    (raw"[!#6;!#1]~[#7]", 0),  # QN
    (raw"[#7]~*~*~[#8]", 0),  # NAAO
    (raw"*1~*~*~*~*~1", 0),  # 5 M ring
    (raw"[#7]~*~*~*~[#8]", 0),  # NAAAO
    (raw"[!#6;!#1]1~*~*~*~*~*~1", 0),  # QAAAAA@1
    (raw"[#6]=[#6]", 0),  # C=C
    (raw"*~[CH2]~[#7]", 0),  # ACH2N
    (
        raw"[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]",
        0
    ),  # 8M Ring or larger. This only handles up to ring sizes of 14
    (raw"[!#6;!#1]~[#8]", 0),  # QO
    (raw"Cl", 0),  # CL
    (raw"[!#6;!#1;!H0]~*~[CH2]~*", 0),  # QHACH2A
    (raw"*@*(@*)@*", 0),  # A$A($A)$A
    (raw"[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]", 0),  # QA(Q)Q
    (raw"[F,Cl,Br,I]~*(~*)~*", 0),  # XA(A)A
    (raw"[CH3]~*~*~*~[CH2]~*", 0),  # CH3AAACH2A
    (raw"*~[CH2]~[#8]", 0),  # ACH2O
    (raw"[#7]~[#6]~[#8]", 0),  # NCO
    (raw"[#7]~*~[CH2]~*", 0),  # NACH2A
    (raw"*~*(~*)(~*)~*", 0),  # AA(A)(A)A
    (raw"[#8]!:*:*", 0),  # Onot%A%A
    (raw"[CH3]~[CH2]~*", 0),  # CH3CH2A
    (raw"[CH3]~*~[CH2]~*", 0),  # CH3ACH2A
    (raw"[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]", 0),  # CH3AACH2A
    (raw"[#7]~*~[#8]", 0),  # NAO
    (raw"[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]", 1),  # ACH2CH2A > 1
    (raw"[#7]=*", 0),  # N=A
    (raw"[!#6;R]", 1),  # Heterocyclic atom > 1 (&...) Spec Incomplete
    (raw"[#7;R]", 0),  # N Heterocycle
    (raw"*~[#7](~*)~*", 0),  # AN(A)A
    (raw"[#8]~[#6]~[#8]", 0),  # OCO
    (raw"[!#6;!#1]~[!#6;!#1]", 0),  # QQ
    (raw"?", 0),  # Aromatic Ring > 1
    (raw"*!@[#8]!@*", 0),  # A!O!A
    (raw"*@*!@[#8]", 1),  # A$A!O > 1 (&...) Spec Incomplete
    (
        raw"[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$(*~[CH2]~*~[R]1@[R]@[CH2;R]1)]",
        0
    ),  # ACH2AAACH2A
    (
        raw"[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[CH2;R]1)]",
        0
    ),  # ACH2AACH2A
    (raw"[!#6;!#1]~[!#6;!#1]", 1),  # QQ > 1 (&...)  Spec Incomplete
    (raw"[!#6;!#1;!H0]", 1),  # QH > 1
    (raw"[#8]~*~[CH2]~*", 0),  # OACH2A
    (raw"*@*!@[#7]", 0),  # A$A!N
    (raw"[F,Cl,Br,I]", 0),  # X (HALOGEN)
    (raw"[#7]!:*:*", 0),  # Nnot%A%A
    (raw"[#8]=*", 1),  # O=A>1
    (raw"[!C;!c;R]", 0),  # Heterocycle
    (raw"[!#6;!#1]~[CH2]~*", 1),  # QCH2A>1 (&...) Spec Incomplete
    (raw"[O;!H0]", 0),  # OH
    (raw"[#8]", 3),  # O > 3 (&...) Spec Incomplete
    (raw"[CH3]", 2),  # CH3 > 2  (&...) Spec Incomplete
    (raw"[#7]", 1),  # N > 1
    (raw"*@*!@[#8]", 0),  # A$A!O
    (raw"*!:*:*!:*", 0),  # Anot%A%Anot%A
    (raw"*1~*~*~*~*~*~1", 1),  # 6M ring > 1
    (raw"[#8]", 2),  # O > 2
    (raw"[$(*~[CH2]~[CH2]~*),$([R]1@[CH2;R]@[CH2;R]1)]", 0),  # ACH2CH2A
    (raw"*~[!#6;!#1](~*)~*", 0),  # AQ(A)A
    (raw"[C;H3,H4]", 1),  # CH3 > 1
    (raw"*!@*@*!@*", 0),  # A!A$A!A
    (raw"[#7;!H0]", 0),  # NH
    (raw"[#8]~[#6](~[#6])~[#6]", 0),  # OC(C)C
    (raw"[!#6;!#1]~[CH2]~*", 0),  # QCH2A
    (raw"[#6]=[#8]", 0),  # C=O
    (raw"*!@[CH2]!@*", 0),  # A!CH2!A
    (raw"[#7]~*(~*)~*", 0),  # NA(A)A
    (raw"[#6]-[#8]", 0),  # C-O
    (raw"[#6]-[#7]", 0),  # C-N
    (raw"[#8]", 1),  # O>1
    (raw"[C;H3,H4]", 0),  #CH3
    (raw"[#7]", 0),  # N
    (raw"a", 0),  # Aromatic
    (raw"*1~*~*~*~*~*~1", 0),  # 6M Ring
    (raw"[#8]", 0),  # O
    (raw"[R]", 0),  # Ring
    (raw"?", 0)  # Fragments  FIX: this can"t be done in SMARTS
]

# pre-compute the queries. must ignore "?". replace with nothing.
if ! Sys.iswindows()
    _maccs_queries = [(smarts_pattern != "?" ? get_qmol(smarts_pattern) : nothing, nb_matches) 
                          for (smarts_pattern, nb_matches) in maccs_queries]
end

"""
returns the MACCS fingerprint from a SMILES string

!!! warning
    
    Not supported in Windows (RDKitMinimalLib incompatibility).
"""
function maccs_fp(smiles::String)::BitVector
    mol = get_mol(smiles)
    x = falses(length(maccs_queries))
    for (i, (pattern, nb_matches_needed)) in enumerate(_maccs_queries)
        if ! isnothing(pattern)
            nb_matches = length(get_substruct_matches(mol, pattern))
            if nb_matches > nb_matches_needed
                x[i] = true
            end
        end
    end
    return x
end
