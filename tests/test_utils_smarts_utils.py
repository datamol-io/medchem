from typing import cast

import datamol as dm

from medchem.utils.smarts import SMARTSUtils


def test_standard_attachment():
    smiles = [
        "C1CC(C)C(O)CC1[1*]",  # ortho
        "c1c(C)cc(O)cc1[*:2]",  # meta + aromatic
        "c1(C)ccc(O)cc1[*]",  # para + aromat
    ]

    for sm in smiles:
        new_sm = SMARTSUtils.standardize_attachment(sm, "[*:1]")
        assert "[*:1]" in new_sm


def test_ortho_meta_para():
    smiles = [
        "C1CC(C)C(O)CC1",  # ortho
        "c1c(C)cc(O)cc1",  # meta + aromatic
        "c1(C)ccc(O)cc1",  # para + aromat
    ]
    mols = [dm.to_mol(x) for x in smiles]
    mols = [cast(dm.Mol, x) for x in mols]

    sm1 = "[#6;!R]"
    sm2 = "[#8]"
    ortho_query = dm.from_smarts(SMARTSUtils.ortho(sm1, sm2))
    ortho_aro_query = dm.from_smarts(SMARTSUtils.ortho(sm1, sm2, aromatic_only=True))
    meta_query = dm.from_smarts(SMARTSUtils.meta(sm1, sm2))
    para_query = dm.from_smarts(SMARTSUtils.para(sm1, sm2))
    expected_output = [True, True, True, False, False]
    output = [
        mols[0].HasSubstructMatch(ortho_query),  # True
        mols[1].HasSubstructMatch(meta_query),  # True
        mols[2].HasSubstructMatch(para_query),  # True
        mols[0].HasSubstructMatch(para_query),  # False
        mols[0].HasSubstructMatch(ortho_aro_query),  # False
    ]

    assert output == expected_output


def test_long_chain():
    smiles = [
        "CCCCCCCC",  # octane
        "CCCC=CCCN",  #  hept-3-en-1-amine
        "CCOCCC",  #  1-ethoxypropane
        "CC(C)(C)CN",  #  2,2-dimethylpropan-1-amine
        "CCC(C)(C)CN",  # 2,2-dimethylbutan-1-amine
        "C1CCCC1CC",  # ethylcyclopentane
    ]
    mols = [dm.to_mol(x) for x in smiles]
    mols = [cast(dm.Mol, x) for x in mols]

    chain1 = dm.from_smarts(
        SMARTSUtils.aliphatic_chain(
            min_size=5,
            unsaturated_bondtype=None,
        ),
    )
    chain2 = dm.from_smarts(
        SMARTSUtils.aliphatic_chain(
            min_size=6,
            unbranched=True,
            unsaturated_bondtype=dm.SINGLE_BOND,
        )
    )
    chain3 = dm.from_smarts(
        SMARTSUtils.aliphatic_chain(
            min_size=6,
            unbranched=True,
            unsaturated_bondtype=dm.DOUBLE_BOND,
            allow_hetero_atoms=False,
        )
    )
    chain1_out = [m.HasSubstructMatch(chain1) for m in mols]
    chain2_out = [m.HasSubstructMatch(chain2) for m in mols]
    chain3_out = [m.HasSubstructMatch(chain3) for m in mols]

    assert chain1_out == [True, True, True, False, True, False]
    assert chain2_out == [True, False, True, False, False, False]
    assert chain3_out == [True, True, False, False, False, False]


def test_atom_in_env():
    mol = dm.to_mol("c1[c:1](OC)cc(F)cc1")
    mol = cast(dm.Mol, mol)

    # we are trying to match the carbon with atom map 1
    # these queries should be equivalent
    atom_id = tuple([atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() == 1])
    query1 = SMARTSUtils.atom_in_env("[#6;r6]", "[*][OD2][C&D1]", "[c]aa[F]", union=False)
    query2 = SMARTSUtils.atom_in_env("[#6;r6][OD2][C&D1]", "[c]aa[F]", union=False)

    res1 = mol.GetSubstructMatch(dm.from_smarts(query1))
    res2 = mol.GetSubstructMatch(dm.from_smarts(query2))

    assert res1 == res2
    assert res1 == atom_id
