from medchem.utils.graph import score_symmetry


def test_score_symmetry():
    test_mols = [
        # high symmetry ~0.857
        ("O=C(O)c1cc(-n2ccnc2)cc(-n2ccnc2)c1", False),
        # fully symmetric
        ("OC(O)c1cc(-n2ccnc2)cc(-n2ccnc2)c1", True),
        # symmetry should be 1, modif of above
        ("CC(C)(C)[C@@H]1COC(C2(C3=N[C@H](C(C)(C)C)CO3)Cc3ccccc3C2)=N1", True),
        ("c1ccc2oc(-c3ccc(-c4nc5ccccc5o4)s3)nc2c1", True),
        # symmetrical star graph
        ("CC(C)(C)", True),
        ("NC(C)(N)C", True),
        # non symmetrical
        ("NC(C)(C)C", False),
        # modification of the above. not symmetrical
        ("CCCC1=CNC(CC)=C1C", False),
        # disconnection does not prevent symmetric detection
        ("Cl.[O-]C(=O)CC([O-])=O.CCCC1=CNC=C1CCC", True),
        # regular salt, should be symmetrical
        ("Cl.CC([O-])=O.CCCC1=CNC=C1CCC", True),
        # not a salt, should not be symmetrical
        ("[O-]C(=O)C1=CC=CC=C1.CCCC1=CNC=C1CCC", False),
    ]

    mols, expected_vals = zip(*test_mols)

    scores = [score_symmetry(x, exclude_self_mapped_edged=False) for x in mols]
    fully_symmetrical = [x == 1 for x in scores]

    assert list(expected_vals) == fully_symmetrical

    # if we don't exclude self_mapped edges, then we will have non symmetrical here
    expected_mol2 = score_symmetry(mols[1], exclude_self_mapped_edged=True)
    assert expected_mol2 < 1
    assert expected_mol2 > 0.8
