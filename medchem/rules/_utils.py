_DESCRIPTOR_LIST = [
    "mw",
    "fsp3",
    "n_lipinski_hba",
    "n_lipinski_hbd",
    "n_rings",
    "n_hetero_atoms",
    "n_heavy_atoms",
    "n_rotatable_bonds",
    "n_radical_electrons",
    "n_NHOH",
    "n_NO",
    "tpsa",
    "qed",
    "clogp",
    "sas",
    "formal_charge",
    "refractivity",
    "n_aliphatic_carbocycles",
    "n_aliphatic_heterocyles",
    "n_aliphatic_rings",
    "n_aromatic_carbocycles",
    "n_aromatic_heterocyles",
    "n_aromatic_rings",
    "n_saturated_carbocycles",
    "n_saturated_heterocyles",
    "n_saturated_rings",
    "n_charged_atoms",
    "n_stereo_center",
    "n_aromatic_atoms",
    "n_rigid_bonds",
    "n_aromatic_atoms_proportion",
]


def _in_range(x, min_val: float = -float("inf"), max_val: float = float("inf")):
    """Check if a value is in a range
    Args:
        x: value to check
        min_val: minimum value
        max_val: maximum value
    """
    return min_val <= x <= max_val


def list_descriptors():
    """List all descriptors available for computation"""
    # EN: this is until datamol release
    return _DESCRIPTOR_LIST
