# from typing import Iterable
# from typing import List
# from typing import Optional
# from typing import Union

# import os
# import copy
# import functools
# import pandas as pd
# import datamol as dm

# from tqdm.auto import tqdm
# from loguru import logger
# from rdkit.Chem import rdchem
# from rdkit.Chem.Descriptors import MolWt, MolLogP, NumHDonors, NumHAcceptors, TPSA
# from medchem.utils.loader import get_data_path
# from medchem.catalog import NamedCatalogs


# class NovartisFilters:
#     """Filtering class for building a screening deck following the novartis filtering process
#     published in https://dx.doi.org/10.1021/acs.jmedchem.0c01332.

#     The output of the filter are explained below:
#     - **status**: one of `["Exclude", "Flag", "Annotations", "Ok"]` (ordered by quality).
#         Generally, you can keep anything without the "Exclude" label, as long as you also apply
#         a maximum severity score for compounds that collects too many flags.
#     - **covalent**: number of potentially covalent motifs contained in the compound
#     - **severity**: how severe are the issues with the molecules:
#         - `0`: compound has no flags, might have annotations;
#         - `1-9`:  number of flags the compound raises;
#         - `>= 10`:  default exclusion criterion used in the paper
#     - **special_mol**: whether the compound/parts of the compound belongs to a special class of molecules
#         (e.g peptides, glycosides, fatty acid). In that case, you should review the rejection reasons.
#     """

#     def __call__(
#         self,
#         mols: Iterable[Union[str, rdchem.Mol]],
#         n_jobs: Optional[int] = None,
#         progress: bool = False,
#     ):
#         """Run alert evaluation on this list of molecule and return the full dataframe

#         Args:
#             mols: input list of molecules
#             n_jobs: number of jobs
#             progress: whether to show progress or not
#         """

#         catalog = NamedCatalogs.nibr()
#         if n_jobs is not None:
#             if isinstance(mols[0], str):
#                 mols = dm.parallelized(dm.to_mol, mols, n_jobs=n_jobs, progress=progress)
#             matches = dm.parallelized(
#                 catalog.GetMatches,
#                 mols,
#                 n_jobs=n_jobs,
#                 progress=progress,
#                 scheduler="threads",
#             )
#         else:
#             mols = [dm.to_mol(x) if isinstance(x, str) else x for x in mols]
#             iter_mols = mols
#             if progress:
#                 iter_mols = tqdm(mols)
#             matches = [catalog.GetMatches(mol) for mol in iter_mols]

#         results = []
#         for i, (mol, entries) in enumerate(zip(mols, matches)):
#             status = "Ok"
#             smiles = None
#             reasons = None
#             co = None
#             sm = None
#             sc = 0
#             try:
#                 smiles = dm.to_smiles(mol)
#                 if len(list(entries)):
#                     # initialize empty lists
#                     names, severity, covalent, special_mol = ([] for _ in range(4))
#                     # get the matches
#                     for entry in entries:
#                         pname = entry.GetDescription()
#                         _, name, sev, cov, m = pname.split("||")
#                         names.append(name)
#                         severity.append(int(sev))
#                         covalent.append(int(cov))
#                         special_mol.append(int(m))
#                     # concatenate all matching filters
#                     reasons = "; ".join(names)
#                     # severity of 2 means EXCLUDE
#                     if severity.count(2):
#                         sc = 10
#                         status = "Exclude"
#                     else:
#                         sc = sum(severity)
#                         if severity.count(1):
#                             status = "Flag"
#                         elif severity.count(0):
#                             status = "Annotations"
#                     # get number of covalent flags and special molecule flags
#                     co = sum(covalent)
#                     sm = sum(special_mol)
#             except Exception as e:
#                 logger.warning(f"Fail on molecule at index {i}")

#             results.append([smiles, status, reasons, sc, co, sm])
#         df = pd.DataFrame(
#             results,
#             columns=[
#                 "_smiles",
#                 "status",
#                 "reasons",
#                 "severity",
#                 "covalent",
#                 "special_mol",
#             ],
#         )
#         return df

