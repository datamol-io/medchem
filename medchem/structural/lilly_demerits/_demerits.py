from typing import Optional
from typing import Union
from typing import Any
from typing import Sequence

import os
import uuid
import re
import shutil
import tempfile
import importlib.resources as importlib_resources

import pandas as pd
import datamol as dm
import numpy as np


from ._lilly import find_lilly_binaries
from ._lilly import parse_output
from ._lilly import run_cmd


BIN2PATH = find_lilly_binaries()


class LillyDemeritsFilters:
    """Lilly MedChem Rules published in:

    [Robert F. Bruns and Ian A. Watson, Rules for Identifying Potentially Reactive or Promiscuous Compounds,
    Journal of Medicinal Chemistry 2012 55 (22), 9763-9772](https://pubs.acs.org/doi/10.1021/jm301008n)


    !!! abstract "Description"
        This is a set of 275 rules, developed over an 18-year period, used to identify compounds that may interfere with biological assays,
        allowing their removal from screening sets. Reasons for rejection include reactivity (e.g., acyl halides),
        interference with assay measurements (fluorescence, absorbance, quenching), activities that damage proteins (oxidizers, detergents),
        instability (e.g., latent aldehydes), and lack of druggability (e.g., compounds lacking both oxygen and nitrogen).


    """

    def __init__(
        self,
        mc_first_pass_options: Optional[str] = None,
        iwd_options: Optional[str] = None,
        stop_after_step: int = 3,
        **run_options: Any,
    ):
        """
        Constructor for the Lilly MedChem Rules

        Args:
            mc_first_pass_options: Initial options to pass to mc_first_pass
            iwd_options: Initial options to pass to iwdemerit
            stop_after_step: Where to stop in the pipeline. Don't change this if you don't know.
            run_options: Additional option to run the pipeline
        """
        self.mc_first_pass_options = mc_first_pass_options
        self.iwd_options = iwd_options
        self.stop_after_step = stop_after_step
        self.run_options = run_options

    def __call__(
        self,
        mols: Sequence[Union[str, dm.Mol]],
        n_jobs: Optional[int] = -1,
        batch_size: int = 5_000,
        progress: bool = False,
        progress_leave: bool = False,
        scheduler: str = "threads",
    ):
        if scheduler == "processes":
            raise ValueError("Demerits filtering does not support processes or loky mode yet.")

        n_splits = max(1, int(np.ceil(len(mols) / batch_size)))

        if n_splits > 1:
            mols_batch_list = np.array_split(mols, n_splits)
            # EN: cannot run this code in processes or loky mode
            results = dm.parallelized(
                self._score,
                mols_batch_list,
                n_jobs=n_jobs,
                progress=progress,
                scheduler=scheduler,
                tqdm_kwargs=dict(
                    desc="Demerits filtering",
                    leave=progress_leave,
                ),
            )
            return pd.concat(results, ignore_index=True)  # type: ignore

        else:
            return self._score(mols)

    def _score(self, mols: Sequence[Union[str, dm.Mol]]):
        """Run lilly medchem scorer on input smile list:

        Args:
            mols: list of smiles
        """

        mc_first_pass_options = self.mc_first_pass_options
        iwd_options = self.iwd_options
        stop_after_step = self.stop_after_step
        run_options = self.run_options

        if mc_first_pass_options is None:
            mc_first_pass_options = ""

        if iwd_options is None:
            iwd_options = ""

        extra_iwdemerit_options = ""
        demerit_cutoff = run_options.get("dthresh", None)
        soft_upper_atom = run_options.get("soft_max_atoms", 30)
        hard_upper_atom = run_options.get("hard_max_atoms", 50)
        max_size_rings = run_options.get("max_size_rings")
        min_num_rings = run_options.get("min_num_rings")
        max_num_rings = run_options.get("max_num_rings")
        max_size_chain = run_options.get("max_size_chain")
        atom_count = run_options.get("min_atoms", 1)
        allow_non_int_atoms = run_options.get("allow_non_interesting", False)

        ring_bond_ratio = run_options.get("ring_bond_ratio", -1)
        okiso = run_options.get("okiso", False)
        if not okiso:
            mc_first_pass_options += " -I 0 "
        if min_num_rings:
            mc_first_pass_options += f"-r {min_num_rings} "
        if max_num_rings:
            mc_first_pass_options += f"-R {max_num_rings} "
        if max_size_rings:
            mc_first_pass_options += f"-Z {max_size_rings} "
            extra_iwdemerit_options += f" -Z {max_size_rings} "
        if max_size_chain:
            extra_iwdemerit_options += f" -z {max_size_chain} "

        if allow_non_int_atoms:
            mc_first_pass_options += "-k "
        mc_first_pass_options += " -A I -A ipp"

        query_files = ["reject1", "reject2", "demerits"]
        for i, query_file in enumerate(query_files):
            query_files[i] = str(importlib_resources.files("medchem.data.queries").joinpath(query_file))

        # output file dir
        files_to_be_deleted = []
        run_id = str(uuid.uuid4())[:8]
        bad_file_dir = tempfile.mkdtemp(suffix=f"_lilly_{run_id}")
        files_to_be_deleted.append(bad_file_dir)
        bad_file_0 = os.path.join(bad_file_dir, "bad0")
        bad_file_1 = os.path.join(bad_file_dir, "bad1")
        bad_file_2 = os.path.join(bad_file_dir, "bad2")
        bad_file_3 = os.path.join(bad_file_dir, "bad3")
        mc_pass_out = os.path.join(bad_file_dir, "mc_pass.smi")
        tsub_out_1 = os.path.join(bad_file_dir, "tsub1.smi")
        tsub_out_2 = os.path.join(bad_file_dir, "tsub2.smi")
        iwd_out = os.path.join(bad_file_dir, "iwd.smi")

        # optional_queries
        optional_queries = ""
        optional_queries += " -q ".join(run_options.get("rej_queries", []))
        optional_queries += " -s ".join(run_options.get("smarts", []))

        # extra iwdemerit options
        nodemerit = run_options.get("nodemerit", False)

        if nodemerit or soft_upper_atom is None:
            soft_upper_atom = hard_upper_atom - 1

        if hard_upper_atom < soft_upper_atom:
            hard_upper_atom = soft_upper_atom + 1

        if demerit_cutoff:
            extra_iwdemerit_options += f" -f {demerit_cutoff} "

        if nodemerit:
            extra_iwdemerit_options += " -r "

        if iwd_options:
            extra_iwdemerit_options += " " + iwd_options

        charge_assigner_path = str(
            importlib_resources.files("medchem.data.charge_assigner").joinpath("queries")
        )
        extra_iwdemerit_options += " -N F:" + charge_assigner_path

        odm = run_options.get("odm", [])

        if odm:
            odm_patterns = [re.compile(odm_p, re.I) for odm_p in odm]
            with open(query_files[2]) as QRY_IN:
                current_queries = [qry.strip() for qry in QRY_IN]
                current_queries = [
                    qry for qry in current_queries if not any(odm_p.search(qry) for odm_p in odm_patterns)
                ]

            with tempfile.NamedTemporaryFile(mode="w+t", suffix=".qry", delete=False) as tmp_file:
                tmp_file.write("\n".join(current_queries))
                query_files[2] = tmp_file.name
                files_to_be_deleted.append(tmp_file.name)

        smiles_file = None
        with tempfile.NamedTemporaryFile(
            mode="w+t",
            suffix=f"_lilly_{run_id}.smi",
            delete=False,
        ) as smiles_tmp_files:
            # Convert the input mols to a list of SMILES
            smiles_list = []
            for mol in mols:
                if isinstance(mol, str):
                    smiles = mol
                else:
                    smiles = dm.to_smiles(mol)

                # Sanity check
                if smiles is None or dm.to_mol(smiles) is None:
                    raise ValueError(f"Invalid SMILES: {smiles}. Demerits does not support invalid mol yet.")

                smiles_list.append(smiles)

            smiles_tmp_files.write(
                "\n".join([f"{sm.strip().split()[0]}\t{i}" for i, sm in enumerate(smiles_list)])
            )
            smiles_file = smiles_tmp_files.name
            files_to_be_deleted.append(smiles_file)

        cmd = [BIN2PATH["mc_first_pass"]]
        if ring_bond_ratio >= 0:
            cmd.extend(["-b", str(ring_bond_ratio)])

        if mc_first_pass_options:
            cmd.extend(mc_first_pass_options.split())

        cmd.extend(["-c", str(atom_count), "-C", str(hard_upper_atom)])
        cmd.extend("-E autocreate -o smi -V -g all -g ltltr -i ICTE".split())
        cmd.extend(["-L", bad_file_0, "-K", "TP1"])
        cmd.extend(["-a", "-u", "-S", mc_pass_out, smiles_file])

        out = []
        out.append(run_cmd(cmd))

        if stop_after_step >= 1:
            cmd = []
            cmd.extend((BIN2PATH["tsubstructure"] + " -E autocreate -b -u -i smi -o smi -A D ").split())
            cmd.extend(("-m " + bad_file_1 + " -m QDT").split())
            cmd.extend(
                (f"-n {tsub_out_1} -q F:" + query_files[0] + optional_queries + f" {mc_pass_out}").split()
            )
            out.append(run_cmd(cmd))

        if stop_after_step >= 2:
            cmd = []
            cmd.extend((BIN2PATH["tsubstructure"] + " -A D -E autocreate -b -u -i smi -o smi ").split())
            cmd.extend(("-m " + bad_file_2 + " -m QDT").split())
            cmd.extend((f"-n {tsub_out_2} -q F:" + query_files[1] + f" {tsub_out_1}").split())
            out.append(run_cmd(cmd))

        if stop_after_step >= 3:
            cmd = []
            cmd.extend(
                (
                    BIN2PATH["iwdemerit"]
                    + " -u -k -x -t "
                    + extra_iwdemerit_options
                    + " -E autocreate -A D -i smi -o smi -q F:"
                    + query_files[2]
                ).split()
            )
            cmd.extend(f"-R {bad_file_3}".split())
            cmd.extend(
                f"-G {iwd_out} -c smax={soft_upper_atom} -c hmax={hard_upper_atom} {tsub_out_2}".split()
            )
            out.append(run_cmd(cmd))

        data_list = []
        for i, bad_file in enumerate([bad_file_0, bad_file_1, bad_file_2]):
            with open(bad_file + ".smi") as IN:
                df_bad = parse_output(IN)
                df_bad["step"] = i + 1
                df_bad["rejected"] = True
                data_list.append(df_bad)

        demerit_extractor = re.compile(r"'([A-Za-z0-9_\s]+)'")
        i = 0
        for dt_file, rej in [(bad_file_3 + ".smi", True), (iwd_out, False)]:
            parseable = []
            demerit_scores = []
            with open(dt_file) as IN:
                for row in IN:
                    in_string, *demerit_string = row.split(" : ")
                    in_string = in_string.strip()
                    demerit_score = 0
                    if not demerit_string:
                        in_string += ' ""'
                    else:
                        demerit_string = str(demerit_string[0]).strip()
                        m = re.match(r"D\(([0-9]+)\)", demerit_string)
                        if m is not None:
                            demerit_score = int(m.group(1))
                        tmp = demerit_extractor.findall(row)
                        demerit_string = ",".join([":".join(reversed(line.split())) for line in tmp])
                        in_string += " " + demerit_string
                    parseable.append(in_string)
                    demerit_scores.append(demerit_score)

            df = parse_output(parseable)
            df["rejected"] = rej
            df["step"] = i + 1
            df["demerit_score"] = demerit_scores
            data_list.append(df)
            i += 1
        results = pd.concat(data_list).sort_values("ID").reset_index(drop=True)

        # Postprocessing
        results["status"] = results["rejected"].apply(lambda x: "Exclude" if x else "Ok")
        results.loc[((results.demerit_score > 0) & (~results.rejected)), "status"] = "Flag"

        results["status"] = results["status"].str.lower()
        results["pass_filter"] = ~results["rejected"]
        results = results.drop(columns=["ID", "rejected"])
        results["mol"] = mols

        # clean
        for to_del in files_to_be_deleted:
            if os.path.isfile(to_del):
                os.remove(to_del)
            else:
                shutil.rmtree(to_del)

        return results
