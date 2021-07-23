from typing import List
from typing import Optional

import glob
import os
import pandas as pd
import re
import shutil
import subprocess
import tempfile
from io import StringIO

this_dir, _ = os.path.split(__file__)
BUILD_DIR = os.path.join(this_dir, "lilly/build")
QUERY_DIR = os.path.join(this_dir, "data/queries")
CHARGE_ASSIGN_DIR = os.path.join(this_dir, "data/charge_assigner")
BIN2PATH = dict((x, os.path.join(BUILD_DIR, x)) for x in os.listdir(BUILD_DIR))


def _rreplace(input_str, old, rep, occurrence):
    """Replace last occurence of 'old' in 'input_str' by 'rep'"""
    tmp = input_str.rsplit(old, occurrence)
    return rep.join(tmp)


def _parse_output(rowlist):
    """Parse content of `rowlist` to dataframe"""
    content = "\n".join(
        [
            _rreplace(
                re.sub(r"\s+\([0-9]*\s+(matches\sto\s)", ' "', line.strip(), 1),
                "')",
                "'\"",
                1,
            ).strip("'")
            for line in rowlist
        ]
    )
    flux = StringIO(content)
    df = pd.read_csv(
        flux, sep="\s+", doublequote=True, names=["_smiles", "ID", "reasons"]
    )
    df["reasons"] = df["reasons"].apply(
        lambda x: x.strip("'") if x and isinstance(x, str) else x
    )
    return df


def run_cmd(cmd, shell=False):
    """Run command"""
    res = subprocess.run(cmd, capture_output=True, shell=shell, check=False)
    if res.returncode != 0:
        print("".join(res.stderr.decode("utf-8")))
        print(" ".join(cmd))
        res.check_returncode()
    return res


def score(
    smiles_list: List,
    mc_first_pass_options: Optional[str] = "",
    iwd_options: Optional[str] = "",
    stop_after_step: Optional[int] = 3,
    **run_options,
):
    """Run scorer on input smile list:

    Args:
        smiles_list: list of smiles
        mc_first_pass_options: Initial options to pass to mc_first_pass
        iwd_options: Initial options to pass to iwdemerit
        stop_after_step: Where to stop in the pipeline. Don't change this if you don't know.
        run_options: Additional option to run the pipeline

    Returns:
        out_df (pd.DataFrame): Dataframe containing the smiles and computed properties:
            (rejected, demerit_score, reason, step)
    """

    extra_iwdemerit_options = ""
    demerit_cutoff = run_options.get("dthresh", None)
    soft_upper_atom = run_options.get("soft_max_atoms")
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
    if not os.path.isdir(QUERY_DIR):
        raise ValueError(
            "Query dir not found ! Make sure package is properly installed."
        )
    for i, qf in enumerate(query_files):
        if not os.path.isfile(os.path.join(QUERY_DIR, qf)):
            raise ValueError(
                f"{qf} not found ! Make sure package is properly installed."
            )
        query_files[i] = os.path.join(QUERY_DIR, qf)

    # output file dir
    files_to_be_deleted = []
    bad_file_dir = tempfile.mkdtemp(suffix="lilly")
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
    if os.path.isfile(os.path.join(CHARGE_ASSIGN_DIR, "queries")):
        extra_iwdemerit_options += " -N F:" + os.path.join(CHARGE_ASSIGN_DIR, "queries")

    odm = run_options.get("odm", [])

    if odm:
        odm_patterns = [re.compile(odm_p, re.I) for odm_p in odm]
        with open(query_files[2]) as QRY_IN:
            current_queries = [qry.strip() for qry in QRY_IN]
            omit_queries = [
                qry
                for odm_p in odm_patterns
                for qry in set(current_queries)
                if odm_p.search(qry) and not current_queries.remove(qry)
            ]

        with tempfile.NamedTemporaryFile(
            mode="w+t", suffix=".qry", delete=False
        ) as tmp_file:
            tmp_file.write("\n".join(current_queries))
            query_files[2] = tmp_file.name
            files_to_be_deleted.append(tmp_file.name)

    smiles_file = None
    with tempfile.NamedTemporaryFile(
        mode="w+t", suffix=".smi", delete=False
    ) as smiles_tmp_files:
        smiles_tmp_files.write(
            "\n".join(
                [f"{sm.strip().split()[0]}\t{i}" for i, sm in enumerate(smiles_list)]
            )
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
        cmd.extend(
            (
                BIN2PATH["tsubstructure"] + " -E autocreate -b -u -i smi -o smi -A D "
            ).split()
        )
        cmd.extend(("-m " + bad_file_1 + " -m QDT").split())
        cmd.extend(
            (
                f"-n {tsub_out_1} -q F:"
                + query_files[0]
                + optional_queries
                + f" {mc_pass_out}"
            ).split()
        )
        out.append(run_cmd(cmd))
    if stop_after_step >= 2:
        cmd = []
        cmd.extend(
            (
                BIN2PATH["tsubstructure"] + " -A D -E autocreate -b -u -i smi -o smi "
            ).split()
        )
        cmd.extend(("-m " + bad_file_2 + " -m QDT").split())
        cmd.extend(
            (f"-n {tsub_out_2} -q F:" + query_files[1] + f" {tsub_out_1}").split()
        )
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
            df_bad = _parse_output(IN)
            df_bad["step"] = i + 1
            df_bad["rejected"] = True
            data_list.append(df_bad)

    demerit_extractor = re.compile(r"'([A-Za-z0-9_\s]+)'")
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
                    try:
                        demerit_score = int(
                            re.match(r"D\(([0-9]+)\)", demerit_string).group(1)
                        )
                    except:
                        pass
                    tmp = demerit_extractor.findall(row)
                    demerit_string = ",".join(
                        [":".join(reversed(l.split())) for l in tmp]
                    )
                    in_string += " " + demerit_string
                parseable.append(in_string)
                demerit_scores.append(demerit_score)
        df = _parse_output(parseable)
        df["rejected"] = rej
        df["step"] = i + 1
        df["demerit_score"] = demerit_scores
        data_list.append(df)
        i += 1
    final_df = pd.concat(data_list).sort_values("ID").reset_index(drop=True)
    final_df["status"] = final_df["rejected"].apply(lambda x: "Exclude" if x else "Ok")
    final_df.loc[
        ((final_df.demerit_score > 0) & (~final_df.rejected)), "status"
    ] = "Flag"
    # clean
    for to_del in files_to_be_deleted:
        if os.path.isfile(to_del):
            os.remove(to_del)
        else:
            try:
                shutil.rmtree(to_del)
            except:
                pass
    return final_df
