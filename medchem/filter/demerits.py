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
    """Replace last occurence of 'old' in 'input_str' by 'rep' """
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
    smiles_list,
    mc_first_pass_options="",
    iwd_options="",
    stop_after_step=3,
    **run_options,
):
    """Run scorer on input smile list:

    Args:
        smiles_list: list of smiles
        mc_first_pass_options: str
            Initial options to pass to mc_first_pass
        iwd_options: str
            Initial options to pass to iwdemerit
        stop_after_step: int, optional
            Where to stop in the pipeline. Don't change this if you don't know.
        run_options: dict
            Additional option to run the pipeline

    Returns:
        out_df: pd.DataFrame
            Dataframe containing the smiles and computed properties: (rejected, demerit_score, reason, step)
    """

    extra_iwdemerit_options = ""
    demerit_cutoff = run_options.get("dthresh", None)
    soft_upper_atom = run_options.get("soft_max_atoms")
    hard_upper_atom = run_options.get("hard_max_atoms", 50)
    max_size_rings = run_options.get("max_size_rings")
    min_num_rings = run_options.get("min_num_rings")
    max_num_rings = run_options.get("max_num_rings")
    max_size_chain = run_options.get("max_size_chain")
    atom_count = run_options.get("min_atoms")
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


if __name__ == "__main__":
    test_config = {
        "output": "test",
        "min_atoms": 7,
        "soft_max_atoms": 30,
        "hard_max_atoms": 50,
        "smarts": [],
        "nodemerit": False,
        "dthresh": 160,
        "odm": [],
        "okiso": False,
        "noapdm": False,
    }
    smiles_list = [
        "Cc1cnc(CNc(cccc2-c3cn(CC(C4)CC4O)c4ncnc(N)c34)c2F)s1"
        "Cc1cnc(CNc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc(cccc2-c3cn(C[C@H](C4)C[C@H]4O)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc(cccc2-c3cn(C[C@H]4C[C@H](CO)C4)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc(cccc2-c3cn(C[C@H]4C[C@@H](CO)C4)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc2c(CO)c(-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)ccc2)s1",
        "CC(C)(C)c1cnc(CNc(cccc2-c3cn(CC(C4)CC4O)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc2nc(-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)ccc2)s1",
        "Cc1nc(CNc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)cs1",
        "NC(c(ccc(NC(c1cccc(C(F)(F)F)c1)=O)c1)c1-c1cc2cnc(NC3CC3)nc2cc1)=O",
        "Cc1cnc(CNc(cccc2-c3cn(Cc4cn(CCN5CCOCC5)nn4)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc2cc(-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)ccc2)s1",
        "Cc1cnc(CNc(cccc2-c3cn(CC4CCOCC4)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc(cccc2-c3cn(C[C@@H]4COCC4)c4ncnc(N)c34)c2F)s1",
        "Cn1nnc(CNc2c(c(-c3cccc(NCc4nccs4)c3F)c[nH]3)c3ncn2)c1",
        "Cn1c(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)cnc1",
        "Cc1cnc(C[C@@H]2NCc3c2cccc3-c2cn(C[C@H](C3)C[C@@H]3O)c3ncnc(N)c23)s1",
        "Cc1ncc(Cn2nnc(CNc3ncnc4c3c(-c(cccc3NCc5nccs5)c3F)c[nH]4)c2)s1",
        "CCc1noc(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)n1",
        "CCn1c(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)cnc1",
        "Fc1cccc(-n2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)c1",
        "CCc1nc(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)no1",
        "Cc(n1CCn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)ncc1[N+]([O-])=O",
        "Fc(cc1F)cc(Br)c1-n1nnc(CNc2ncnc3c2c(-c(cccc2NCc4nccs4)c2F)c[nH]3)c1",
        "Fc1cc(-n2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)cc(F)c1",
        "Cc1cnc(CNc(cc2)cc(-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)s1",
        "C[C@@H](c1cn(-c2cc(C(F)(F)F)cc(Br)c2)nn1)Nc1ncnc2c1c(-c(cccc1NCc3nccs3)c1F)c[nH]2",
        "Nc1c(c(-c2cccc(NCc3cscn3)c2F)cn2C[C@H](C3)C[C@@H]3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(Cc4cncc(F)c4)nn3)c2ncn1",
        "C[C@@H](c1cn(Cc2cncn2C2CC2)nn1)Nc1c(c(-c2cccc(NCc3nccs3)c2F)c[nH]2)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3cscn3)c2F)cn2C[C@H](C3)C[C@H]3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(-c4cc([N+]([O-])=O)ccc4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C[C@@H]4OCCc5c4cccc5)nn3)c2ncn1",
        "CC(c1cn(CC2C(CC3)CC3C2)nn1)N(Cc1nccs1)c(cccc1-c2c[nH]c3ncnc(N)c23)c1F",
        "Nc1c(c(-c2cccc(NCc3ncc4n3CCC4)c2F)cn2C[C@H](C3)C[C@@H]3O)c2ncn1",
        "NC1CC(Cn2c3ncnc(N)c3c(-c3cccc(NCc4nc(CO)cs4)c3F)c2)C1",
        "Nc1c(c(-c2cccc(NCc3ncccc3)c2F)nn2CC(C3)CC3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CC5(COCC5)OCC4)nn3)c2ncn1",
        "C[C@@H](c1cn(CC2CCCCC2)nn1)N(Cc1nccs1)c(cccc1-c2c[nH]c3ncnc(N)c23)c1F",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CC4C(CC5)CC5C4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CC5(CSCC5)OCC4)nn3)c2ncn1",
        "CCn1ncnc1Cn1nnc(CNc2c(c(-c3cccc(NCc4nccs4)c3F)c[nH]3)c3ncn2)c1",
        "Nc1c(c(-c(cccc2NCc3nccs3)c2F)cn2Cc3cn([C@@H]4CC5(CCSCC5)OCC4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn([C@@H]4CC5(CCCCC5)OCC4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn([C@@H]4CC5(CCC5)OCC4)nn3)c2ncn1",
        "CCn1c(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)ncc1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn([C@@H]4CC5(CCOCC5)OCC4)nn3)c2ncn1",
        "CC1OC(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)CC1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C[C@@H](C4)c5c4cccc5)nn3)c2ncn1",
        "C[C@@H](c1cn(Cc2nnnn2C2CC2)nn1)Nc1c(c(-c2cccc(NCc3nccs3)c2F)c[nH]2)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCN4CCOCC4)nn3)c2ncn1",
        "Cc1cnc(C(Nc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)=O)s1",
        "Nc1c(c(-c(cccc2NCc3nccs3)c2F)cn2Cc3cn(-c(cc4)cc5c4scc5)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCC4CC4)nn3)c2ncn1",
        "CC(C)Cn1ncnc1Cn1nnc(CNc2c(c(-c3cccc(NCc4nccs4)c3F)c[nH]3)c3ncn2)c1",
        "Nc1c(c(-c2cccc(NCc3nc(CO)cs3)c2F)cn2CC3CCC3)c2ncn1",
        "CC(c1cn(C2CC3(CCCC3)OCC2)nn1)Nc1c(c(-c2cccc(NCc3nccs3)c2F)c[nH]2)c2ncn1",
        "Nc1c(c(-c2cccc(NCC(C(NC=C3)=O)=C3c3ncccc3)c2F)cn2C[C@H](C3)C[C@@H]3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CCOCC4)nn3)c2ncn1",
        "CCN1CC(Cn2nnc(C(C)Nc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)OCC1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(-c4cc(F)cc(F)c4)nn3)c2ncn1",
        "CC(c1cn(C2CC3(CCOCC3)OCC2)nn1)Nc1c(c(-c2cccc(NCc3nccs3)c2F)c[nH]2)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CCCCCC4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCN4CCCC4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCC(C(NC=C3)=O)=C3c3ccccc3)c2F)cn2C[C@H](C3)C[C@@H]3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CCCCC4)nn3)c2ncn1",
        "Cc(n1CCn2nnc(Cn3c4ncnc(N)c4c(-c4cccc(NCc5nccs5)c4F)c3)c2)ncc1[N+]([O-])=O",
        "Nc1c(c(-c2cccc(NCc3ncc(-c4ccccc4)s3)c2F)cn2C[C@H](C3)C[C@@H]3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CCCC4)nn3)c2ncn1",
    ]
    res = score(smiles_list, **test_config)
    print(res)
