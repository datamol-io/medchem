#!/usr/bin/env python
import click
import os
import pandas as pd
from medchem.demerits import score
from medchem.filter import lead


@click.command(
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.pass_context
@click.option(
    "--min-atoms", "-n", type=int, help="Use lower atom count cutoff", default=7
)
@click.option(
    "--soft-max-atoms",
    "-Cs",
    type=int,
    help="Apply soft upper atom count cuttof",
    default=30,
)
@click.option(
    "--hard-max-atoms",
    "-Ch",
    type=int,
    help="Apply hard upper atom count cuttof",
    default=50,
)
@click.option("--min-num-rings", type=int, help="Minimum number of rings accepted")
@click.option("--max-num-rings", type=int, help="Maximum number of rings accepted")
@click.option(
    "--max-size-rings", type=int, help="Maximum ring size (number of atom per ring)"
)
@click.option(
    "--max-size-chain", type=int, default=7, help="Threshold for long carbon chain (7)"
)
@click.option(
    "--smarts", multiple=True, help="Optional smarts to reject. File or smarts string"
)
@click.option(
    "--nodemerit",
    is_flag=True,
    help="Use hard rejections only, do not apply any demerits",
)
@click.option(
    "--dthresh",
    type=int,
    help="Demerit threshold. For relaxed rules, use 160 demerit cutoff",
)
@click.option("--output", "-o", required=True, help="output file where to write result")
@click.option(
    "--odm", multiple=True, help="Optional demerits to omit to apply. File or demerit"
)
@click.option("--okiso", is_flag=True, help="Allow isotopic atoms to pass through")
@click.option("--noapdm", is_flag=True, help="Do not append demerit reasons")
@click.option("--smcol", help="Name of the smiles columns")
@click.option(
    "--allow-non-interesting",
    is_flag=True,
    help="Allow molecules with non interesting atoms only to pass",
)
@click.option(
    "--input-file",
    "-i",
    required=True,
    help="Input csv or smi files. Header expected and first column should always be the smiles if smiles column name is not provided.",
)
@click.option(
    "--alert-filter",
    is_flag=True,
    help="Whether to run lead filtering (alerts) on the molecules",
)
@click.option(
    "--alerts",
    multiple=True,
    default=["BMS"],
    help="alerts to use for lead filtering",
)
def run(
    ctx,
    min_atoms,
    soft_max_atoms,
    hard_max_atoms,
    min_num_rings,
    max_num_rings,
    max_size_rings,
    max_size_chain,
    smarts,
    nodemerit,
    dthresh,
    output,
    odm,
    okiso,
    noapdm,
    smcol,
    allow_non_interesting,
    input_file,
    alert_filter,
    alerts,
):
    if (hard_max_atoms and soft_max_atoms) and hard_max_atoms < soft_max_atoms:
        raise ValueError("--hard_max_atoms should be greater than --soft_max_atoms")

    def _populate(params):
        if isinstance(params, str):
            params = [params]
        res = []
        for p in params:
            if os.path.isfile(p):
                res.extend([x.strip() for x in open(p)])
            else:
                res.append(p)
        return res

    ctx.params["smarts"] = _populate(smarts)
    ctx.params["odm"] = _populate(odm)
    ctx.params.pop("input_file", None)
    ctx.params.pop("smcol", None)
    df = pd.read_csv(input_file)
    if smcol:
        smiles_list = df[smcol].values
    else:
        smiles_list = df.ix[:, 0].values
    if alert_filter:
        df["leadlike"] = lead.alert_filter(
            smiles_list, alerts=alerts, n_jobs=os.cpu_count()
        )
    results = score(smiles_list, **ctx.params)

    cols = ["rejected"]
    if not nodemerit:
        cols += ["reasons", "demerit_score"]
    df[cols] = results[cols]
    df.to_csv(output, index=False)


if __name__ == "__main__":
    run()