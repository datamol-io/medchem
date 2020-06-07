Lilly-Medchem-Rules
===================


# Preamble
This is an implementation of Eli Lilly Medchem Rules published under "Rules for Identifying Potentially Reactive or Promiscuous Compounds" by Robert F. Bruns and Ian W. Watson, J. Med. Chem. 2012, 55, 9763--9772 as ACS Author choice, i.e. open access at [doi 10.1021/jm301008n](https://doi.org/10.1021/jm301008n).


To quote the abstract:

```
"[This approach] describes a set of 275 rules, developed over an 18-year period, used to identify compounds that may interfere with biological assays, allowing their removal from screening sets. Reasons for rejection include reactivity (e.g., acyl halides), interference with assay measurements (fluorescence, absorbance, quenching), activities that damage proteins (oxidizers, detergents), instability (e.g., latent aldehydes), and lack of druggability (e.g., compounds lacking both oxygen and nitrogen). The structural queries were profiled for frequency of occurrence in druglike and nondruglike compound sets and were extensively reviewed by a panel of experienced medicinal chemists. As a means of profiling the rules and as a filter in its own right, an index of biological promiscuity was developed. The 584 gene targets with screening data at Lilly were assigned to 17 subfamilies, and the number of subfamilies at which a compound was active was used as a promiscuity index."
```

# Installation

## Requirements
This package requires : `gcc` and `g++` for compilation. Use your OS package manager or conda:

```bash
conda install -c conda-forge c-compiler cxx-compiler
# conda install gcc_linux-64 
# conda install gxx_linux-64
```
## Source
Clone the repo and install it locally
```bash
git clone https://bitbucket.org/invivoai/medchem.git
cd medchem 
pip install . # Alternatively you can install a develop version
```
This will build (compile the C source) and install the package. If you are having trouble with pip, use setup.py:

```bash
python setup.py install # "python setup.py build" should not be necessary

```

## pip
You can install from our pypi repo on anaconda

```bash
pip install --index-url https://pypi.anaconda.org/t/$TOKEN/invivoai/medchem
``` 

# Running

## Command line
You can use the provided binary: ```chemfilter --help```
```
Usage: chemfilter [OPTIONS]

Options:
  -n, --min-atoms INTEGER        Use lower atom count cutoff
  -Cs, --soft-max-atoms INTEGER  Apply soft upper atom count cuttof
  -Ch, --hard-max-atoms INTEGER  Apply hard upper atom count cuttof
  --min-num-rings INTEGER        Minimum number of rings accepted
  --max-num-rings INTEGER        Maximum number of rings accepted
  --max-size-rings INTEGER       Maximum ring size (number of atom per ring)
  --max-size-chain INTEGER       Threshold for long carbon chain (7)
  --smarts TEXT                  Optional smarts to reject. File or smarts
                                 string

  --nodemerit                    Use hard rejections only, do not apply any
                                 demerits

  --dthresh INTEGER              Demerit threshold. For relaxed rules, use 160
                                 demerit cutoff

  -o, --output TEXT              output file where to write result  [required]
  --odm TEXT                     Optional demerits to omit to apply. File or
                                 demerit

  --okiso                        Allow isotopic atoms to pass through
  --noapdm                       Do not append demerit reasons
  --smcol TEXT                   Name of the smiles columns
  --allow-non-interesting        Allow molecules with non interesting atoms
                                 only to pass

  -i, --input-file TEXT          Input csv or smi files. Header expected and
                                 first column should always be the smiles if
                                 smiles column name is not provided.
                                 [required]

  --help                         Show this message and exit.
```

## Module
You can import the package instead 
```python
from medchem.filter import score
test_config = {
    'output': 'test',
    'min_atoms': 7,
    'soft_max_atoms': 30,
    'hard_max_atoms': 50,
    'smarts': [],
    'nodemerit': False,
    'dthresh': 160,
    'odm': [],
    'okiso': False,
    'noapdm': False
}
smiles_list = [
    'Cc1cnc(CNc(cccc2-c3cn(CC(C4)CC4O)c4ncnc(N)c34)c2F)s1'
    'Cc1cnc(CNc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)s1',
    'Cc1cnc(CNc(cccc2-c3cn(C[C@H](C4)C[C@H]4O)c4ncnc(N)c34)c2F)s1',
    'Cc1cnc(CNc(cccc2-c3cn(C[C@H]4C[C@H](CO)C4)c4ncnc(N)c34)c2F)s1',
    'Cc1cnc(CNc(cccc2-c3cn(C[C@H]4C[C@@H](CO)C4)c4ncnc(N)c34)c2F)s1',
    'Cc1cnc(CNc2c(CO)c(-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)ccc2)s1',
    'CC(C)(C)c1cnc(CNc(cccc2-c3cn(CC(C4)CC4O)c4ncnc(N)c34)c2F)s1',
    'Cc1cnc(CNc2nc(-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)ccc2)s1',
    'Cc1nc(CNc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)cs1',
    'NC(c(ccc(NC(c1cccc(C(F)(F)F)c1)=O)c1)c1-c1cc2cnc(NC3CC3)nc2cc1)=O',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CC5(COCC5)OCC4)nn3)c2ncn1',
    'C[C@@H](c1cn(CC2CCCCC2)nn1)N(Cc1nccs1)c(cccc1-c2c[nH]c3ncnc(N)c23)c1F',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CC4C(CC5)CC5C4)nn3)c2ncn1',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CC5(CSCC5)OCC4)nn3)c2ncn1',
    'CC1OC(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)CC1',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C[C@@H](C4)c5c4cccc5)nn3)c2ncn1',
    'C[C@@H](c1cn(Cc2nnnn2C2CC2)nn1)Nc1c(c(-c2cccc(NCc3nccs3)c2F)c[nH]2)c2ncn1',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCN4CCOCC4)nn3)c2ncn1',
    'Cc1cnc(C(Nc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)=O)s1',
    'Nc1c(c(-c(cccc2NCc3nccs3)c2F)cn2Cc3cn(-c(cc4)cc5c4scc5)nn3)c2ncn1',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCC4CC4)nn3)c2ncn1',
    'CC(C)Cn1ncnc1Cn1nnc(CNc2c(c(-c3cccc(NCc4nccs4)c3F)c[nH]3)c3ncn2)c1',
    'Nc1c(c(-c2cccc(NCc3nc(CO)cs3)c2F)cn2CC3CCC3)c2ncn1',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CCCC4)nn3)c2ncn1'
]
res = run_scorer(smiles_list, **test_config)
print(res)
```