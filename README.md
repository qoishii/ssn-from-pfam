# SSN-from-PFam

Construct a protein sequence similarity network from a UniProt family query, using BLASTP for EFI-like similarity and Graphviz/NetworkX for the rendered network.

Install dependencies either with conda or with Homebrew plus pip.

## What it does

1. Fetches a UniProt family or InterPro/Pfam query.
2. Filters the resulting protein set.
3. Runs all-vs-all BLASTP locally using the BLAST+ installation on your machine.
4. Builds the network with NetworkX.
5. Exports SVG, XGMML, GraphML, DOT, and TSV summaries.

## Install

### Option 1: conda

Use this option if you want all project dependencies isolated inside a dedicated conda environment. This is helpful for bioinformatics users who work across multiple projects, shared computers, clusters, or systems where installing command-line tools globally is inconvenient.

The `environment.yml` file is a conda environment recipe. It belongs in the root of this repository, beside `README.md`, and should be uploaded to GitHub with the rest of the project. Users can recreate the expected software environment from it:

```bash
conda env create -f environment.yml
conda activate ssn-from-pfam
```

This option is not especially relevant if you do not use conda, if you already manage dependencies with Homebrew/pip, or if your institution provides BLAST+ and Graphviz through another module system.

### Option 2: Homebrew plus pip

```bash
brew install python blast graphviz
python3 -m pip install -r requirements.txt
```

## Run

```bash
python3 work/run_ssn_workflow.py \
  --query 'xref:interpro-IPR024001' \
  --name IPR024001
```

## Example output

The reference IPR024001 outputs are in `outputs/uniprot_IPR024001_network/`.

## Files

* `work/` contains the scripts.
* `outputs/uniprot_IPR024001_network/` contains the final example outputs.
* `environment.yml` defines an optional conda environment with Python, BLAST+, Graphviz, and NetworkX.
* `requirements.txt` lists the Python dependency.
* `networkx` must be installed by the user, either through conda or pip.
* NCBI BLAST+ must be installed locally through conda, Homebrew, or another system package manager; the workflow looks for `blastp` and `makeblastdb` on your PATH.

