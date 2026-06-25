# SSN-from-PFam

Construct a protein sequence similarity network from a UniProt family query, using BLASTP for EFI-like similarity and Graphviz/NetworkX for the rendered network.

Install Python dependencies with `python3 -m pip install -r requirements.txt`.

## What it does

1. Fetches a UniProt family or InterPro/Pfam query.
2. Filters the resulting protein set.
3. Runs all-vs-all BLASTP locally using the BLAST+ installation on your machine.
4. Builds the network with NetworkX.
5. Exports SVG, XGMML, GraphML, DOT, and TSV summaries.

## Install

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
* `requirements.txt` lists the Python dependency.
* `networkx` must be installed by the user with pip.
* NCBI BLAST+ must be installed locally with Homebrew; the workflow looks for `blastp` and `makeblastdb` on your PATH.

