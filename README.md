## Install

Requires python >= 3.11

Dependencies:

- pandas
- numpy
- plotly
- dash
- dash_cytoscape==0.3.0
- dash-bootstrap-components==1.6.0

## Run in browser

1. Install application -- from project root directory, run pip install .

2. Ensure your CSV has the following columns (case insensitive)

   \*\* asterisks indicate required columns

- \*Path (ligand, em, receptor, target columns are generated from the Path column)
- \*Sender
- \*Receiver
- \*sigprob_x, where x can be any suffix representing experimental group
- \*sigprob_y, where y can be any suffix representing experimental group
- p_value_x, where x can be any suffix representing experimental group
- p_value_y, where y can be any suffix representing experimental group
- tprs
- prs
- umap1 (required for umap display)
- umap2(required for umap display)

3. Open run.sh and update --clusters and --pathways to appropriate paths

Example run.sh:

```
incytr-viz --clusters <clusters path> --pathways <pathways path>

```

4. Run 'bash run.sh' and navigate to http://127.0.0.1:8000/ in your browser

## Run in jupyter notebook (beta)

1. Update filepaths in notebook.ipynb
2. Run notebook

## Notes

Warning! Dash cytoscape > 0.3.0 (e.g. 1.0.0 and above) have bugs that affect graph updates as of July 2024
