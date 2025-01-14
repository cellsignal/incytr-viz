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

1. Install dependencies if necessary (pip install -r requirements.txt)
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

3. Open run.sh and update --group_a_populations (exp group) --group_b_populations (wt group) and --pathways to appropriate file paths (see cluster file naming convention must match group names in sigprob columns)

Example run.sh:

```
python app.py --group_a_populations data/covid/BL_clusters.csv --group_b_populations data/covid/HC_clusters.csv --pathways data/covid/hc_bl_incytr_heginput_p_rnascore_ligand-target.tsv

BL_clusters.csv --> sigprob_BL
HC_clusters.csv --> sigprob_HC

```

4. Run 'bash run.sh' and navigate to http://127.0.0.1:8050/ in your browser

## Run in jupyter notebook (beta)

1. Update filepaths in notebook.ipynb
2. Run notebook

## Notes

Warning! Dash cytoscape > 0.3.0 (e.g. 1.0.0 and above) have bugs that affect graph updates as of July 2024
