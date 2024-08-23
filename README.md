## Install

Requires python >= 3.11

Dependencies:

- pandas
- numpy
- plotly
- dash
- dash_cytoscape==0.3.0

## Run in browser

1. Install dependencies if necessary (pip install -r requirements.txt)
2. Ensure your CSV has the following columns (currently case sensitive):

- Path (ligand, em, receptor, target columns are generated from the Path column)
- Ligand
- Receptor
- EM
- Target
- final_score
- SigWeight_EXP
- SigWeight_WT
- RNA_score
- Sender
- Receiver
- adjlog2FC

3. Open main.py and update CLUSTERS_A_FILE (exp group) CLUSTERS_B_FILE (wt group) and PATHWAYS_FILE to appropriate file paths. Examples in docs/ folder
4. Run 'python main.py'

## Run in jupyter notebook

1. See steps 2-3 above
2. Run notebook

## Notes

Warning! Dash cytoscape > 0.3.0 (e.g. 1.0.0 and above) have bugs that affect graph updates as of July 2024
