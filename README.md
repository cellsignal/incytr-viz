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
2. Ensure your CSV has the following columns (currently case insensitive):

- Path (ligand, em, receptor, target columns are generated from the Path column)
- Sender
- Receiver
- Ligand
- Receptor
- EM
- Target
- SigWeight_x, where x can be any suffix representing experimental group
- SigWeight_y, where y can be any suffix representing experimental group
- p_value_x, where x can be any suffix representing experimental group
- p_value_y, where y can be any suffix representing experimental group
- RNA_score (optional)
- final_score (optional)

3. Open main.py and update CLUSTERS_A_FILE (exp group) CLUSTERS_B_FILE (wt group) and PATHWAYS_FILE to appropriate file paths

4. Run 'python main.py'

## Run in jupyter notebook

1. See steps 2-3 above
2. Run notebook

## Notes

Warning! Dash cytoscape > 0.3.0 (e.g. 1.0.0 and above) have bugs that affect graph updates as of July 2024
