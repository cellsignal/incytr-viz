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

3. Open run.sh and update --group_a_popluations (exp group) --group_b_populations (wt group) and --pathways to appropriate file paths

4. Run 'bash run.sh' and navigate to http://127.0.0.1:8050/ in your browser

## Run in jupyter notebook

1. Update filepaths in notebook.ipynb
2. Run notebook

## Notes

Warning! Dash cytoscape > 0.3.0 (e.g. 1.0.0 and above) have bugs that affect graph updates as of July 2024
