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
2. Open main.py and update CLUSTERS_FILE and PATHWAYS_FILE to appropriate file paths. Examples in docs/ folder
3. Run 'python main.py'

## Run in jupyter notebook

1. Open notebook.ipynb and update CLUSTERS_FILE and PATHWAYS_FILE to appropriate file paths. Examples in docs/ folder
2. Run notebook

## Notes

Warning! Dash cytoscape > 0.3.0 (e.g. 1.0.0 and above) have bugs that affect graph updates as of July 2024
