# Incytr Visualization

## Install

Requires Python >= 3.10

*It is recommended to install this program in a virtual environment using tools such as [venv](https://docs.python.org/3/library/venv.html) or [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)


```pip install --upgrade incytr-viz```


## Demo

This will download small demo files from [Zenodo](https://zenodo.org/records/14775408) drawn from the 5XFAD Alzheimer's model dataset and run the program automatically.

1) Run ```incytr-viz-demo``` in your console (no arguments required). If this does not work for you, you can download the file ```incytr_tutorial.zip``` manually from (TBD zenodo repository link), unzip the file, and proceed with the quickstart instructions below

2) After the pathways have loaded, navigate to http://127.0.0.1:8000/ in your web browser. Instructions on using the tool can be found in the "Help" section on the loaded page.


## Quickstart

1) See Input Format section below for details on clusters and pathways files

2) Run incytr: ```incytr-viz --clusters path/to/clusters.csv --pathways path/to/pathways.csv```

3) After pathways have loaded, navigate to http://127.0.0.1:8000/ in your web browser. Instructions on using the tool can be found in the "Help" section on the loaded page. Large datasets (~1M pathways) may take longer on initial load.


## Use

Instructions on using the tool can be found in the "Help" section on the loaded page. You can also find the same information in this repository at src/incytr_viz/assets/help.md


## Input Format

### A. Clusters File

CSV or TSV with the names of experimental conditions and cell populations analyzed by incytr. Column names are case-insensitive

*Important: Ensure that the entries in the "Condition" column in the clusters file match the condition names provided to incytr. Note the Sigprob_5X and SigProb_WT column headers in pathways file (below) align with "5X" and "WT" conditions in the clusters file.*

#### Required columns:

- condition -- Name of experimental condition.
- type -- Name of cell population. Not all population types need to be present in both conditions.

#### Optional columns:
- population: number of cells for that condition and type. This is used only for sizing nodes in the network view

*Note: running ```incytr-viz-demo``` will download an appropriately-formatted example input at incytr_viz_demo/clusters.csv*

#### Example:

```
      +-----------+----------------------+--------------------+
      | Condition |         Type         |     Population     |
      +-----------+----------------------+--------------------+
      |    5X     |  Excitatory.neurons  |         8665       |
      |    5X     | Medium.spiny.neurons |         1590       |
      |    5X     |   Oligodendrocytes   |         2845       |
      |    5X     |  Endothelial.cells   |         311        |
      |    5X     |         OPCs         |         478        |
      |    5X     |      Microglia       |         1486       |
      |    WT     |  Excitatory.neurons  |         11718      |
      |    WT     |     Interneurons     |         2086       |
      |    WT     |      Astrocytes      |         1199       |
      |    WT     |  Endothelial.cells   |          480       |
      |    WT     |         OPCs         |          683       |
      |    WT     |      Microglia       |          803       |
      +-----------+----------------------+--------------------+

```

### B. Pathways File

This file should be a CSV or TSV output previously generated by the incytr analysis package. Some columns are optional and are dependent on the data modalities available

*Important: Ensure that the entries in the "Condition" column in the clusters file match the condition names provided to incytr. Note the Sigprob_5X and SigProb_WT column headers in pathways file (below) align with "5X" and "WT" conditions in the clusters file.*

#### Required columns:
- path -- 4-step network with components separated by asterisks (e.g. A*B*C*D)
- sender -- sending cell population
- receiver -- receiving cell population
- sigprob_X -- signaling probability for experimental (positive aFC) condition ('X' will vary)
- sigprob_Y -- signaling probability for control (negative aFC) condition ('Y' will vary)
- afc -- adjusted fold change

*The left-most sigprob column should be the *

#### Optional columns:
- p_value_X  ('X' will vary)
- p_value_Y  ('Y' will vary)
- TPDS (transcriptomics-based pathway differential score)
- PPDS (proteomics-based pathway differential score)
- sik_r_of_em -- Signaling involved kinase relationship between receptor/effector (receptor is kinase, effector is substrate)
- sik_r_of_t -- Signaling involved kinase relationship between receptor/target gene
- sik_em_of_t -- Signaling involved kinase relationship between effector/target gene
- sik_em_of_r -- Signaling involved kinase relationship between effector/receptor
- sik_t_of_r -- Signaling involved kinase relationship between target gene/receptor
- sik_t_of_em -- Signaling involved kinase relationship between target gene/effector
- umap1 -- optional first 2-d umap coordinate for pathway
- umap2 -- optional second 2-d umap coordinate for pathway

*Note: running ```incytr-viz-demo``` will also download an appropriately-formatted example input at incytr_viz_demo/pathways.csv*

#### Example (with all required and optional columns):

```
+---+------------------------+--------------------+--------------------+-------------+-------------+--------------+      \
|   |          Path          |       Sender       |      Receiver      | SigProb_5X  | SigProb_WT  |     aFC      |      \
+---+------------------------+--------------------+--------------------+-------------+-------------+--------------+      \
| 0 |  Cntn4*App*Aak1*Dock7  | Excitatory.neurons | Excitatory.neurons | 0.819899544 | 0.465168271 | 0.816354872  |      \
| 1 |  Cntn4*App*Aak1*Etl4   | Excitatory.neurons | Excitatory.neurons | 0.886444905 | 0.890006529 | -0.005778447 |      \
| 2 |  Cntn4*App*Aak1*Stk39  | Excitatory.neurons | Excitatory.neurons | 0.293952701 | 0.243354205 | 0.098756719  |      \
| 3 |   Cntn4*App*Aak1*Syp   | Excitatory.neurons | Excitatory.neurons | 0.909707784 | 0.922833896 | -0.020645266 |      \
| 4 |  Cntn4*App*Aak1*Basp1  | Excitatory.neurons | Excitatory.neurons | 0.863494378 | 0.878796052 | -0.025312546 |      \
| 5 |  Cntn4*App*Aak1*Rph3a  | Excitatory.neurons | Interneurons | 0.826835011 | 0.76565119  | 0.110772934  |            \
| 6 |  Cntn4*App*Aak1*Mark2  | Excitatory.neurons | Interneurons | 0.42778563  | 0.359672825 | 0.159537491  |            \
| 7 | Cntn4*App*Aak1*Pitpnm2 | Excitatory.neurons | Interneurons | 0.824016257 | 0.494738382 | 0.734843587  |            \
| 8 | Cntn4*App*Aak1*Slc2a3  | Excitatory.neurons | Interneurons | 0.489517379 | 0.752109522 | -0.618555434 |            \
+---+------------------------+--------------------+--------------------+-------------+-------------+--------------+      \


+---+------------+------------+-------------+------------+-------------+-------------+------------+-------------+ \
|   | p_value_5X | p_value_WT | SiK_R_of_EM | SiK_R_of_T | SiK_EM_of_T | SiK_EM_of_R | SiK_T_of_R | SiK_T_of_EM | \
+---+------------+------------+-------------+------------+-------------+-------------+------------+-------------+ \
| 0 |     0      |     0      |             |            |    Aak1     |             |            |             | \
| 1 |     0.01   |     .1     |             |            |    Aak1     |             |            |             | \
| 2 |     0      |     0      |             |            |    Aak1     |             |            |    Stk39    | \
| 3 |     0      |     0      |             |            |    Aak1     |             |            |             | \
| 4 |     1      |     1      |             |            |    Aak1     |             |            |             | \
| 5 |     0      |     0      |             |            |    Aak1     |             |            |             | \
| 6 |     0      |     0      |             |            |    Aak1     |             |            |    Mark2    | \
| 7 |     0      |     0      |             |            |    Aak1     |             |            |             | \
| 8 |     0      |     0      |             |            |    Aak1     |             |            |             | \
+---+------------+------------+-------------+------------+-------------+-------------+------------+-------------+ \

+---+--------------+--------------+-----------+-----------+
|   |     TPDS     |     PPDS     |   umap1   |   umap2   |
+---+--------------+--------------+-----------+-----------+
| 0 | 0.673081017  | -0.077709086 | 12.823026 | -1.528055 |
| 1 | -0.005778383 | -0.167344206 | 6.683978  | 6.683978  |
| 2 | 0.098436913  | -0.081214482 | -1.552083 | -1.552083 |
| 3 | -0.020642333 | -0.107245137 | -1.578296 | -1.578296 |
| 4 | -0.025307141 | -0.211122838 | -1.457079 | -1.457079 |
| 5 | 0.110322062  | -0.095464928 | -2.996066 | -2.996066 |
| 6 | 0.158197604  | -0.088897533 | 0.496085  | 0.496085  |
| 7 | 0.626019667  | -0.078997138 | -1.435241 | -1.435241 |
| 8 | -0.550121437 | -0.15471897  | 6.417764  | 6.417764  |
+---+--------------+--------------+-----------+-----------+
```




