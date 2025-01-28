## Install

*It is recommended to install this program in a virtual environment, using tools such as venv, virtualenv, or conda

```conda create -n incytr-viz python=3.11```

```pip install incytr-viz```





## Input Format

Column names are case-insensitive.

#### A. Clusters File

      TSV/CSV file with the names of experimental conditions and cell populations analyzed by incytr.

      Required columns:

      - condition -- Name of experimental condition.
      - type -- Name of cell population. Not all population types need to be present in both conditions.

      #### Example

      ```
            Condition,Type,Population
            5X,Interneurons,256
            5X,Astrocytes,50
            5X,Microglia,300
            WT,Excitatory.neurons,399
            WT,Interneurons,234
            WT,Astrocytes,124
            WT,OPCs,235
            WT,Microglia,3500

      ```

#### B. Pathway File

This file should be a CSV or TSV that is the main output of the incytr analysis program.

*Notes on Input Format*
-- Incytr has a flexible input format 


#### Required:
- path
- sender
- receiver
- sigprob_\<control group name>\
- sigprob_ \<experimental group name>\
- afc

#### Optional:
- p_value_\<control group name>\
- p_value_ \<experimental group name>\
- tprs
- prs
- sik_r_of_em
- sik_r_of_t
- sik_em_of_t
- sik_em_of_r
- sik_t_of_r
- sik_t_of_em
- umap1
- umap2



## Run

```incytr-viz --clusters <path/to/clusters.csv> --pathways <path/to/pathways.csv>```

Navigate to http://127.0.0.1:8000/ in your browser after the program has loaded. Large datasets (~1M pathways) may take longer on initial load.



