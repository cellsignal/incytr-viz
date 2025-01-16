from incytr_viz.app import create_app

clusters = "/home/icossentino/code/incytr-viz/data/5x/pop.csv"
pathways = (
    "/home/icossentino/code/incytr-viz/data/5x/umap_embedding_60_1_5x_all_columns.csv"
)
# clusters = "/home/icossentino/code/incytr-viz/data/mc38_011325/pop.csv"
# pathways = "/home/icossentino/code/incytr-viz/data/mc38_011325/mc38_hegs_degs_proteomics_kldb_jan2025_CH.csv"


create_app(pathways, clusters).run(debug=True)
