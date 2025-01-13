from incytr_viz.app import create_app

group_a_populations = "/home/icossentino/code/incytr-viz/data/5x/5X_pop.csv"
group_b_populations = "/home/icossentino/code/incytr-viz/data/5x/WT_pop.csv"
pathways = (
    "/home/icossentino/code/incytr-viz/data/5x/umap_embedding_60_1_5x_all_columns.csv"
)

create_app(pathways, group_a_populations, group_b_populations).run(debug=True)
