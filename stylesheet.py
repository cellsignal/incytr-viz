cytoscape_styles = [
    {
        "selector": "node",
        "style": {
            "label": "data(id)",
            "text-wrap": "ellipsis",
            "text-valign": "top",
            "text-halign": "right",
            "font-size": "20px",
            "height": "data(height)",
            "width": "data(width)",
            "background-color": "data(background_color)",
        },
    },
    {
        "selector": "edge",
        "style": {
            "curve-style": "bezier",
            "target-arrow-shape": "vee",
            "arrow-scale": ".75",
            "label": "data(label)",
            "width": "data(width)",
            "line-color": "data(line_color)",
            "target-arrow-color": "data(line_color)",
        },
    },
]
