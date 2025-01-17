import argparse
import subprocess
import os
from incytr_viz.util import create_logger

logger = create_logger(__name__)


def jupyter(pathways, clusters):

    os.environ["INCYTR_PATHWAYS"] = pathways
    os.environ["INCYTR_CLUSTERS"] = clusters

    # load after environment variables are set
    from incytr_viz.app import app

    app.run(debug=True)


def main():
    parser = argparse.ArgumentParser(description="Run the InCytr visualization app.")
    parser.add_argument(
        "--clusters",
        type=str,
        required=True,
        help="Path to clusters A CSV file",
    )
    parser.add_argument(
        "--pathways", type=str, required=True, help="Path to pathways CSV file"
    )

    args = parser.parse_args()

    PATHWAYS = args.pathways
    CLUSTERS = args.clusters

    os.environ["INCYTR_PATHWAYS"] = PATHWAYS
    os.environ["INCYTR_CLUSTERS"] = CLUSTERS

    # load after environment variables are set
    from incytr_viz.app import app

    app.run(debug=True)

    # subprocess.run(["gunicorn", "incytr_viz.app:server"])


if __name__ == "__main__":
    main()
