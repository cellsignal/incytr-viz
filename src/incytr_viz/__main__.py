import argparse
import subprocess
from incytr_viz.util import create_logger

logger = create_logger(__name__)


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

    logger.info("Running Incytr Viz using gunicorn web server")
    logger.info("Serving at http://localhost:8000")
    subprocess.run(
        [
            "gunicorn",
            f"incytr_viz.app:create_app(pathways_file='{PATHWAYS}', clusters_file='{CLUSTERS}')",
        ]
    )


if __name__ == "__main__":
    main()
