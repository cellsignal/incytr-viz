import argparse
import io
import os
import sys
import time
import zipfile

import requests

from incytr_viz.app import create_dash_app
from incytr_viz.util import create_logger

logger = create_logger(__name__)


def run_wsgi(pathways, clusters):
    if os.name == "nt":
        from incytr_viz.wsgi_windows import run_waitress

        run_waitress(pathways, clusters)
    else:
        from incytr_viz.wsgi_posix import run_gunicorn

        run_gunicorn(pathways, clusters)


def main():
    parser = argparse.ArgumentParser(description="Run the InCytr visualization app.")
    parser.add_argument(
        "--clusters",
        type=str,
        required=True,
        help="cell clusters filepath",
    )
    parser.add_argument("--pathways", type=str, required=True, help="pathways filepath")

    args = parser.parse_args()

    PATHWAYS = args.pathways
    CLUSTERS = args.clusters

    run_wsgi(PATHWAYS, CLUSTERS)


def develop():
    parser = argparse.ArgumentParser(description="Run the InCytr visualization app.")
    parser.add_argument(
        "--clusters",
        type=str,
        required=True,
        help="cell clusters filepath",
    )
    parser.add_argument("--pathways", type=str, required=True, help="pathways filepath")

    args = parser.parse_args()

    PATHWAYS = args.pathways
    CLUSTERS = args.clusters

    logger.info("Running Incytr Viz using gunicorn web server")
    app = create_dash_app(pathways_file=PATHWAYS, clusters_file=CLUSTERS)

    app.run(debug=True)


def demo():
    """
    Downloads a zip file from Zenodo, unzips it, and returns the filepaths
    of the extracted files.

    Args:
        zenodo_url: The URL of the zip file on Zenodo.
        extract_dir: The directory to extract the files to. Defaults to "incytr_viz_demo".

    Returns: None
    """
    zenodo_url = (
        "https://zenodo.org/records/14775408/files/incytr_viz_tutorial.zip?download=1"
    )

    extract_dir = "incytr_viz_demo"

    try:
        print("Downloading the demo zip file from Zenodo...")
        response = requests.get(zenodo_url, stream=True)
        response.raise_for_status()

        os.makedirs(extract_dir, exist_ok=True)

        print(f"Extracting the demo zip file into directory {extract_dir}...")
        time.sleep(0.5)
        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            file_paths = []
            for file_info in z.infolist():
                filename = os.path.basename(file_info.filename)
                if filename:
                    filepath = os.path.join(extract_dir, filename)
                    with open(filepath, "wb") as f:
                        f.write(z.read(file_info))
                    file_paths.append(filepath)

        PATHWAYS = next(
            f for f in file_paths if os.path.basename(f).startswith("pathways")
        )
        CLUSTERS = next(
            f for f in file_paths if os.path.basename(f).startswith("clusters")
        )

        run_wsgi(PATHWAYS, CLUSTERS)

    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
        return None
    except zipfile.BadZipFile as e:
        print(f"Error unzipping file: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None


if __name__ == "__main__":
    main()
