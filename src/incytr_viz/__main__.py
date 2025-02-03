import argparse
import io
import os
import sys
import time
import zipfile

import requests

from incytr_viz.app import create_app
from incytr_viz.util import ascii, create_logger

logger = create_logger(__name__)


import gunicorn.app.base
from gunicorn.arbiter import Arbiter


class CustomArbiter(Arbiter):

    def handle_hup(self):
        self.handle_term()


class StandaloneApplication(gunicorn.app.base.BaseApplication):

    def __init__(self, app, options=None):
        self.options = options or {}
        self.application = app
        super().__init__()

    def load_config(self):
        config = {
            key: value
            for key, value in self.options.items()
            if key in self.cfg.settings and value is not None
        }
        for key, value in config.items():
            self.cfg.set(key.lower(), value)

    def load(self):
        return self.application

    def run(self):
        bind = self.cfg.settings.get("bind").value

        # should always be pointing to a single port on localhost
        if isinstance(bind, list) and len(bind) == 1:
            location = "http://" + bind[0]
        else:
            location = bind

        logger.info(f"Running incytr-viz with gunicorn at {location}")

        try:
            CustomArbiter(self).run()
        except RuntimeError as e:
            print("\nError: %s\n" % e, file=sys.stderr)
            sys.stderr.flush()
            sys.exit(1)


def run_gunicorn(pathways, clusters):

    print(ascii())
    time.sleep(1)
    app = create_app(pathways_file=pathways, clusters_file=clusters)

    g_app = StandaloneApplication(app=app)

    try:
        g_app.cfg.settings["loglevel"].value = "warning"
    except:
        logger.warning("Could not set gunicorn loglevel")

    g_app.run()


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

    run_gunicorn(PATHWAYS, CLUSTERS)


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
    app = create_app(pathways_file=PATHWAYS, clusters_file=CLUSTERS)

    app.run(debug=True)


def demo():
    """
    Downloads a zip file from Zenodo, unzips it, and returns the filepaths
    of the extracted files.

    Args:
        zenodo_url: The URL of the zip file on Zenodo.
        extract_dir: The directory to extract the files to. Defaults to "incytr_viz_demo".

    Returns:
        A list of filepaths of the unzipped files, or None if an error occurs.
    """

    # Example usage:
    zenodo_url = "https://zenodo.org/api/records/14775408/draft/files/incytr_viz_demo.zip/content?download=1&token=eyJhbGciOiJIUzUxMiIsImlhdCI6MTczODI2OTE0MywiZXhwIjoxNzQwNzAwNzk5fQ.eyJpZCI6IjM3MjkwMjhhLTQ4ZDktNDYyOC1hMDc1LThiNjMxZGQ3MDJiNyIsImRhdGEiOnt9LCJyYW5kb20iOiI0NjI1NWFkNWI2NGMyYTc1OTkyMWY0MDc1NjlhNmU2NCJ9.pR2eHubISH-4aC8ozbBbDBAdxKM7EJ9mGSi7lGt6C_oYo5reSe24ADS2hSyw08zwf5eKtq5mYC7w-Bi1zTDK-g"  # Replace with the actual Zenodo URL
    extract_dir = "incytr_viz_demo"

    try:
        # 1. Download the zip file
        print("Downloading the demo zip file from Zenodo...")
        response = requests.get(zenodo_url, stream=True)
        response.raise_for_status()

        os.makedirs(extract_dir, exist_ok=True)

        print(f"Extracting the demo zip file into directory {extract_dir}...")
        time.sleep(0.5)
        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            file_paths = []
            for file_info in z.infolist():  # Iterate through each file in the archive
                # Sanitize file name to avoid path traversal vulnerabilities
                filename = os.path.basename(file_info.filename)  # extract the file name
                if filename:  # if the file name is not empty
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

        run_gunicorn(PATHWAYS, CLUSTERS)

    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
        return None
    except zipfile.BadZipFile as e:
        print(f"Error unzipping file: {e}")
        return None
    except Exception as e:  # Catch other potential errors (e.g., file I/O)
        print(f"An unexpected error occurred: {e}")
        return None


if __name__ == "__main__":
    main()
