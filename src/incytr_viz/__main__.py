import argparse
import sys
from incytr_viz.app import create_app
from incytr_viz.util import create_logger

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
        logger.info(
            f"Running incytr-viz with gunicorn at {self.cfg.settings['bind'].value}"
        )
        try:
            CustomArbiter(self).run()
        except RuntimeError as e:
            print("\nError: %s\n" % e, file=sys.stderr)
            sys.stderr.flush()
            sys.exit(1)


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
    app = create_app(pathways_file=PATHWAYS, clusters_file=CLUSTERS)
    StandaloneApplication(app=app).run()


if __name__ == "__main__":
    main()
