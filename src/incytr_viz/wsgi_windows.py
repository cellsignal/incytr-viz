import waitress

from incytr_viz.app import create_app
from incytr_viz.util import create_logger

logger = create_logger(__name__)


def run_waitress(pathways, clusters):

    port = 8000
    app = create_app(pathways_file=pathways, clusters_file=clusters)
    logger.info(f"Running with waitress wsgi at http://127.0.0.1:{port}")
    waitress.serve(app, port=port)
