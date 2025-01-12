import argparse
import gunicorn
import subprocess
import gunicorn.app
from incytr_viz.app import create_app


def main():
    parser = argparse.ArgumentParser(description="Run the InCytr visualization app.")
    parser.add_argument(
        "--group_a_populations",
        type=str,
        required=True,
        help="Path to clusters A CSV file",
    )
    parser.add_argument(
        "--group_b_populations",
        type=str,
        required=True,
        help="Path to clusters B CSV file",
    )
    parser.add_argument(
        "--pathways", type=str, required=True, help="Path to pathways CSV file"
    )

    args = parser.parse_args()

    PATHWAYS = args.pathways
    CLUSTERS_A = args.group_a_populations
    CLUSTERS_B = args.group_b_populations

    app_string = f"incytr_viz.app:create_app(pathways='{PATHWAYS}', clusters_a='{CLUSTERS_A}', clusters_b='{CLUSTERS_B}')"

    subprocess.run(["gunicorn", app_string])


if __name__ == "__main__":
    main()
