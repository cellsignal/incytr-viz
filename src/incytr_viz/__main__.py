import argparse
import subprocess


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

    app_string = f"incytr_viz.app:get_server(raw_pathways='{PATHWAYS}', raw_clusters='{CLUSTERS}')"

    subprocess.run(["gunicorn", app_string])


if __name__ == "__main__":
    main()
