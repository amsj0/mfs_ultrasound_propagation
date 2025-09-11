"""
Main script to run mesh-free solution, analysis, and time series.
- Runs mfsolution always.
- Runs analysis if output_path is provided.
- Runs tseries if both output_path and config_file are provided.

author: Agesinaldo Silva
date: June, 2024
"""


from mfsolution import *
from analyse_concurrent import *
from run_tseries_parallel import *
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Run mfsolution, optionally analysis, and optionally tseries."
    )
    p.add_argument("input_file", help="YAML config for mfsolution.")
    # Optional positional: if provided (non-empty), analysis runs
    p.add_argument("output_path", nargs="?", default="",
                   help="Output directory. If provided, analysis will run.")
    # Optional positional: if provided (together with output_path), tseries runs
    p.add_argument("config_file", nargs="?", default=None,
                   help="YAML config for tseries (requires output_path).")
    return p.parse_args()

def main():
    args = parse_args()

    # 1) Always run mfsolution
    print("Running mfsolution code")
    ST = mfsolution(args.input_file, args.output_path, __name__)

    # 2) Run analysis iff output_path was provided (semantics unchanged)
    if args.output_path:
        print("Running analysis code")
        ST = analyse_concurrent(
            store=ST,
            input_file=args.input_file,
            output_path=args.output_path,
            task_name=__name__,
        )

    # 3) Run tseries iff BOTH output_path and config_file were provided
    #    (matches original `if len(sys.argv) == 4:` behavior)
    if args.output_path and args.config_file is not None:
        print("Running tseries code")
        run_tseries_parallel(
            store=ST,
            config_file=args.config_file,
            output_path=args.output_path,
            task_name=__name__,
        )

if __name__ == "__main__":
    main()
