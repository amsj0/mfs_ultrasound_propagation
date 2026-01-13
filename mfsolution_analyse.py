import sys

from mfsolution import *
from analyse_concurrent import *
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Run mfsolution and optionally analysis."
    )
    p.add_argument("input_file", help="YAML config for mfsolution.")
    # Optional positional: if provided (non-empty), analysis runs
    p.add_argument("output_path", nargs="?", default="",
                   help="Output directory. If provided, analysis will run.")
    # Optional positional: if provided (together with output_path), tseries runs
    return p.parse_args()

def main():
    args = parse_args()

    # 1) Always run mfsolution
    print("Running mfsolution code")
    ST = mfsolution(args.input_file, args.output_path, "mfsolution_analyse")

    # 2) Run analysis iff output_path was provided (semantics unchanged)
    if args.output_path:
        print("Running analysis code")
        analyse_concurrent(
            store=ST,
            input_file=args.input_file,
            output_path=args.output_path,
            task_name="mfsolution_analyse",
        )


if __name__ == "__main__":
    main()