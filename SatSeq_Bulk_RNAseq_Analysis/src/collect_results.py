import argparse
import os
from glob import glob

import pandas as pd


def collect_boot_results(results_dir: str, output_file: str):
    """
    Collect all _boot.csv files from experiment subdirectories and combine them
    into a multi-sheet Excel file.

    :param results_dir: Path to the results directory containing experiment subdirectories
    :param output_file: Path to the output Excel file
    """
    # Find all boot.csv files
    boot_files = glob(os.path.join(results_dir, '*', '*_boot.csv'))

    if not boot_files:
        print(f"No _boot.csv files found in {results_dir}")
        return

    print(f"Found {len(boot_files)} boot files")

    # Create Excel writer
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for boot_file in sorted(boot_files):
            # Extract experiment name from the directory name
            experiment_name = os.path.basename(os.path.dirname(boot_file))

            # Excel sheet names have a 31 character limit
            sheet_name = experiment_name[:31]

            print(f"  Adding {experiment_name}")

            # Read the CSV and write to Excel
            df = pd.read_csv(boot_file, low_memory=False)
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"Saved combined results to {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect SatSeq boot results into a multi-sheet Excel file.')
    parser.add_argument('results_dir', nargs='?', default='results',
                        help='Path to the results directory (default: results)')
    parser.add_argument('-o', '--output', default='collected_boot_results.xlsx',
                        help='Output Excel file path (default: collected_boot_results.xlsx)')
    args = parser.parse_args()

    collect_boot_results(args.results_dir, args.output)
