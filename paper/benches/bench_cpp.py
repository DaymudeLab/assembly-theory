"""
bench_cpp: Benchmark assemblycpp-v5 on the assembly-theory reference datasets
"""

import argparse
from collections import defaultdict
import csv
import os
import os.path as osp
import subprocess as sp


if __name__ == "__main__":
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-E",
        "--exe_path",
        type=str,
        default=osp.join("build", "bin", "assembly"),
        help="Path to assemblycpp-v5 executable",
    )
    parser.add_argument(
        "-W",
        "--warmup",
        type=int,
        default=5,
        help="Number of warmup samples to run but not benchmark",
    )
    parser.add_argument(
        "-S",
        "--samples",
        type=int,
        default=20,
        help="Number of benchmark samples per reference dataset",
    )
    args = parser.parse_args()

    # Define reference datasets.
    datasets = ["gdb13_1201", "gdb17_200", "checks", "coconut_55"]

    # Run the benchmark on each reference dataset.
    sample_times = defaultdict(list)
    for dataset in datasets:
        # Get the paths for all .mol files, minus their extensions.
        print(f"{dataset}:")
        path = osp.join("..", "..", "data", dataset)
        mol_files = [
            osp.splitext(osp.join(path, f))[0]
            for f in os.listdir(osp.join(path))
            if osp.isfile(osp.join(path, f)) and osp.splitext(f)[1] == ".mol"
        ]

        # Now that things are warm, do the actual benchmark.
        for s in range(args.warmup + args.samples):
            if s < args.warmup:
                print(f"\tWarming up on {s + 1} of {args.warmup} samples...", end="\r")
            else:
                print(
                    f"\tBenchmarking {s - args.warmup + 1} of {args.samples} samples...",
                    end="\r",
                )

            total_time = 0
            for mol_file in mol_files:
                # Our modified version of assemblycpp-v5 prints only the
                # number of nanoseconds required to run the improvedBnB search
                # function. So we can capture and parse this output to get the
                # time per molecule and add it to the dataset time.
                result = sp.run([args.exe_path, mol_file], capture_output=True)
                total_time += int(result.stdout.strip())

            # Record the benchmark time if we are out of the warmup phase.
            if s >= args.warmup:
                sample_times[dataset].append(total_time)

        # Print benchmark result for this dataset.
        print(f"\n\tMean sample: {sum(sample_times[dataset]) / args.samples} ns")

    # Write the sample benchmark times to a .csv file.
    with open("joss_bench.csv", "w") as f:
        fw = csv.writer(f)
        fw.writerow(datasets)
        for s in range(args.samples):
            fw.writerow([sample_times[d][s] for d in datasets])
