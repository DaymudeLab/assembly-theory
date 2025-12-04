"""
bench_stats: Compute benchmark statistics from assembly-theory (Rust, us),
             assembly_go (Go), and assemblycpp-v5 (C++) outputs.
"""

import argparse
from collections import defaultdict
import json
import numpy as np
import os.path as osp
import pandas as pd
import scipy.stats as sps


if __name__ == "__main__":
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    def_at_path = osp.join("..", "target", "criterion", "joss")
    parser.add_argument(
        "-A",
        "--at_path",
        type=str,
        default=def_at_path,
        help="Path to assembly-theory benchmark output",
    )
    def_ago_path = osp.join("assembly_go", "cmd", "app", "joss_bench.tsv")
    parser.add_argument(
        "-G",
        "--ago_path",
        type=str,
        default=def_ago_path,
        help="Path to assembly_go benchmark output",
    )
    args = parser.parse_args()

    # Define reference datasets.
    datasets = ["gdb13_1201", "gdb17_200", "checks", "coconut_55"]

    # Create results dict of the form {bounds: {dataset: (mean, 95% conf.)}}.
    results = defaultdict(dict)

    # Load assembly_go benchmark outputs.
    ago_df = pd.read_csv(
        args.ago_path,
        sep="\t",
        header=None,
        names=["bench", "reps", "time"],
        engine="python",
        skiprows=4,
        skipfooter=2,
    ).dropna()
    ago_df["dataset"] = ago_df.bench.apply(
        lambda x: x.strip().split("/")[1].split("-")[0]
    )
    ago_df.time = ago_df.time.apply(lambda x: float(x.split()[0]))

    # Compute assembly_go means and 95% confidence intervals.
    for dataset in datasets:
        # Skip very slow benchmarks that we don't do multiple samples for.
        if dataset == "coconut_55":
            continue

        times = np.array(ago_df.loc[ago_df["dataset"] == dataset].time)
        mean_time = times.mean()
        conf = sps.t.interval(0.95, len(times) - 1, loc=mean_time, scale=sps.sem(times))
        conf_perc = (conf[1] - conf[0]) / 2 / mean_time * 100
        results["assembly_go"][dataset] = (mean_time, conf_perc)

    # Do the same for assembly-theory's different bounds.
    for dataset in datasets:
        with open(osp.join(args.at_path, dataset, "new", "estimates.json"), "r") as f:
            stats = json.load(f)
            mean_time = stats["mean"]["point_estimate"]
            conf_low = stats["mean"]["confidence_interval"]["lower_bound"]
            conf_high = stats["mean"]["confidence_interval"]["upper_bound"]
            conf_perc = (conf_high - conf_low) / 2 / mean_time * 100
            results["assembly-theory"][dataset] = (mean_time, conf_perc)

    # Print results.
    results_df = pd.DataFrame(results)
    results_df = results_df.map(
        lambda x: f"{x[0] / 1e9:.3f} s ± {x[1]:.2f}%" if not pd.isna(x) else "NaN"
    )
    print(results_df)
