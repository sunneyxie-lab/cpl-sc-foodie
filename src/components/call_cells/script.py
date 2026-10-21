import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from kneed import KneeLocator

## VIASH START
par = {
    "input": "grouped.bam",
    "id": "sample",
    "output_cell_list": "cell_list.tsv",
    "output_cell_counts": "cell_counts.tsv",
    "output_cell_qc": "cell_qc.tsv",
    "output_knee_plot": "knee_plot.png",
    "drop_lowest": 1000,
    "max_cells": 2000,
    "fallback_cells": 1000,
}
meta = {}
## VIASH END


def count_groups(path):
    groups = defaultdict(set)
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if not read.has_tag("CB"):
                sys.exit(f"call_cells: read {read.query_name} has no CB tag")
            mis = groups[read.get_tag("CB")]
            if read.has_tag("MI"):
                mis.add(read.get_tag("MI"))
    df = pd.DataFrame(
        {"barcode": list(groups), "count": [len(s) for s in groups.values()]},
        columns=["barcode", "count"],
    ).astype({"barcode": str, "count": "int64"})
    df = df.sort_values(["count", "barcode"], kind="stable").reset_index(drop=True)
    return pd.Series(df["count"].values, index=df["barcode"].values, dtype="int64")


def find_knee(remaining):
    n = len(remaining)
    if n == 0:
        return None
    try:
        elbow = KneeLocator(
            np.arange(n), remaining.values, S=1.0, curve="convex", direction="increasing", online=True
        ).elbow
    except Exception:
        return None
    if elbow is None or not 0 <= elbow < n:
        return None
    return int(elbow)


def write_tsv(series, path):
    series.to_csv(path, sep="\t", header=False)


def hist(ax, counts, log_scale):
    if log_scale:
        counts = counts[counts > 0]
    if len(counts) == 0:
        return False
    sns.histplot(counts.values, kde=counts.nunique() > 1, log_scale=log_scale, fill=False, bins=100, ax=ax)
    return True


def plot(path, remaining, kept, not_kept, knee_idx, knee_count, title):
    fig = plt.figure(figsize=(14, 8))

    ax = fig.add_subplot(2, 2, 1)
    if len(remaining):
        ax.scatter(np.arange(len(remaining)), remaining.values, alpha=0.3, s=5)
        if knee_idx is not None:
            ax.axvline(knee_idx, color="red", linestyle="--", label=f"knee index: {knee_idx}")
            ax.axhline(knee_count, color="red", linestyle="--", label=f"count at knee: {knee_count}")
            ax.legend()
        ax.set_xlabel("Barcode rank (ascending, after drop)")
        ax.set_ylabel("Group ids per barcode")

    ax = fig.add_subplot(2, 2, 2)
    groups = pd.DataFrame(
        {
            "count": np.concatenate([kept.values, not_kept.values]),
            "group": ["kept"] * len(kept) + ["not_kept"] * len(not_kept),
        }
    )
    groups = groups[groups["count"] > 0]
    if len(groups):
        order = [g for g in ["kept", "not_kept"] if (groups["group"] == g).any()]
        sns.boxplot(groups, x="group", y="count", order=order, showfliers=False, color="black", fill=False, gap=0.65, ax=ax)
        sns.stripplot(groups, x="group", y="count", order=order, size=4, alpha=0.2, jitter=True, ax=ax)
        for i, g in enumerate(order):
            median = groups.loc[groups["group"] == g, "count"].median()
            ax.text(i, median, f"{median:g}", ha="center", va="bottom")
        if knee_count is not None:
            ax.axhline(knee_count, color="red", linestyle="--")
        ax.set_yscale("log")
        ax.set_ylabel("Group ids per barcode (log)")

    ax = fig.add_subplot(2, 2, 3)
    if hist(ax, kept, log_scale=True):
        if knee_count is not None and knee_count > 0:
            ax.axvline(knee_count, color="red", linestyle="--")
        ax.set_xlabel("Group ids per kept barcode")

    ax = fig.add_subplot(2, 2, 4)
    if not_kept.nunique() > 1:
        hist(ax, not_kept, log_scale=False)
        ax.set_xlabel("Group ids per non-kept barcode")

    fig.suptitle(title, fontsize=16)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


counts = count_groups(par["input"])
write_tsv(counts, par["output_cell_counts"])

remaining = counts.iloc[par["drop_lowest"]:]
knee_idx = find_knee(remaining)
knee_count = None
if knee_idx is not None:
    knee_count = int(remaining.iloc[knee_idx])
    knee_cells = remaining.iloc[knee_idx:]
    if len(knee_cells) < par["max_cells"]:
        kept, state = knee_cells, "Normal"
    else:
        kept, state = counts.iloc[-par["fallback_cells"]:], "Abnormal"
else:
    kept, state = counts.iloc[-par["fallback_cells"]:], "Abnormal"
write_tsv(kept, par["output_cell_list"])

nonzero = counts[counts > 0]
not_kept = nonzero[~nonzero.index.isin(kept.index)]
kept_count, total_count = int(kept.sum()), int(counts.sum())
kept_share = f"{kept_count / total_count:.2%}" if total_count else "NA"
fmt = lambda v: "NA" if v is None else str(v)

qc = pd.DataFrame(
    [[par["id"], fmt(knee_idx), fmt(knee_count), len(kept), len(nonzero), kept_count, total_count, kept_share, state]],
    columns=["id", "knee_index", "knee_count", "kept_cells", "nonzero_barcodes", "kept_count", "total_count", "kept_share", "state"],
)
qc.to_csv(par["output_cell_qc"], sep="\t", index=False)

title = "\n".join(
    [
        par["id"],
        f"knee index: {fmt(knee_idx)}; count at knee: {fmt(knee_count)}; kept cells: {len(kept)}; non-zero barcodes: {len(nonzero)}",
        f"kept/all group ids: {kept_count}/{total_count} = {kept_share}; {state}",
    ]
)
plot(par["output_knee_plot"], remaining, kept, not_kept, knee_idx, knee_count, title)
