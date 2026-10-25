"""Shared helpers: parse condition / time point / image number from Workflow file names.

File names are expected to look like
    <time>_notch_<condition>_<anything>[_Image_<n>]
"_Image_<n>" is optional (TIF inputs can be named freely); if any name lacks it, all images
are numbered 1, 2, 3, ... within each time point x condition, in name order.
e.g. 12h_notch_act_cd31AF555_Hes1AF647_Image_1. Adjust the patterns below if yours differ.
"""
import re
import pandas as pd

TIME_ORDER = ["12h", "48h", "day5"]
COND_ORDER = ["control", "act", "inh"]
TIME_PATTERN = r"^(12h|48h|day5)"          # case-insensitive
COND_PATTERN = r"notch_([A-Za-z]+)_"       # condition is the word after "notch_"


def add_design_columns(df, source_col="SourceFile"):
    """Add categorical 'time', 'cond' and integer 'img' columns parsed from the file name."""
    src = df[source_col].astype(str)
    time = src.str.extract(TIME_PATTERN, flags=re.IGNORECASE)[0].str.lower()
    cond = src.str.extract(COND_PATTERN)[0].str.lower()
    df = df.copy()
    df["time"] = pd.Categorical(time, TIME_ORDER, ordered=True)
    df["cond"] = pd.Categorical(cond, COND_ORDER, ordered=True)
    img = pd.to_numeric(src.str.extract(r"Image_(\d+)", flags=re.IGNORECASE)[0], errors="coerce")
    if img.isna().any():
        # some names have no "Image_<n>" (e.g. TIFs named freely): number ALL images
        # 1, 2, 3, ... within each time point x condition, in name order (mixing
        # explicit and generated numbers could give two images the same number)
        order = (df.assign(_src=src)[["time", "cond", "_src"]].drop_duplicates()
                 .sort_values("_src"))
        order["_n"] = order.groupby(["time", "cond"], observed=True).cumcount() + 1
        fallback = df.assign(_src=src).merge(order, on=["time", "cond", "_src"], how="left")["_n"].values
        img = pd.Series(fallback, index=img.index)
    df["img"] = img.astype(int)
    unparsed = df["time"].isna() | df["cond"].isna()
    if unparsed.any():
        bad = df.loc[unparsed, source_col].unique()[:5]
        raise ValueError(f"Could not parse time/condition from: {list(bad)} - edit the patterns in common.py")
    return df


def mean_sd(series_by_image, groups):
    """Mean and SD over images, per group."""
    g = series_by_image.groupby(groups, observed=True)
    return pd.DataFrame({"mean": g.mean(), "sd": g.std(), "n_images": g.size()})


def histogram_table(values, bins, clip_low=True, clip_high=True):
    """Bin-by-bin histogram of 'values' exactly as plotted in the distribution PNGs.

    Values outside the bin range are counted in the first / last bin (when clip_low /
    clip_high), so those bins are open-ended: their bin_start / bin_end is -inf / inf.
    Returns a DataFrame with bin_start, bin_end, bin_center, count and pct (% of values).
    """
    import numpy as np
    x = pd.Series(values).dropna()
    lo = bins[0] if clip_low else None
    hi = bins[-1] - 1e-6 if clip_high else None
    x = x.clip(lower=lo, upper=hi)
    counts, _ = np.histogram(x, bins=bins)
    starts = np.array(bins[:-1], dtype=float)
    ends = np.array(bins[1:], dtype=float)
    centers = starts + np.diff(bins) / 2
    if clip_low:
        starts[0] = -np.inf
    if clip_high:
        ends[-1] = np.inf
    n = len(x)
    return pd.DataFrame({"bin_start": starts, "bin_end": ends, "bin_center": centers,
                         "count": counts, "pct": 100 * counts / n if n else np.nan})


def save_histogram_tables(long_df, group_cols, out_base):
    """Write the histogram rows as <out_base>.csv (long: one row per group and bin) and
    <out_base>_wide.csv (one row per bin, one % column per group, e.g. for Excel/Prism)."""
    long_df.to_csv(out_base + ".csv", index=False)
    lab = long_df[group_cols].astype(str).agg(" ".join, axis=1) + " (n=" + long_df["n"].astype(str) + ")"
    wide = (long_df.assign(group=lab)
            .pivot_table(index=["bin_start", "bin_end", "bin_center"], columns="group",
                         values="pct", sort=False)
            .reset_index())
    wide.columns.name = None
    wide.to_csv(out_base + "_wide.csv", index=False)
