"""Cell elongation analysis of combined Workflow1 tables, optionally comparing Cellpose models.

Give one combined table per model (the All_Cells_Shape_slim.csv that Workflow1 writes into each
Results_<model> folder). The model name is taken from --names, or
from the ResultsFolder column.

Outputs (in --out):
  elongation_image_level.csv     per image and model: n cells, median aspect ratio, % cells AR >= 3
  elongation_group_summary.csv   mean +- SD over images per time point x condition x model
  elongation_percentiles.csv     per-cell percentiles of the aspect ratio
  model_agreement.csv            cell matching between the first model and each other model
  elongation_distributions.png   per-cell histograms (first model)
  elongation_distributions.csv   the histograms bin by bin (% of cells per 0.25 bin), for every model;
                                 *_wide.csv has one column per group, e.g. for Excel or Prism

Nuclei-only filter (optional): Cellpose sometimes segments only the nucleus of a cell. Such segments
are removed when their area is at or below a size threshold, given directly (--min-area, um2) or taken
as a percentile of the nuclear areas measured by Workflow2 times a factor (--nuclei-table,
--nuclei-percentile, --nuclei-factor; default 1.2 x the 95th percentile, to also catch a nucleus
plus a little cytoplasm).
  filter_summary.csv             % of segments removed per image and model

Usage:
  python elongation_analysis.py cpsam.csv cpsam_v2.csv cpdino.csv --names cpsam cpsam_v2 cpdino --out results_shape
  python elongation_analysis.py cpsam.csv --names cpsam --nuclei-table All_Cells_Hes_slim.csv --out results_filtered
  python elongation_analysis.py cpsam.csv --names cpsam --min-area 362 --exclude day5_notch_inh_CDH5AF647_Image_3 --out r
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from common import add_design_columns, mean_sd, histogram_table, save_histogram_tables, TIME_ORDER, COND_ORDER

AR = "Ellipse.Elong"
COLORS = {"control": "#9a9a9a", "act": "#2b6cb0", "inh": "#e07b39"}
NAMES = {"control": "Control", "act": "Activated", "inh": "Inhibited"}


def match_cells(a, b, max_dist_frac=0.3, max_area_change=0.3):
    """Match cells of model a to model b: mutual nearest centroid, distance < max_dist_frac x
    equivalent radius of the cell in a, and area ratio within +-max_area_change."""
    out = []
    for sf, ga in a.groupby("SourceFile"):
        gb = b[b.SourceFile == sf]
        if len(gb) == 0:
            continue
        pa, pb = ga[["Centroid.X", "Centroid.Y"]].values, gb[["Centroid.X", "Centroid.Y"]].values
        da, ia = cKDTree(pb).query(pa)
        _, ib = cKDTree(pa).query(pb)
        radius = np.sqrt(ga.Area.values / np.pi)
        area_ratio = gb.Area.values[ia] / ga.Area.values
        ok = (ib[ia] == np.arange(len(ga))) & (da < max_dist_frac * radius) & \
             (area_ratio > 1 - max_area_change) & (area_ratio < 1 / (1 - max_area_change))
        out.append(pd.DataFrame({"time": ga.time.values, "cond": ga.cond.values,
                                 "ar_a": ga[AR].values, "ar_b": gb[AR].values[ia], "matched": ok}))
    return pd.concat(out, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tables", nargs="+")
    ap.add_argument("--names", nargs="*")
    ap.add_argument("--out", default="results_shape")
    ap.add_argument("--min-area", type=float, help="remove segments with Area <= this value (um2)")
    ap.add_argument("--nuclei-table", help="Workflow2 table; the threshold is a percentile of its nuclear Area")
    ap.add_argument("--nuclei-percentile", type=float, default=95.0)
    ap.add_argument("--nuclei-factor", type=float, default=1.2,
                    help="multiplies the nuclear-area percentile, to also catch a nucleus plus a little cytoplasm (default 1.2)")
    ap.add_argument("--exclude", nargs="*", default=[], help="image names (SourceFile) to leave out entirely")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    frames = []
    for i, path in enumerate(args.tables):
        t = add_design_columns(pd.read_csv(path))
        t["model"] = args.names[i] if args.names else t["ResultsFolder"].astype(str)
        frames.append(t)
    d = pd.concat(frames, ignore_index=True)
    if args.exclude:
        before = d.SourceFile.nunique()
        d = d[~d.SourceFile.isin(args.exclude)].reset_index(drop=True)
        print(f"Excluded {before - d.SourceFile.nunique()} image(s): {', '.join(args.exclude)}")

    # nuclei-only filter
    thr = args.min_area
    if args.nuclei_table:
        nuc = pd.read_csv(args.nuclei_table)
        pct_area = float(np.nanpercentile(nuc["Area"], args.nuclei_percentile))
        thr = pct_area * args.nuclei_factor
        print(f"Nuclei-only threshold: {args.nuclei_factor:g} x {args.nuclei_percentile:g}th percentile of "
              f"{len(nuc):,} nuclear areas ({pct_area:.0f}) = {thr:.0f}")
    if thr is not None:
        d["nuclei_only"] = d["Area"] <= thr
        fs = d.groupby(["model", "SourceFile"], observed=True).agg(
            n_segments=("Area", "size"), n_removed=("nuclei_only", "sum")).reset_index()
        fs["pct_removed"] = 100 * fs.n_removed / fs.n_segments
        fs.to_csv(os.path.join(args.out, "filter_summary.csv"), index=False)
        print(f"Removed {d.nuclei_only.sum():,} of {len(d):,} segments with Area <= {thr:.0f} (nuclei-only)")
        d = d[~d.nuclei_only].reset_index(drop=True)

    models = list(dict.fromkeys(d.model))
    d["model"] = pd.Categorical(d.model, models)
    print(d.groupby("model", observed=True).agg(cells=(AR, "size"), images=("SourceFile", "nunique")))

    keys = ["model", "time", "cond", "img"]
    img = d.groupby(keys, observed=True).agg(
        SourceFile=("SourceFile", "first"), n_cells=(AR, "size"), median_ar=(AR, "median"), mean_ar=(AR, "mean"),
        pct_ar_ge3=(AR, lambda s: 100 * (s >= 3).mean()), median_area=("Area", "median")).reset_index()
    img.to_csv(os.path.join(args.out, "elongation_image_level.csv"), index=False)

    rows = []
    for col in ["median_ar", "pct_ar_ge3", "n_cells", "median_area"]:
        s = mean_sd(img.set_index(["model", "time", "cond"])[col], ["time", "cond", "model"]).reset_index()
        s.insert(0, "measure", col)
        rows.append(s)
    summary = pd.concat(rows)
    summary.to_csv(os.path.join(args.out, "elongation_group_summary.csv"), index=False)
    print(summary[summary.measure == "median_ar"].pivot_table(index=["time", "cond"], columns="model", values="mean", observed=True).round(2))

    pct = d.groupby(["time", "cond", "model"], observed=True)[AR].quantile([.1, .25, .5, .75, .9]).unstack()
    pct.to_csv(os.path.join(args.out, "elongation_percentiles.csv"))

    # model agreement (first model vs each other model)
    agree = []
    ref = d[d.model == models[0]]
    for m in models[1:]:
        mt = match_cells(ref, d[d.model == m])
        mm = mt[mt.matched]
        bins = pd.cut(mt.ar_a, [1, 1.5, 2, 3, 4, 6, np.inf])
        for b, grp in mt.groupby(bins, observed=False):
            agree.append({"reference": models[0], "other": m, "ar_bin": str(b), "n_cells": len(grp),
                          "pct_unmatched": 100 * (1 - grp.matched.mean()) if len(grp) else np.nan})
        print(f"{models[0]} vs {m}: {100 * mt.matched.mean():.0f}% matched; Spearman AR of matched cells "
              f"{mm.ar_a.corr(mm.ar_b, method='spearman'):.3f}; median |dAR| {np.median(np.abs(mm.ar_a - mm.ar_b)):.3f}")
    if agree:
        pd.DataFrame(agree).to_csv(os.path.join(args.out, "model_agreement.csv"), index=False)
    else:
        print("One model only - model comparison skipped.")

    # distributions for the first model
    bins = np.arange(1, 6.26, 0.25)
    times = [t for t in TIME_ORDER if (ref.time == t).any()]
    fig, axs = plt.subplots(1, len(times), figsize=(3.4 * len(times), 3), sharey=True, squeeze=False)
    axs = axs[0]
    for ax, t in zip(axs, times):
        for c in COND_ORDER:
            x = ref.loc[(ref.time == t) & (ref.cond == c), AR].dropna()
            if len(x) == 0:
                continue
            h = histogram_table(x, bins, clip_low=False)   # aspect ratio >= 1: no open lower bin
            ax.plot(h.bin_center, h.pct, color=COLORS[c], label=f"{NAMES[c]} ({len(x):,})")
        ax.set_title(f"{t} ({models[0]})", loc="left"); ax.set_xlabel("Aspect ratio (last bin: 6+)")
        ax.legend(frameon=False, fontsize=7)
    axs[0].set_ylabel("% of cells")
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "elongation_distributions.png"), dpi=200); plt.close(fig)

    # the same histograms bin by bin, for every model (the PNG shows the first model)
    hist_rows = []
    for m in models:
        for t in times:
            for c in COND_ORDER:
                x = d.loc[(d.model == m) & (d.time == t) & (d.cond == c), AR].dropna()
                if len(x) == 0:
                    continue
                h = histogram_table(x, bins, clip_low=False)
                hist_rows.append(h.assign(model=m, time=t, cond=c, n=len(x)))
    if hist_rows:
        long_df = pd.concat(hist_rows, ignore_index=True)[
            ["model", "time", "cond", "n", "bin_start", "bin_end", "bin_center", "count", "pct"]]
        save_histogram_tables(long_df, ["model", "time", "cond"], os.path.join(args.out, "elongation_distributions"))
    print(f"Results written to {args.out}/")


if __name__ == "__main__":
    main()
