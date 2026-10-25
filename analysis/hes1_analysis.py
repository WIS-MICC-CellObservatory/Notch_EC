"""Per-nucleus Hes1 analysis of the combined Workflow2 table (All_Cells_Hes_slim.csv).

Main measure (as in the paper's methods): background-corrected nuclear Hes1 = nuclear mean minus
the median of a ring 3-8 px around the nucleus (column Mean_minus_BgRing from Workflow2),
averaged per image, then mean +- SD over images per time point x condition.
The global-background correction and the Hes1-positive analysis are also computed, as extras.

Outputs (in --out):
  hes1_image_level.csv        per image: n nuclei, mean ring-corrected / global-corrected / raw Hes1, % positive
  hes1_group_summary.csv      mean +- SD over images per time point x condition
  hes1_threshold_sweep.csv    % positive for a range of fixed thresholds (both measures)
  hes1_distributions_*.png    per-nucleus histograms per time point (ratio, ring- and global-corrected)
  hes1_distributions_*.csv    the plotted histograms bin by bin (% of nuclei per bin, per time x condition);
                              *_wide.csv has one column per group, e.g. for Excel or Prism

Usage:
  python hes1_analysis.py All_Cells_Hes_slim.csv --out results_hes1
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from common import add_design_columns, mean_sd, histogram_table, save_histogram_tables, TIME_ORDER, COND_ORDER

COLORS = {"control": "#9a9a9a", "act": "#2b6cb0", "inh": "#e07b39"}
NAMES = {"control": "Control", "act": "Activated", "inh": "Inhibited"}


def gmm_boundary(log_ratio, seed=0):
    """Boundary between the two components of a 2-component Gaussian mixture (posterior = 0.5)."""
    x = log_ratio.reshape(-1, 1)
    g = GaussianMixture(2, random_state=seed, n_init=5).fit(x)
    lo, hi = np.sort(g.means_.ravel())
    grid = np.linspace(lo, hi, 4000).reshape(-1, 1)
    p_hi = g.predict_proba(grid)[:, np.argmax(g.means_.ravel())]
    return float(np.exp(grid[np.argmin(np.abs(p_hi - 0.5)), 0]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--out", default="results_hes1")
    ap.add_argument("--inh-percentile", type=float, default=95.0)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = add_design_columns(pd.read_csv(args.table))
    has_ring = {"BgRing_Median", "Mean_minus_BgRing"} <= set(d.columns)
    has_global = "Mean_minus_BgGlobal" in d.columns
    if not (has_ring or has_global):
        raise SystemExit("No background-corrected columns found (expected Mean_minus_BgRing and/or Mean_minus_BgGlobal).")
    has_inh = (d.cond == "inh").any()
    print(f"{len(d):,} nuclei in {d.SourceFile.nunique()} images; "
          f"ring background: {'yes' if has_ring else 'no'}, global background: {'yes' if has_global else 'no'}")

    # positivity (optional extras)
    if has_ring:
        d["ratio"] = d["Mean"] / d["BgRing_Median"]
        ratios = d["ratio"].dropna()
        ratios = ratios[ratios > 0]
        if len(ratios) > 10:
            print(f"Ratio boundary from 2-component GMM: {gmm_boundary(np.log(ratios.values)):.3f} (1.2 used in the paper)")
        d["pos_ratio"] = d["ratio"] >= 1.2
    thr_global = None
    if has_global and has_inh:
        thr_global = float(np.nanpercentile(d.loc[d.cond == "inh", "Mean_minus_BgGlobal"], args.inh_percentile))
        print(f"Global-corrected threshold ({args.inh_percentile:g}th pct of inhibited, pooled): {thr_global:.0f} a.u.")
        d["pos_global"] = d["Mean_minus_BgGlobal"] > thr_global

    # image level
    keys = ["time", "cond", "img"]
    aggs = {"SourceFile": ("SourceFile", "first"), "n_nuclei": ("Mean", "size"), "hes1_raw_mean": ("Mean", "mean")}
    if has_ring:
        aggs["hes1_ring_corrected_mean"] = ("Mean_minus_BgRing", "mean")
        aggs["pct_pos_ratio"] = ("pos_ratio", lambda s: 100 * s.mean())
    if has_global:
        aggs["hes1_corrected_mean"] = ("Mean_minus_BgGlobal", "mean")
    if thr_global is not None:
        aggs["pct_pos_global"] = ("pos_global", lambda s: 100 * s.mean())
    img = d.groupby(keys, observed=True).agg(**aggs).reset_index()
    img.to_csv(os.path.join(args.out, "hes1_image_level.csv"), index=False)

    order = ["hes1_ring_corrected_mean", "hes1_corrected_mean", "hes1_raw_mean", "pct_pos_ratio", "pct_pos_global", "n_nuclei"]
    rows = []
    for col in [c for c in order if c in img.columns]:
        s_ = mean_sd(img.set_index(["time", "cond"])[col], ["time", "cond"]).reset_index()
        s_.insert(0, "measure", col)
        rows.append(s_)
    summary = pd.concat(rows)
    summary.to_csv(os.path.join(args.out, "hes1_group_summary.csv"), index=False)
    main_col = "hes1_ring_corrected_mean" if has_ring else "hes1_corrected_mean"
    print(f"{'Ring' if has_ring else 'Global'} background-corrected nuclear Hes1 (a.u.), mean +- SD over images:")
    print(summary[summary.measure == main_col].round(1).to_string(index=False))

    # threshold sweep
    sweeps = []
    if has_ring:
        sweeps.append(("ratio", "ratio", [1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 2.0]))
    if has_global:
        ths = [300, 400, 500, 700, 900] + ([round(thr_global)] if thr_global is not None else [])
        sweeps.append(("global_corrected", "Mean_minus_BgGlobal", sorted(set(ths))))
    sweep = []
    for measure, col, ths in sweeps:
        for th in ths:
            p = (d.assign(p=d[col] >= th).groupby(keys, observed=True).p.mean() * 100).groupby(["time", "cond"], observed=True).mean()
            for (t, c), v in p.items():
                sweep.append({"measure": measure, "threshold": th, "time": t, "cond": c, "pct_positive": v})
    pd.DataFrame(sweep).to_csv(os.path.join(args.out, "hes1_threshold_sweep.csv"), index=False)

    # distributions
    plots = []
    if has_ring:
        plots.append(("ratio", "ratio", np.arange(0.8, 2.65, 0.1), 1.2, "Nuclear mean / ring background"))
        plots.append(("ring", "Mean_minus_BgRing", np.arange(-400, 1901, 100), None, "Nuclear mean - ring background (a.u.)"))
    if has_global:
        plots.append(("global", "Mean_minus_BgGlobal", np.arange(-200, 2101, 100), thr_global, "Nuclear mean - global background (a.u.)"))
    times = [t for t in TIME_ORDER if (d.time == t).any()]
    for measure, col, bins, cut, xlabel in plots:
        fig, axs = plt.subplots(1, len(times), figsize=(3.4 * len(times), 3), sharey=True, squeeze=False)
        hist_rows = []
        for ax, t in zip(axs[0], times):
            for c in COND_ORDER:
                x = d.loc[(d.time == t) & (d.cond == c), col].dropna()
                if len(x) == 0:
                    continue
                h = histogram_table(x, bins)
                hist_rows.append(h.assign(measure=measure, time=t, cond=c, n=len(x)))
                ax.plot(h.bin_center, h.pct, color=COLORS[c], label=f"{NAMES[c]} ({len(x):,})")
            if cut is not None:
                ax.axvline(cut, ls="--", color="#777777", lw=0.8)
            ax.set_title(t, loc="left"); ax.set_xlabel(xlabel); ax.legend(frameon=False, fontsize=7)
        axs[0][0].set_ylabel("% of nuclei")
        fig.tight_layout(); fig.savefig(os.path.join(args.out, f"hes1_distributions_{measure}.png"), dpi=200); plt.close(fig)
        if hist_rows:
            long_df = pd.concat(hist_rows, ignore_index=True)[
                ["measure", "time", "cond", "n", "bin_start", "bin_end", "bin_center", "count", "pct"]]
            save_histogram_tables(long_df, ["time", "cond"], os.path.join(args.out, f"hes1_distributions_{measure}"))
    print(f"Results written to {args.out}/")


if __name__ == "__main__":
    main()
