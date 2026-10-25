# Python analysis

Two scripts analyse the per-cell tables written by the Fiji workflows:

| Script | Input | Analysis |
|---|---|---|
| `hes1_analysis.py` | `All_Cells_Hes_slim.csv` (Workflow 2) | Nuclear Hes1 per image and per nucleus, background-corrected; Hes1-positive nuclei |
| `elongation_analysis.py` | `All_Cells_Shape_slim.csv` (Workflow 1), one per Cellpose model | Cell elongation per image and per cell; optional nuclei-only filter; comparison between models |

`common.py` holds shared helpers (reading time point, condition and image number from file names; histograms) and `requirements.txt` lists the Python packages.

## Input files

Each Fiji workflow writes `All_Cells_<type>.csv` and `All_Cells_<type>_slim.csv` into its Results folder at the end of every run. Use the `_slim.csv` file. To compare Cellpose models, give the elongation script one table per model (one per Results folder).

Time point and condition are read from the image names, which must look like `<time>_notch_<condition>_..._Image_<n>`, e.g. `day5_notch_act_cd31AF555_Hes1AF647_Image_2`:

- time point: `12h`, `48h` or `day5` at the start of the name
- condition: the word after `notch_` (`control`, `act` or `inh`; capitals are ignored)
- `_Image_<n>` is optional; if any name lacks it, images are numbered 1, 2, 3, ... within each time point and condition, in name order

To use other names or conditions, edit the patterns at the top of `common.py`.

## Setup in Visual Studio Code (Windows)

Done once:

1. **Install Python** 3.10 or newer from [python.org](https://www.python.org/downloads/). In the installer, tick *Add python.exe to PATH*.
2. **Install VS Code** and, in it, the *Python* extension by Microsoft (Extensions panel, `Ctrl+Shift+X`).
3. **Open the repository folder**: *File › Open Folder*, and choose the folder that contains this repository's main `README.md` (not just a single file).
4. **Trust the folder**: if the status bar shows *Restricted Mode*, click it and choose *Trust*. Otherwise VS Code keeps the Python extension switched off.
5. **Use a Windows terminal**: press `Ctrl+Shift+P`, run *Terminal: Select Default Profile* and choose *PowerShell* (or *Command Prompt*). Windows Subsystem for Linux (WSL) is not needed.
6. **Create an environment and install the packages**. Open a terminal (*Terminal › New Terminal*) and run:

   ```
   py -m venv .venv
   .venv\Scripts\activate
   python -m pip install -r analysis\requirements.txt
   ```

   The prompt now starts with `(.venv)`. If activation fails with *running scripts is disabled on this system*, run `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` once and activate again, or use a Command Prompt terminal and `.venv\Scripts\activate.bat`. If `py` is not recognised, use `python`.

   Alternatively, with the Python extension active: `Ctrl+Shift+P` › *Python: Create Environment* › *Venv*, and tick `analysis\requirements.txt`.

Each time you open a new terminal later, activate the environment again with `.venv\Scripts\activate`.

## Running the analysis

Run the commands from the repository folder in the VS Code terminal, with `(.venv)` active. Put paths that contain spaces in quotes; in Windows Explorer, Shift + right-click › *Copy as path* copies a file's path.

**Hes1**

```
python analysis\hes1_analysis.py "path\to\Results_Hes\All_Cells_Hes_slim.csv" --out results_hes1
```

| Option | Meaning |
|---|---|
| `--out` | Output folder (created if needed) |
| `--inh-percentile` | Percentile of inhibited nuclei used as the global-corrected positivity threshold (default 95) |

**Elongation**

```
python analysis\elongation_analysis.py "path\to\Results_cpsam\All_Cells_Shape_slim.csv" --names cpsam --out results_shape
```

| Option | Meaning |
|---|---|
| `--names` | Label per input table, in the same order (default: the `ResultsFolder` column) |
| `--out` | Output folder (created if needed) |
| `--min-area` | Remove segments with area ≤ this value (µm²) as nuclei-only, e.g. `--min-area 434` (1.2 × the 95th percentile of nuclear area in our data) |
| `--nuclei-table` | Instead of `--min-area`: take the threshold from the nuclear areas in a Workflow 2 table (percentile × factor, see below) |
| `--nuclei-percentile` | Percentile for `--nuclei-table` (default 95) |
| `--nuclei-factor` | Multiplies that percentile, to also catch a nucleus with a little cytoplasm (default 1.2; use 1 for the plain percentile) |
| `--exclude` | Image names (`SourceFile`) to leave out entirely, e.g. images with failed segmentation |

To compare Cellpose models, give one table per model:

```
python analysis\elongation_analysis.py "path\to\Results_cpsam\All_Cells_Shape_slim.csv" "path\to\Results_cpdino\All_Cells_Shape_slim.csv" --names cpsam cpdino --out results_models
```

Each script prints a short summary in the terminal and ends with *Results written to ...*. With only one image per group, standard deviations are empty, and some analyses are skipped when the data don't allow them (for example the inhibited-based threshold without inhibited images, or the model comparison with a single model).

## Output files

All tables are CSV files. Groups are identified by `time` and `cond` (and `model` for elongation). Images are identified by `img` (image number) and `SourceFile` (image name). Image-level values are summaries over all cells or nuclei of an image; group values are the mean ± SD over the images of a group.

### `hes1_analysis.py`

| File | Content |
|---|---|
| `hes1_image_level.csv` | One row per image: `n_nuclei`; mean nuclear Hes1 over the image's nuclei, raw (`hes1_raw_mean`), minus ring background (`hes1_ring_corrected_mean`) and minus global background (`hes1_corrected_mean`); % Hes1-positive nuclei with ratio ≥ 1.2 to the ring background (`pct_pos_ratio`) and above the inhibited-based threshold (`pct_pos_global`) |
| `hes1_group_summary.csv` | One row per measure and group: `mean`, `sd` and `n_images` over the images of the group, for each column of the image-level table |
| `hes1_threshold_sweep.csv` | % positive nuclei (`pct_positive`) per group for a range of fixed thresholds (`threshold`), on the ratio and on the global-corrected measure (`measure`) |
| `hes1_distributions_<measure>.png` | Per-nucleus histograms, one panel per time point, one curve per condition. `<measure>` is `ratio` (nuclear mean ÷ ring background), `ring` (minus ring background) or `global` (minus global background); the dashed line is the positivity threshold |
| `hes1_distributions_<measure>.csv` | The histograms of the PNG bin by bin: `n` (nuclei in the group), `bin_start`, `bin_end`, `bin_center`, `count` and `pct` (% of the group's nuclei in the bin) |
| `hes1_distributions_<measure>_wide.csv` | The same `pct` values with one column per group (header includes `n`), for plotting in Excel or Prism |

The terminal also prints the two positivity thresholds computed from the data: the global-corrected threshold (95th percentile of inhibited nuclei) and the boundary between negative and positive nuclei from a two-component Gaussian mixture fitted to the ratios.

### `elongation_analysis.py`

Elongation is the aspect ratio of the equivalent ellipse (`Ellipse.Elong`; 1 = round).

| File | Content |
|---|---|
| `elongation_image_level.csv` | One row per image and model: `n_cells`, `median_ar` and `mean_ar` (aspect ratio), `pct_ar_ge3` (% of cells with aspect ratio ≥ 3), `median_area` |
| `elongation_group_summary.csv` | One row per measure, group and model: `mean`, `sd` and `n_images` over the images, for `median_ar`, `pct_ar_ge3`, `n_cells` and `median_area` |
| `elongation_percentiles.csv` | Per-cell aspect-ratio percentiles (10th, 25th, 50th, 75th, 90th) per group and model, pooling all cells |
| `elongation_distributions.png` | Per-cell aspect-ratio histograms of the first model, one panel per time point, one curve per condition; the last bin includes 6 and above |
| `elongation_distributions.csv` | The histograms bin by bin for **every** model: `n`, `bin_start`, `bin_end`, `bin_center`, `count`, `pct` |
| `elongation_distributions_wide.csv` | The same `pct` values with one column per model and group, for Excel or Prism |
| `model_agreement.csv` | Only with two or more models: per aspect-ratio range (`ar_bin`) of the first model's cells, the share (`pct_unmatched`) that has no matching cell in the other model (`other`). Cells match when their centroids are mutual nearest neighbours within 0.3× the cell radius and their areas are within 30% |
| `filter_summary.csv` | Only with `--min-area` or `--nuclei-table`: per image and model, segments before filtering (`n_segments`), removed as nuclei-only (`n_removed`) and `pct_removed` |

In the distribution tables, a bin with `bin_start` of `-inf` or `bin_end` of `inf` also holds all values beyond the plotted range (shown as e.g. *6+* in the figures).

## Troubleshooting

| Message | Fix |
|---|---|
| `python` / `py` *is not recognized* | Python is not installed or not on the PATH: re-run the installer with *Add python.exe to PATH* ticked |
| No `(.venv)` in the prompt, or *No module named pandas* | Activate the environment (`.venv\Scripts\activate`) or select it via `Ctrl+Shift+P` › *Python: Select Interpreter* |
| *Could not parse time/condition from: ...* | An image name lacks the time point or `notch_<condition>_`; rename it or edit the patterns in `common.py` |
| *No such file or directory* | Check the path to the CSV file, and put it in quotes if it contains spaces |
| A terminal window asks to install or update *Windows Subsystem for Linux* | Press Esc and set PowerShell as the default terminal (setup step 5) |
