// @String(label="Run mode", choices=("Segment", "Update"), style="list", value="Segment") runMode
// @String(label="Image selection", choices=("Single open image", "Single file (.lif or .tif)", "Whole folder"), style="list") processMode
// @String(label="File extension(s), comma-separated (for Single file / Whole folder), e.g. .lif or .tif, .tiff", value=".lif") fileTypeChoice
// @File(label="Image file, .lif or .tif (ignored for 'Single open image'; for 'Whole folder', pick any file inside the target folder)", style="open", required=false) file_name
// @File(label="Cellpose conda env folder (required for Segment mode)", style="directory", required=false, value="C:\\ProgramData\\anaconda3\\envs\\cellpose3") cellposeEnvPath
// @File(label="Cellpose model_path folder (required for Segment mode)", style="directory", required=false, value="C:\\Fiji\\fiji-latest-win64-jdk\\Fiji") cellposeModelPath
//
// =============================================================================
//  Workflow 2 : Hes images  -  Nuclear segmentation + Hes intensity quantification
// =============================================================================
//
//  WHAT THIS SCRIPT DOES (per image):
//    1. Duplicates the Hes + DAPI channels from a 5-channel .lif or .tif image
//    2. Builds a 2-channel composite with fixed, comparable display ranges
//    3. Segmentation - one of two modes (chosen from a dialog at start):
//         - Segment: runs Cellpose ("nuc" model, DAPI channel only) and
//           converts the resulting nuclei label image to ROIs, then SAVES
//           those ROIs as <prefix>_RoiSet.zip so they can be manually
//           corrected later.
//         - Update: skips Cellpose entirely. Looks for
//           <prefix>_RoiSet_manual.zip first (a manually-corrected copy you
//           made yourself); if it doesn't exist, falls back to the original
//           <prefix>_RoiSet.zip. A label image is then reconstructed from
//           those ROIs so the intensity color-coding step still works.
//           MANUAL CORRECTION STEP (done by you, in Fiji, between runs):
//             a) Run in Segment mode once to produce <prefix>_RoiSet.zip
//             b) Open that zip into the ROI Manager (More >> Open...)
//             c) Correct/add/delete ROIs as needed
//             d) ROI Manager > More >> Save... and save it as EXACTLY
//                <prefix>_RoiSet_manual.zip in the same Results folder
//             e) Re-run this script in Update mode - it will pick up your
//                corrected ROIs automatically for that image.
//    4. Measures Hes (not DAPI) channel intensity inside every nuclear ROI
//       (Mean, StdDev, Min/Max, Integrated Density, Median, Area) via plain
//       ROI Manager "Measure" - no shape descriptors are needed here, so
//       there is no MorphoLibJ "Analyze Regions" step (unlike Workflow 1).
//    4b. Estimates the Hes BACKGROUND outside the nuclei in two ways and adds
//        background-corrected columns to the per-nucleus table:
//         - GLOBAL (one value per image, main correction): median (and mean)
//           of all Hes pixels outside the nuclei, after growing every nucleus
//           by bgGapPx pixels so blurry edges / slightly under-segmented
//           masks don't leak nuclear signal into the background.
//         - RING (one value per nucleus, extra comparison column): median (and
//           mean) of the Hes pixels in a band around each nucleus, from
//           bgGapPx to (bgGapPx + bgRingWidthPx) pixels outside its edge, with
//           ALL nuclei (+ the same gap) excluded, so neighbouring nuclei never
//           count as background. Rings with fewer than bgRingMinPixels usable
//           pixels (e.g. in very dense areas) are reported as NaN.
//        The MEDIAN is used for the subtraction (robust to bright debris and
//        missed nuclei); the MEAN is recorded alongside for reference.
//        New per-nucleus columns:
//          BgGlobal_Median, BgGlobal_Mean,
//          Mean_minus_BgGlobal       = Mean - BgGlobal_Median
//          RawIntDen_minus_BgGlobal  = RawIntDen - BgGlobal_Median * nPixels
//          BgRing_Median, BgRing_Mean, BgRing_NPixels,
//          Mean_minus_BgRing         = Mean - BgRing_Median
//          RawIntDen_minus_BgRing    = RawIntDen - BgRing_Median * nPixels
//        Negative corrected values are kept (not clipped to 0) on purpose, so
//        the distribution of dim nuclei is not biased upward.
//    5. Saves:
//         A) <prefix>_Cells.csv        - per-nucleus Hes intensity measurements
//         B) <prefix>_Composite.tif    - Hes+DAPI composite (fixed LUT/range)
//            <prefix>_Composite_flat.tif - flattened (RGB) version of the same
//         C) <prefix>_Overlay.tif      - ROIs (no labels) flattened on composite
//         E) <prefix>_HesOnlyOverlay.tif - same ROIs flattened on the Hes channel
//              alone (extra QC image, no DAPI channel to distract from it)
//         D) <prefix>_HesIntensity.tif + <prefix>_HesIntensity_flat.tif -
//              Hes-intensity color-coded label image (Fire LUT + calibration
//              bar), flattened version included
//         F) <prefix>_BgQC.tif - background QC image on the Hes channel:
//              shaded = nuclei + gap (EXCLUDED from all background estimates),
//              yellow outlines = boundary of the pixels actually used for the
//              RING backgrounds (neighbouring nuclei + gap are cut out, so the
//              outlines never cross another nucleus). Everything unshaded =
//              pixels used for the GLOBAL background.
//         (Segment mode only) <prefix>_RoiSet.zip - the ROI set itself
//    6. Appends one row (mean/min/max of the columns listed in summaryColumns,
//       NaN values ignored) to a single run-wide summary table, plus per-image
//       background values: BgGlobal_Median, BgGlobal_Mean, BgGlobal_NPixels,
//       BgGlobal_AreaFraction (fraction of the image used as global background
//       - low values = dense image) and N_RingsEmpty.
//
//  NOTE ON DENSE vs SPARSE IMAGES:
//    In sparse images the GLOBAL background is dominated by cell-free dish;
//    in dense images it is mostly cytoplasm, which can carry some Hes signal.
//    BgGlobal_AreaFraction is recorded so you can check whether BgGlobal_Median
//    drifts with density. The RING background is always perinuclear and is
//    therefore more consistent between dense and sparse fields.
//
//  At the end of the run, the script also writes:
//    - Summary_Cells.csv          (one row per image, whole run)
//    - All_Cells_Hes.csv + All_Cells_Hes_slim.csv  (every <image>_Cells.csv in the
//      Results folder combined, with SourceFile + ResultsFolder columns;
//      ready for the Python scripts in analysis/)
//    - Workflow2_Hes_Parameters.txt  (all parameters used + timestamp)
//
//  RUN MODES (all chosen from ONE dialog when you click Run):
//    A) Single already-open image   (leave "Image file" blank)
//    B) Single file (.lif or .tif)   (pick the file in the dialog; all series
//                                     inside it are opened and processed; a .tif
//                                     is one image, named after its file)
//    C) Whole folder                 (pick ANY matching *Hes*.lif or .tif file that
//                                     lives in the target folder - the script
//                                     then processes every matching file in
//                                     that file's parent folder)
//
//  The SAME dialog also lets you choose Segment vs Update mode, and (for
//  Segment mode) the Cellpose conda-env and model_path folders, so nothing
//  needs to be edited in the code just to run the script day to day. Fields
//  not needed for your choice (e.g. the image file in "Single open image" mode,
//  or the Cellpose folders in "Update" mode) can be left blank.
//
//  RESULTS LOCATION:
//    A "Results_Hes" (variable: resultsSubfolderName) subfolder is created
//    next to the source file(s) (mode A falls back to a folder you choose if
//    the open image has no file on disk).
//
//  REQUIRED DEPENDENCIES (enable via Help > Update... > Manage update sites):
//    - BIOP update site        -> Cellpose wrapper ("Cellpose ...") and
//                                  "Label image to ROIs"
//    - IJPB-plugins (MorphoLibJ) update site -> LabelToValuePlugin
//                                  ("Assign Measure to Label")
//    - Bio-Formats              -> bundled with Fiji, used to read .lif files (.tif files open directly)
//    - A local Cellpose Python/conda environment (ending in "cellpose3") with
//      the "nuc" (nuclei) model installed and reachable at cellposeEnvPath.
//      Test it once via Plugins > BIOP > Cellpose before running this script.
//
//  NOTE ON DIFFERENCES FROM WORKFLOW 1:
//    - Cellpose here segments a SINGLE channel (DAPI) with the "nuc" model,
//      so ch2 is always 0 in the Cellpose call (no second channel).
//    - No shape descriptors are computed (no MorphoLibJ "Analyze Regions") -
//      Workflow 1 needed it to fix an unreliable ROI-based elongation metric,
//      but plain intensity statistics from ROI Manager "Measure" are reliable
//      as-is, so the simpler/faster path is used here.
//    - Cellpose IS called with a diameter parameter (cellposeDiameter): it
//      still divides by this value internally, so 0 crashes with a
//      ZeroDivisionError. Nuclei are usually smaller than whole cells -
//      adjust to your typical nucleus diameter in pixels if segmentation
//      quality is off.
//
// =============================================================================


// =============================================================================
//  1) USER-EDITABLE PARAMETERS
// =============================================================================

scriptName            = "Workflow2_Hes";         // used for the parameters file name
resultsSubfolderName  = "Results_Hes";           // output subfolder name (no spaces)
fileSuffix            = "Hes";                   // used to find files in folder mode (must appear in filename)
// Input file types come from the "File extension(s)" field in the dialog
// (fileTypeChoice), e.g. ".lif", ".tif, .tiff" or ".lif, .tif". Case-insensitive;
// a missing leading dot is added. A .lif is processed series by series; any other
// file (e.g. .tif) is opened as one image.
fileExtensions = parseExtensions(fileTypeChoice);

// Combined per-cell tables written into the Results folder at the end of every run
// (all <image>_Cells.csv files in that folder, input for the Python scripts in analysis/)
saveCombinedTables    = true;
combinedTableName     = "All_Cells_Hes";
combinedSlimCols      = newArray("Label", "Area", "Mean", "Median", "RawIntDen",
    "BgGlobal_Median", "BgGlobal_Mean", "Mean_minus_BgGlobal", "RawIntDen_minus_BgGlobal",
    "BgRing_Median", "BgRing_Mean", "BgRing_NPixels", "Mean_minus_BgRing", "RawIntDen_minus_BgRing");

// ---- Channel roles (1-based, as in the ORIGINAL 5-channel image) -----------
CH_HES  = 3;
CH_DAPI = 5;

// ---- Display settings for the 2-channel composite ---------------------------
hesColor  = "Green";
hesMin    = 0;
hesMax    = 5000;

dapiColor = "Blue";
dapiMin   = 0;
dapiMax   = 20000;

// ---- Cellpose settings -------------------------------------------------------
// cellposeEnvPath and cellposeModelPath are set from the dialog (see the @File
// parameters at the very top of this file), not here.
cellposeEnvType      = "conda";
cellposeModel        = "nuc";  // nuclei model - segments the DAPI channel only
cellposeDiameter     = 40.0;   // cellpose still divides by this value internally (0 causes
                                // a ZeroDivisionError) - default value used;
                                // adjust to your typical nucleus diameter in pixels, but do
                                // not set it to 0.
cellposeExtraArgs    = "";     // extra "key=value " flags appended verbatim if needed

// ---- Per-nucleus intensity measurements (plain ROI Manager "Measure" on the Hes
// channel - no shape descriptors needed for this workflow) ----------------------
measurementsString = "area mean standard min integrated median display redirect=None decimal=3";

// ---- Hes BACKGROUND estimation (see header, step 4b) ----------------------------
bgGapPx          = 3;  // pixels each nucleus is grown by before it is excluded from the
                       // background (protects against blurry edges / under-segmentation)
bgRingWidthPx    = 5;  // width (pixels) of the per-nucleus ring, starting bgGapPx outside
                       // the nuclear edge
bgRingMinPixels  = 10; // rings with fewer usable pixels than this are reported as NaN
saveBgQcImage    = true; // save <prefix>_BgQC.tif showing the background regions
bgQcExcludedFill = "#5500FFFF"; // ARGB fill for the excluded (nuclei + gap) zone in the QC image
bgQcRingColor    = "yellow";    // outline color of the ring regions in the QC image

// ---- Columns included in the run-wide SUMMARY table (mean/min/max of each) -----
// The per-nucleus CSV always keeps every measured column; only THIS list controls
// what also gets aggregated into Summary_Cells.csv. Edit freely - names must match
// the per-cell CSV's column headers exactly.
summaryColumns = newArray("Area", "Mean", "StdDev", "Min", "Max", "IntDen", "RawIntDen", "Median",
    "Mean_minus_BgGlobal", "RawIntDen_minus_BgGlobal",
    "BgRing_Median", "BgRing_Mean", "Mean_minus_BgRing", "RawIntDen_minus_BgRing");

// ---- Hes intensity color-coding settings ---------------------------------------
intensityColorColumn = "Mean"; // column from the per-nucleus table used for the
                                // color-coded image below (e.g. "Mean_minus_BgGlobal"
                                // for background-corrected values - note these can be
                                // negative, so you may want hesColorMin below 0)
hesColorMin = 0;
hesColorMax = 5000; // adjust to your typical Hes intensity range

calBarLocation    = "Upper Right";
calBarFill        = "White";
calBarLabelColor  = "Black";
calBarNumberLabels = 5;
calBarDecimal     = 0;
calBarFontSize    = 12;
calBarZoom        = 4;   // size of the calibration bar on the color-coded images (2 = previous size)

// ---- Overlay images (cell / nucleus outlines) --------------------------------
outlineWidth      = 4;          // line width (pixels) of the ROI outlines in the overlay images
// Scale bar burned into the overlay images (lower right). Length is in the image's
// calibrated units (um); it is skipped, with a warning, if the image has no calibration.
scaleBarLength    = 100;        // um
scaleBarThickness = 10;         // pixels
scaleBarFontSize  = 40;         // pixels; set 0 to hide the length label
scaleBarColor     = "White";
scaleBarLocation  = "Lower Right";

// ---- Output file naming -------------------------------------------------------
cellsCsvSuffix         = "_Cells.csv";
compositeTifSuffix     = "_Composite.tif";
compositeFlatTifSuffix = "_Composite_flat.tif";
overlayTifSuffix       = "_Overlay.tif";
hesOnlyOverlayTifSuffix = "_HesOnlyOverlay.tif";
hesColorTifSuffix      = "_HesIntensity.tif";
hesColorFlatTifSuffix  = "_HesIntensity_flat.tif";
bgQcTifSuffix          = "_BgQC.tif";
roiZipSuffix           = "_RoiSet.zip";        // saved automatically after every Segment run
roiZipManualSuffix     = "_RoiSet_manual.zip"; // YOU create this: open roiZipSuffix in the ROI
                                                // Manager, correct it, and "Save" it under this
                                                // exact name (same folder, same prefix). Update
                                                // mode looks for this file first.
summaryCsvName       = "Summary_Cells.csv";
parametersTxtName    = scriptName + "_Parameters.txt";


// =============================================================================
//  2) GLOBALS USED DURING A RUN (do not edit)
// =============================================================================

summaryRowIndex = 0;
runResultsDir   = "";


// =============================================================================
//  3) ENTRY POINT - dispatch based on the dialog shown at Run time
// =============================================================================

print("\\Clear"); // clear the Log window so each run starts with a clean log
getDateAndTime(y0, mo0, dow0, d0, h0, mi0, s0, ms0);
print("=== " + scriptName + " - starting run " +
      y0 + "-" + pad2(mo0 + 1) + "-" + pad2(d0) + " " +
      pad2(h0) + ":" + pad2(mi0) + ":" + pad2(s0) + " ===");
print("Run mode (Segment/Update): " + runMode);
print("Image selection: " + processMode);
print("File extension(s): " + String.join(fileExtensions, ", "));
print("file_name from dialog: [" + file_name + "]");

processingMode = runMode; // "Segment" or "Update", harvested from the dialog above

if (processingMode == "Segment") {
    if (isUnset(cellposeEnvPath) || isUnset(cellposeModelPath)) {
        exit("Segment mode requires both 'Cellpose conda env folder' and " +
             "'Cellpose model_path folder' to be set in the dialog.");
    }
}

// --- start-of-run cleanup -------------------------------------------------------
if (!isOpen("ROI Manager")) {
    run("ROI Manager...");
}
roiManager("reset");
run("Clear Results");
if (isOpen("Summary")) {
    close("Summary");
}
Table.create("Summary");
summaryRowIndex = 0;

if (processMode != "Single open image") {
    // Only safe to close every open image when we are NOT reusing an
    // already-open image as the input (that would defeat "Single open image" mode).
    if (nImages > 0) {
        print("Closing " + nImages + " currently open image(s) before starting...");
        close("*");
    }
}

if (processMode == "Single open image") {
    runSingleOpenImageMode();
} else if (processMode == "Single file (.lif or .tif)") {
    runSingleFileMode(file_name);
} else {
    runFolderBatchMode(file_name);
}

print("=== " + scriptName + " - run complete ===");


// =============================================================================
//  4) MODE FUNCTIONS
// =============================================================================

function runSingleOpenImageMode() {
    if (nImages == 0) {
        exit("No image is open. Open the raw Hes image first, or choose a different Image selection mode.");
    }
    imgTitle = getTitle();

    // Try to recover the source file location for naming/output purposes.
    dir = getDirectory("image");
    fname = getInfo("image.filename");
    if (dir == "" || fname == "") {
        dir = getDirectory("Choose a folder to save results into");
        fname = imgTitle;
    }
    lifBaseName = File.getNameWithoutExtension(fname);
    lifSep = indexOf(lifBaseName, ".lif - ");     // TIF exported from a .lif series
    if (lifSep > 0) lifBaseName = substring(lifBaseName, 0, lifSep);

    runResultsDir = ensureResultsDir(dir);
    writeParametersFile(runResultsDir);

    baseName = getBaseName(imgTitle, lifBaseName);
    processSingleImage(imgTitle, lifBaseName, baseName, runResultsDir);

    finalizeSummary(runResultsDir);
    writeCombinedCellsTables(runResultsDir);
}

function runSingleFileMode(path) {
    print("DEBUG: file_name from dialog = [" + path + "]");
    if (isUnset(path) || !File.exists(path)) {
        print("  No valid image file came from the dialog - asking you to pick one now...");
        path = File.openDialog("Select a " + fileSuffix + " image file (" + String.join(fileExtensions, ", ") + ")");
    }
    if (isUnset(path) || !File.exists(path)) exit("No valid image file was selected.");
    if (!hasImageExtension(path)) {
        exit("Unsupported file type: " + File.getName(path) + "\nSupported: " + String.join(fileExtensions, ", "));
    }

    dir = File.getParent(path) + File.separator;
    if (isUnset(dir) || dir == File.separator) exit("Could not determine the parent folder of: [" + path + "]");
    runResultsDir = ensureResultsDir(dir);
    writeParametersFile(runResultsDir);

    processImageFile(path, runResultsDir);

    finalizeSummary(runResultsDir);
    writeCombinedCellsTables(runResultsDir);
}

function runFolderBatchMode(path) {
    print("DEBUG: file_name from dialog = [" + path + "]");
    if (isUnset(path) || !File.exists(path)) {
        print("  No valid file came from the dialog - asking you to pick the target folder now...");
        dir = getDirectory("Choose the folder containing the " + fileSuffix + " image files");
        if (isUnset(dir)) exit("No folder selected.");
    } else {
        dir = File.getParent(path) + File.separator;
        if (isUnset(dir) || dir == File.separator) exit("Could not determine the parent folder of: [" + path + "]");
    }

    runResultsDir = ensureResultsDir(dir);
    writeParametersFile(runResultsDir);

    list = getFileList(dir);
    nProcessed = 0;
    nLif = 0;
    nTif = 0;
    matchingFiles = newArray(0);
    for (i = 0; i < list.length; i++) {
        fname = list[i];
        if (hasImageExtension(fname) && indexOf(fname, fileSuffix) >= 0) {
            matchingFiles = Array.concat(matchingFiles, fname);
            if (endsWith(toLowerCase(fname), ".lif")) nLif++;
            else nTif++;
        }
    }
    print("Found " + matchingFiles.length + " matching *" + fileSuffix + "* file(s) in " + dir +
          "  (" + nLif + " .lif, " + nTif + " other)");
    if (nLif > 0 && nTif > 0) {
        print("WARNING: this folder has both .lif and other image files. If those were exported from");
        print("         these .lif files, the same images are processed twice and their results");
        print("         overwrite each other. Keep only one of the two kinds in the folder.");
    }
    for (i = 0; i < matchingFiles.length; i++) {
        print("--- File " + (i + 1) + " of " + matchingFiles.length + ": " + matchingFiles[i] + " ---");
        processImageFile(dir + matchingFiles[i], runResultsDir);
        nProcessed++;
    }
    if (nProcessed == 0) {
        print("WARNING: no *" + fileSuffix + "* files with extension " + String.join(fileExtensions, ", ") + " found in " + dir);
    }

    finalizeSummary(runResultsDir);
    writeCombinedCellsTables(runResultsDir);
}


// =============================================================================
//  5) CORE PROCESSING
// =============================================================================

// True if 'name' ends with one of the supported extensions (fileExtensions).
function hasImageExtension(name) {
    lname = toLowerCase(name);
    for (k = 0; k < fileExtensions.length; k++) {
        if (endsWith(lname, toLowerCase(fileExtensions[k]))) return true;
    }
    return false;
}

// Opens one input file (.lif with all its series, or a single .tif/.tiff image)
// and processes every image in it.
function processImageFile(path, resultsDir) {
    if (endsWith(toLowerCase(path), ".lif")) {
        processLifFile(path, resultsDir);
    } else {
        processTifFile(path, resultsDir);
    }
}

// Opens a TIFF image (e.g. one series exported from a .lif with File > Save As > Tiff,
// keeping the channel order and calibration) and processes it. The file name
// (without extension) is used as the image prefix in all outputs.
function processTifFile(path, resultsDir) {
    fileBaseName = File.getNameWithoutExtension(path);
    // Fiji names an exported .lif series "<lif name>.lif - Image 2.tif": use the same
    // prefix as when the .lif itself is processed ("<lif name>_Image_2")
    lifSep = indexOf(fileBaseName, ".lif - ");
    seriesName = "";
    if (lifSep > 0) {
        seriesName = substring(fileBaseName, lifSep + lengthOf(".lif - "));
        fileBaseName = substring(fileBaseName, 0, lifSep);
    }
    print("Processing file: " + path);
    titlesBefore = getList("image.titles");
    open(path);
    titlesAfter = getList("image.titles");
    openedTitles = arrayDiff(titlesAfter, titlesBefore);
    if (openedTitles.length == 0) {
        print("  WARNING: could not open " + path + " - skipped.");
        return;
    }
    for (i = 0; i < openedTitles.length; i++) {
        imgTitle = openedTitles[i];
        if (seriesName != "") baseName = sanitizeName(fileBaseName + "_" + seriesName);
        else baseName = sanitizeName(fileBaseName);
        if (openedTitles.length > 1) baseName = baseName + "_" + (i + 1);
        print("  [" + (i + 1) + "/" + openedTitles.length + "] " + imgTitle);
        processSingleImage(imgTitle, fileBaseName, baseName, resultsDir);
        if (isOpen(imgTitle)) close(imgTitle);
    }
}

// Opens every series of a .lif file, processes each individually, closes them.
function processLifFile(path, resultsDir) {
    lifBaseName = File.getNameWithoutExtension(path);
    print("Processing file: " + path);

    titlesBefore = getList("image.titles");
    run("Bio-Formats Importer",
        "open=[" + path + "] autoscale color_mode=Default rois_import=[ROI manager] " +
        "view=Hyperstack stack_order=XYCZT open_all_series");
    titlesAfter = getList("image.titles");

    openedTitles = arrayDiff(titlesAfter, titlesBefore);
    print("  " + openedTitles.length + " image(s)/series found in this file.");

    for (i = 0; i < openedTitles.length; i++) {
        imgTitle = openedTitles[i];
        baseName = getBaseName(imgTitle, lifBaseName);
        print("  [" + (i + 1) + "/" + openedTitles.length + "] " + imgTitle);
        processSingleImage(imgTitle, lifBaseName, baseName, resultsDir);
        if (isOpen(imgTitle)) close(imgTitle);
    }
}

// Full per-image pipeline: segmentation, measurement, and all exports.
function processSingleImage(imgTitle, lifBaseName, baseName, resultsDir) {
    print("  -> Image: " + imgTitle + "   (prefix: " + baseName + ")");

    run("Clear Results");
    roiManager("reset");

    // --- 1) duplicate ALL channels, then keep+reorder DAPI+Hes via Arrange Channels ---
    // NOTE: we deliberately do NOT use Duplicate's "channels=X-Y" range here. That
    // syntax keeps every channel BETWEEN X and Y inclusive - for adjacent channels
    // (like Workflow 1's 3-4) that happens to be just two channels, but for
    // non-adjacent channels like CH_HES=3/CH_DAPI=5 it would wrongly keep channels
    // 3, 4 AND 5 (three channels, with channel 4 landing in the middle where DAPI
    // was expected). "Arrange Channels..." with new=<DAPI><HES> instead selects
    // and reorders in one step: the output ALWAYS has exactly 2 channels, with
    // DAPI as channel 1 and Hes as channel 2, regardless of their original indices.
    selectWindow(imgTitle);
    cleanBaseName = sanitizeName(baseName);
    dupTitle = cleanBaseName + "_dup";
    if (isOpen(dupTitle)) close(dupTitle);
    run("Duplicate...", "duplicate title=" + dupTitle);
    selectWindow(dupTitle);
    newChannelOrder = d2s(CH_DAPI, 0) + d2s(CH_HES, 0); // e.g. CH_DAPI=5, CH_HES=3 -> "53"
    run("Arrange Channels...", "new=" + newChannelOrder);
    Stack.getDimensions(dupWidth, dupHeight, dupChannels, dupSlices, dupFrames);
    if (dupChannels != 2) {
        print("ERROR: Arrange Channels produced " + dupChannels + " channel(s) instead of 2 for " +
              imgTitle + " (new=" + newChannelOrder + ") - skipping this image.");
        if (isOpen(dupTitle)) close(dupTitle);
        return;
    }
    dapiIdx = 1;
    hesIdx  = 2;

    // --- 2) composite with fixed display ranges --------------------------------
    selectWindow(dupTitle);
    run("Make Composite");
    Stack.setChannel(hesIdx);
    run(hesColor);
    setMinAndMax(hesMin, hesMax);
    Stack.setChannel(dapiIdx);
    run(dapiColor);
    setMinAndMax(dapiMin, dapiMax);

    // --- B) save composite -------------------------------------------------------
    saveAs("Tiff", resultsDir + baseName + compositeTifSuffix);
    rename(dupTitle); // saveAs changes the window title to the saved filename; restore it

    // --- B2) save a flattened (RGB) version of the same composite ------------------
    run("Select None"); // never burn a stray selection into the flattened copy
    run("Flatten");
    saveAs("Tiff", resultsDir + baseName + compositeFlatTifSuffix);
    close(); // close flattened copy
    selectWindow(dupTitle);
    print("     Composite saved.");

    // --- 3) Segmentation: either run Cellpose, or load existing/manual ROIs -------
    roiManager("reset");
    roiZipPath       = resultsDir + baseName + roiZipSuffix;
    roiZipManualPath = resultsDir + baseName + roiZipManualSuffix;

    if (processingMode == "Segment") {
        selectWindow(dupTitle);
        titlesBeforeCP = getList("image.titles");
        // "nuc" model segments a SINGLE channel (DAPI); ch2=0 means "no second channel".
        cellposeArgs = "env_path=" + cellposeEnvPath + " env_type=" + cellposeEnvType +
                       " model=" + cellposeModel + " model_path=" + cellposeModelPath +
                       " diameter=" + cellposeDiameter +
                       " additional_flags=" + cellposeExtraArgs +
                       " ch1=" + dapiIdx + " ch2=0";
        run("Cellpose ...", cellposeArgs);
        titlesAfterCP = getList("image.titles");
        labelTitle = firstNewTitle(titlesAfterCP, titlesBeforeCP);
        if (labelTitle == "") {
            print("ERROR: Cellpose did not produce a new label image for " + imgTitle + " - skipping this image.");
            if (isOpen(dupTitle)) close(dupTitle);
            roiManager("reset");
            run("Clear Results");
            return;
        }

        // Check whether Cellpose actually found any nuclei BEFORE calling "Label
        // image to ROIs" - on an all-background (max=0) label image, that command
        // pops up a blocking "ROI Manager: The list is empty" dialog that would
        // otherwise stall/halt an unattended batch run.
        selectWindow(labelTitle);
        getStatistics(dummyArea, dummyMean, labelMin, labelMax);
        if (labelMax == 0) {
            print("WARNING: Cellpose found 0 nuclei for " + imgTitle + " - recording 0 cells and moving to next image.");
            if (isOpen(dupTitle)) close(dupTitle);
            if (isOpen(labelTitle)) close(labelTitle);
            roiManager("reset");
            run("Clear Results");
            appendZeroCellSummaryRow(lifBaseName, imgTitle, resultsDir);
            return;
        }

        // --- 4) label image -> ROIs ---------------------------------------------
        selectWindow(labelTitle);
        roiManager("reset");
        run("Label image to ROIs", "rm=[RoiManager[visible=true]]");

        // save the fresh ROI set so it can be manually corrected later in Update mode
        roiManager("Save", roiZipPath);
        print("     Cellpose segmentation done, ROIs saved to: " + roiZipPath);

    } else {
        // --- Update mode: load manually-corrected ROIs if present, else the original ---
        if (File.exists(roiZipManualPath)) {
            roiManager("Open", roiZipManualPath);
            print("     Update mode: loaded manually-corrected ROIs from " + roiZipManualPath);
        } else if (File.exists(roiZipPath)) {
            roiManager("Open", roiZipPath);
            print("     Update mode: no manual ROI file found, loaded original ROIs from " + roiZipPath);
        } else {
            print("ERROR: Update mode: no ROI file found for " + baseName +
                  " (looked for " + roiZipManualPath + " and " + roiZipPath + ") - skipping this image.");
            if (isOpen(dupTitle)) close(dupTitle);
            roiManager("reset");
            run("Clear Results");
            return;
        }

        nRoisLoaded = roiManager("count");
        if (nRoisLoaded == 0) {
            print("WARNING: Loaded ROI set has 0 ROIs for " + imgTitle + " - recording 0 cells and moving to next image.");
            if (isOpen(dupTitle)) close(dupTitle);
            roiManager("reset");
            run("Clear Results");
            appendZeroCellSummaryRow(lifBaseName, imgTitle, resultsDir);
            return;
        }

        // rebuild a label image from the loaded ROIs (label value = ROI index + 1, in the
        // same order ROIs are stored/measured) so the intensity color-coding step
        // downstream has a label image to work with, matching Cellpose's own output.
        selectWindow(dupTitle);
        width = getWidth(); height = getHeight();
        labelTitle = cleanBaseName + "_labels";
        if (isOpen(labelTitle)) close(labelTitle);
        newImage(labelTitle, "16-bit black", width, height, 1);
        for (r = 0; r < nRoisLoaded; r++) {
            roiManager("select", r);
            setColor(r + 1);
            fill();
        }
        run("Select None");
    }
    nCells = roiManager("count");

    // --- 5) measure Hes intensity inside every nuclear ROI --------------------------
    selectWindow(dupTitle);
    Stack.setChannel(hesIdx); // make sure the Hes channel is the active one being measured
    roiManager("Show All");
    run("Set Measurements...", measurementsString);
    roiManager("Measure");
    nCells = nResults;
    print("     Measured Hes intensity for " + nCells + " nucleus/nuclei.");

    // --- 5b) Hes background outside the nuclei (global + per-nucleus ring) ----------
    // Adds the Bg*/..._minus_Bg* columns to the Results table and returns the
    // per-image values [median, mean, nPixels, areaFraction, nRingsEmpty].
    bgRes = computeHesBackground(dupTitle, hesIdx, cleanBaseName);
    print("     Hes background (global, outside nuclei+" + bgGapPx + "px): median=" +
          d2s(bgRes[0], 2) + "  mean=" + d2s(bgRes[1], 2) +
          "  (" + d2s(100 * bgRes[3], 1) + "% of image)");
    if (bgRes[4] > 0) {
        print("     WARNING: " + bgRes[4] + " nucleus/nuclei had fewer than " + bgRingMinPixels +
              " usable ring pixels - their BgRing_* / *_minus_BgRing values are NaN.");
    }

    // --- A) save per-nucleus results -------------------------------------------------
    saveAs("Results", resultsDir + baseName + cellsCsvSuffix);
    print("     Saved: " + baseName + cellsCsvSuffix);

    // --- F) background QC image (after the CSV, so a QC problem never loses data) -------------------------------------------------------
    if (saveBgQcImage) {
        saveBackgroundQcImage(dupTitle, hesIdx, cleanBaseName, resultsDir + baseName + bgQcTifSuffix);
        print("     Saved: " + baseName + bgQcTifSuffix);
    }

    // --- append to run-wide summary table -----------------------------------------
    appendSummaryRow(lifBaseName, imgTitle, nCells, "Results", bgRes);
    finalizeSummary(resultsDir); // checkpoint: save Summary_Cells.csv after every image, so
                                  // results are never lost even if a later image/file fails

    // --- C) overlay (ROIs, no labels) flattened on composite -----------------------
    saveFlattenedWithOutlines(dupTitle, resultsDir + baseName + overlayTifSuffix);
    print("     Saved: " + baseName + overlayTifSuffix);

    // --- E) overlay (ROIs, no labels) flattened on the Hes channel ALONE ------------
    // Extra QC image: same ROI outlines, but on just the Hes signal (no DAPI channel
    // competing for attention), so segmentation quality against the actual Hes
    // staining can be checked directly.
    selectWindow(dupTitle);
    roiManager("Show None");
    roiManager("Deselect");
    run("Select None");             // no selection, or Duplicate would crop to it
    Stack.setChannel(hesIdx);
    hesOnlyTitle = cleanBaseName + "_hesOnly";
    if (isOpen(hesOnlyTitle)) close(hesOnlyTitle);
    run("Duplicate...", "title=" + hesOnlyTitle); // duplicates only the current channel/slice
    selectWindow(hesOnlyTitle);
    run(hesColor);
    setMinAndMax(hesMin, hesMax);
    saveFlattenedWithOutlines(hesOnlyTitle, resultsDir + baseName + hesOnlyOverlayTifSuffix);
    if (isOpen(hesOnlyTitle)) close(hesOnlyTitle);
    print("     Saved: " + baseName + hesOnlyOverlayTifSuffix);

    // --- D) Hes intensity color-coded label image -----------------------------------
    selectWindow(labelTitle);
    titlesBeforeHI = getList("image.titles");
    call("inra.ijpb.plugins.LabelToValuePlugin.process", "Table=Results",
         "Column=" + intensityColorColumn, "Min=" + hesColorMin, "Max=" + hesColorMax);
    titlesAfterHI = getList("image.titles");
    hiTitle = firstNewTitle(titlesAfterHI, titlesBeforeHI);
    if (hiTitle == "") {
        print("WARNING: " + intensityColorColumn + " color-coded image was not created for " +
              imgTitle + " - Cells.csv/Overlay were still saved, but this image is missing.");
    } else {
        selectWindow(hiTitle);
        setMinAndMax(hesColorMin, hesColorMax);
        run("Fire");
        run("Calibration Bar...", "location=[" + calBarLocation + "] fill=" + calBarFill +
            " label=" + calBarLabelColor + " number=" + calBarNumberLabels +
            " decimal=" + calBarDecimal + " font=" + calBarFontSize +
            " zoom=" + calBarZoom + " overlay");
        saveAs("Tiff", resultsDir + baseName + hesColorTifSuffix);
        rename(hiTitle);

        run("Select None"); // never burn a stray selection into the flattened copy
        run("Flatten");
        saveAs("Tiff", resultsDir + baseName + hesColorFlatTifSuffix);
        close(); // close flattened copy
        print("     Saved: " + baseName + hesColorTifSuffix + " / " + baseName + hesColorFlatTifSuffix);
    }

    // --- cleanup ------------------------------------------------------------------
    if (isOpen(dupTitle)) close(dupTitle);
    if (isOpen(labelTitle)) close(labelTitle);
    if (isOpen(hiTitle)) close(hiTitle);
    roiManager("reset");
    run("Clear Results");
    print("  <- Done: " + imgTitle);
}


// =============================================================================
//  6) SUMMARY TABLE (one row per image, whole run)
// =============================================================================

function appendSummaryRow(lifBaseName, imgTitle, nCells, tableName, bgRes) {
    Table.set("File", summaryRowIndex, lifBaseName, "Summary");
    Table.set("Image", summaryRowIndex, imgTitle, "Summary");
    Table.set("N_Cells", summaryRowIndex, nCells, "Summary");
    // per-image background values (from computeHesBackground)
    Table.set("BgGlobal_Median", summaryRowIndex, bgRes[0], "Summary");
    Table.set("BgGlobal_Mean", summaryRowIndex, bgRes[1], "Summary");
    Table.set("BgGlobal_NPixels", summaryRowIndex, bgRes[2], "Summary");
    Table.set("BgGlobal_AreaFraction", summaryRowIndex, bgRes[3], "Summary");
    Table.set("N_RingsEmpty", summaryRowIndex, bgRes[4], "Summary");

    // Build the list of columns actually present in this table, so we can check
    // membership before calling Table.getColumn (avoids relying on try/catch,
    // which is not supported by every Fiji/ImageJ macro interpreter version).
    headingsStr = Table.headings(tableName);
    rawCols = split(headingsStr, "\t");
    availableCols = newArray(rawCols.length);
    for (j = 0; j < rawCols.length; j++) {
        availableCols[j] = String.trim(rawCols[j]);
    }

    for (i = 0; i < summaryColumns.length; i++) {
        colName = summaryColumns[i];
        if (arrayContains(availableCols, colName)) {
            // NOTE: an array returned by a built-in (Table.getColumn) must be stored in
            // a variable first - passing it straight into a user-defined function
            // fails with "Array or string expected".
            rawValues = Table.getColumn(colName, tableName);
            values = removeNaN(rawValues);
            if (values.length > 0) {
                Array.getStatistics(values, minV, maxV, meanV, stdV);
            } else {
                minV = NaN; maxV = NaN; meanV = NaN;
            }
            Table.set(colName + "_Mean", summaryRowIndex, meanV, "Summary");
            Table.set(colName + "_Min", summaryRowIndex, minV, "Summary");
            Table.set(colName + "_Max", summaryRowIndex, maxV, "Summary");
        } else {
            print("WARNING: summary column '" + colName + "' not found in " + tableName + " - skipping.");
        }
    }
    summaryRowIndex++;
}

function finalizeSummary(resultsDir) {
    Table.update("Summary");
    Table.save(resultsDir + summaryCsvName, "Summary");
    print("Summary saved to: " + resultsDir + summaryCsvName);
}

// Lightweight summary row for images where 0 cells were detected/loaded - just
// records File/Image/N_Cells=0 and leaves the measurement columns blank for
// that row, then checkpoints the summary file immediately.
function appendZeroCellSummaryRow(lifBaseName, imgTitle, resultsDir) {
    Table.set("File", summaryRowIndex, lifBaseName, "Summary");
    Table.set("Image", summaryRowIndex, imgTitle, "Summary");
    Table.set("N_Cells", summaryRowIndex, 0, "Summary");
    summaryRowIndex++;
    finalizeSummary(resultsDir);
}


// =============================================================================
//  6b) HES BACKGROUND (outside the nuclei)
// =============================================================================

// Computes the Hes background for the current image from the ROIs in the ROI
// Manager (same order as the rows of the "Results" table), writes the per-nucleus
// background + corrected columns into "Results", and returns
//   newArray(globalMedian, globalMean, globalNPixels, globalAreaFraction, nRingsEmpty)
//
// How it works: the Hes channel is copied to a 32-bit image, and every pixel
// inside (nuclei grown by bgGapPx) is set to NaN. ImageJ ignores NaN pixels when
// measuring 32-bit images, so:
//   - measuring the whole image gives the GLOBAL background, and
//   - measuring each nucleus grown by (bgGapPx + bgRingWidthPx) gives that
//     nucleus' RING (its own nucleus and all neighbours are already NaN).
function computeHesBackground(dupTitle, hesIdx, cleanBaseName) {
    result = newArray(NaN, NaN, 0, NaN, 0);
    nRois = roiManager("count");
    if (nRois == 0) return result;

    // --- 32-bit copy of the Hes channel ----------------------------------------
    selectWindow(dupTitle);
    roiManager("Show None");
    run("Select None");
    getPixelSize(pxUnit, pxW, pxH);
    pxArea = pxW * pxH;
    bgTitle = cleanBaseName + "_hesBg";
    if (isOpen(bgTitle)) close(bgTitle);
    run("Duplicate...", "title=" + bgTitle + " duplicate channels=" + hesIdx);
    selectWindow(bgTitle);
    run("32-bit");
    imgPixels = getWidth() * getHeight();

    // --- blank out nuclei + gap ----------------------------------------------------
    selectNucleiExclusionZone();
    run("Set...", "value=NaN");
    run("Select None");
    roiManager("Deselect");

    // --- GLOBAL background -------------------------------------------------------------
    List.setMeasurements;
    globalNPix = round(List.getValue("Area") / pxArea);
    if (globalNPix > 0) {
        globalMedian = List.getValue("Median");
        globalMean   = List.getValue("Mean");
    } else {
        globalMedian = NaN;
        globalMean   = NaN;
        print("WARNING: no pixels left outside the nuclei (+" + bgGapPx +
              "px) - global background is NaN for this image.");
    }
    List.clear();

    // --- per-nucleus RING background + corrected columns ------------------------------
    nRows = nResults;
    if (nRows != nRois) {
        print("WARNING: Results has " + nRows + " rows but ROI Manager has " + nRois +
              " ROIs - background columns are filled for the first " + minOf(nRows, nRois) + " only.");
    }
    nLoop = minOf(nRows, nRois);
    ringOuterPx = bgGapPx + bgRingWidthPx;
    nRingsEmpty = 0;
    for (r = 0; r < nLoop; r++) {
        selectWindow(bgTitle);
        roiManager("select", r);
        run("Enlarge...", "enlarge=" + ringOuterPx + " pixel");
        List.setMeasurements;
        ringNPix = round(List.getValue("Area") / pxArea);
        if (ringNPix >= bgRingMinPixels) {
            ringMedian = List.getValue("Median");
            ringMean   = List.getValue("Mean");
        } else {
            ringMedian = NaN;
            ringMean   = NaN;
            nRingsEmpty++;
        }
        List.clear();

        nucMean = getResult("Mean", r);
        nucRaw  = getResult("RawIntDen", r);
        nucNPix = round(getResult("Area", r) / pxArea);

        setResult("BgGlobal_Median", r, globalMedian);
        setResult("BgGlobal_Mean", r, globalMean);
        setResult("Mean_minus_BgGlobal", r, nucMean - globalMedian);
        setResult("RawIntDen_minus_BgGlobal", r, nucRaw - globalMedian * nucNPix);
        setResult("BgRing_Median", r, ringMedian);
        setResult("BgRing_Mean", r, ringMean);
        setResult("BgRing_NPixels", r, ringNPix);
        setResult("Mean_minus_BgRing", r, nucMean - ringMedian);
        setResult("RawIntDen_minus_BgRing", r, nucRaw - ringMedian * nucNPix);
    }
    updateResults();

    selectWindow(bgTitle);
    run("Select None");
    roiManager("Deselect");
    close(bgTitle);

    result[0] = globalMedian;
    result[1] = globalMean;
    result[2] = globalNPix;
    result[3] = globalNPix / imgPixels;
    result[4] = nRingsEmpty;
    return result;
}

// Makes a selection, on the CURRENT image, covering every ROI in the ROI Manager
// grown by bgGapPx pixels (= the zone excluded from all background estimates).
function selectNucleiExclusionZone() {
    n = roiManager("count");
    if (n == 1) {
        roiManager("select", 0);
    } else {
        allIdx = Array.getSequence(n);
        roiManager("select", allIdx);
        roiManager("Combine");
    }
    if (bgGapPx > 0) run("Enlarge...", "enlarge=" + bgGapPx + " pixel");
}

// Saves a flattened QC image of the Hes channel showing the excluded zone
// (shaded) and the outer edge of every ring region (outlines).
function saveBackgroundQcImage(dupTitle, hesIdx, cleanBaseName, outPath) {
    selectWindow(dupTitle);
    roiManager("Show None");
    run("Select None");
    qcTitle = cleanBaseName + "_bgQC";
    if (isOpen(qcTitle)) close(qcTitle);
    run("Duplicate...", "title=" + qcTitle + " duplicate channels=" + hesIdx);
    selectWindow(qcTitle);
    run(hesColor);
    setMinAndMax(hesMin, hesMax);

    selectNucleiExclusionZone();
    Overlay.addSelection("", 0, bgQcExcludedFill);
    run("Select None");
    roiManager("Deselect");

    // Draw the ring pixels that are ACTUALLY measured: the union of all rings
    // with the exclusion zone (every nucleus + gap) cut out - exactly the pixels
    // that are not NaN inside a ring in computeHesBackground(). Drawing each ring's
    // outer outline instead would cross neighbouring nuclei and wrongly suggest
    // they are part of the ring.
    w = getWidth();
    h = getHeight();
    ringMaskTitle = cleanBaseName + "_ringMask";
    if (isOpen(ringMaskTitle)) close(ringMaskTitle);
    newImage(ringMaskTitle, "8-bit black", w, h, 1);
    selectWindow(ringMaskTitle);
    ringOuterPx = bgGapPx + bgRingWidthPx;
    setColor(255);
    n = roiManager("count");
    for (r = 0; r < n; r++) {
        roiManager("select", r);
        run("Enlarge...", "enlarge=" + ringOuterPx + " pixel");
        fill();
    }
    selectNucleiExclusionZone();
    setColor(0);
    fill();
    run("Select None");
    roiManager("Deselect");
    setThreshold(255, 255);
    run("Create Selection");
    hasRingPixels = (selectionType() != -1);
    // Hand the ring selection over to the QC image through the ROI Manager
    // (temporarily added as the LAST entry, removed again right after).
    // "Restore Selection" is not reliable for this across images.
    tmpIdx = -1;
    if (hasRingPixels) {
        roiManager("add");
        tmpIdx = roiManager("count") - 1;
    }
    run("Select None");
    resetThreshold();

    selectWindow(qcTitle);
    if (tmpIdx >= 0) {
        roiManager("select", tmpIdx);
        if (selectionType() != -1) {
            Overlay.addSelection(bgQcRingColor, 1);
        } else {
            print("     WARNING: could not draw the ring outlines on the BgQC image (QC image only - measurements unaffected).");
        }
        roiManager("select", tmpIdx);
        roiManager("delete");   // keep the ROI Manager = nuclei only
        run("Select None");
        roiManager("Deselect");
    } else {
        print("     Note: no ring pixels to draw on the BgQC image.");
    }
    if (isOpen(ringMaskTitle)) close(ringMaskTitle);
    selectWindow(qcTitle);

    run("Select None"); // never burn a stray selection into the flattened copy
    run("Flatten");
    saveAs("Tiff", outPath);
    close(); // close flattened copy
    if (isOpen(qcTitle)) close(qcTitle);
}

// Returns a copy of 'arr' without NaN entries.
function removeNaN(arr) {
    clean = newArray(0);
    for (k = 0; k < arr.length; k++) {
        if (!isNaN(arr[k])) clean = Array.concat(clean, arr[k]);
    }
    return clean;
}


// Saves a flattened RGB copy of image 'srcTitle' with the outlines of ALL ROIs in the
// ROI Manager drawn in yellow (width outlineWidth), no labels, plus a scale bar. The outlines are drawn straight into the
// pixels of the flattened copy from the stored ROI coordinates (roiManager "Draw"),
// instead of flattening the on-screen "Show All" display. Flattening the display
// could also burn in the currently selected ROI, partly drawn at the window's zoom
// level, which showed up as stray triangles with one corner at half the image
// coordinates of a real cell.
function saveFlattenedWithOutlines(srcTitle, outPath) {
    flatID = makeFlatCopyWithScaleBar(srcTitle);
    run("Line Width...", "line=" + outlineWidth);
    setForegroundColor(255, 255, 0);
    roiManager("Deselect");         // nothing selected -> "Draw" draws every ROI
    roiManager("Draw");
    saveAndCloseFlatCopy(flatID, outPath, srcTitle);
}

// Makes a plain RGB copy of image 'srcTitle' (no overlay, no selection) and burns a
// scale bar into it. Returns the copy's image ID; the copy is the active image.
// The source must NOT carry an overlay when it is flattened: for multi-channel
// images Fiji's Flatten can then convert the source in place instead of making a
// copy, after which the source window is gone.
function makeFlatCopyWithScaleBar(srcTitle) {
    selectWindow(srcTitle);
    srcID = getImageID();
    getPixelSize(fcUnit, fcPw, fcPh);
    roiManager("Show None");
    roiManager("Deselect");
    run("Select None");
    run("Remove Overlay");
    run("Flatten");                 // plain RGB copy of what is displayed
    flatID = getImageID();
    if (flatID == srcID) {
        exit("Flatten did not create a new image for '" + srcTitle + "'.\n" +
             "Please send this message and the Log window to the macro maintainer.");
    }
    setVoxelSize(fcPw, fcPh, 1, fcUnit); // make sure the copy carries the calibration
    addScaleBarBurned();
    selectImage(flatID);
    run("Select None");
    return flatID;
}

// Saves the RGB copy 'flatID' to outPath, closes it, and goes back to the source.
function saveAndCloseFlatCopy(flatID, outPath, srcTitle) {
    run("Line Width...", "line=1"); // restore the default line width
    selectImage(flatID);
    run("Select None");
    saveAs("Tiff", outPath);
    selectImage(flatID);
    close();
    if (isOpen(srcTitle)) {
        selectWindow(srcTitle);
        run("Select None");
    }
}

// Burns a scale bar into the current (RGB) image, using the scaleBar* settings.
// Skipped with a warning if the image is not spatially calibrated.
function addScaleBarBurned() {
    getPixelSize(sbUnit, sbPw, sbPh);
    if (sbUnit == "pixel" || sbUnit == "pixels" || sbUnit == "") {
        print("     WARNING: image has no spatial calibration - no scale bar added.");
        return;
    }
    sbOpts = "width=" + scaleBarLength + " height=" + scaleBarThickness +
             " thickness=" + scaleBarThickness + " font=" + scaleBarFontSize +
             " color=" + scaleBarColor + " background=None location=[" + scaleBarLocation + "] horizontal bold";
    if (scaleBarFontSize <= 0) sbOpts = sbOpts + " hide";
    run("Scale Bar...", sbOpts);
}

// =============================================================================
//  6c) COMBINED PER-CELL TABLES (all images in the Results folder)
// =============================================================================

// Writes <combinedTableName>.csv (all columns) and <combinedTableName>_slim.csv
// (combinedSlimCols only) into resultsDir, combining EVERY <image>_Cells.csv in that
// folder (not only the images of this run, so the tables stay complete after partial
// re-runs). Two columns are added in front of each row: SourceFile (the image prefix)
// and ResultsFolder (the Results folder name). Columns are matched by name; a column
// missing from a file is written as NaN.
function writeCombinedCellsTables(resultsDir) {
    if (!saveCombinedTables) return;
    folderPath = resultsDir;
    if (endsWith(folderPath, File.separator) || endsWith(folderPath, "/")) folderPath = substring(folderPath, 0, lengthOf(folderPath) - 1);
    folderName = File.getName(folderPath);

    list = getFileList(resultsDir);
    files = newArray(0);
    for (i = 0; i < list.length; i++) {
        nm = list[i];
        if (endsWith(nm, cellsCsvSuffix) && !startsWith(nm, "Summary") && !startsWith(nm, "All_Cells")) {
            files = Array.concat(files, nm);
        }
    }
    if (files.length == 0) {
        print("Combined tables: no *" + cellsCsvSuffix + " files in " + resultsDir + " - skipped.");
        return;
    }
    Array.sort(files);

    firstLines = combReadLines(resultsDir + files[0]);
    if (firstLines.length < 1) {
        print("Combined tables: " + files[0] + " is empty - skipped.");
        return;
    }
    rawHeader = split(firstLines[0], ",");
    refCols = combCleanHeader(rawHeader);
    slimIdx = newArray(0);
    for (k = 0; k < combinedSlimCols.length; k++) {
        p = combIndexOf(refCols, combinedSlimCols[k]);
        if (p >= 0) slimIdx = Array.concat(slimIdx, p);
    }

    outFull = resultsDir + combinedTableName + ".csv";
    outSlim = resultsDir + combinedTableName + "_slim.csv";
    totalRows = 0;
    // ImageJ macros can only have ONE output file open at a time -> two passes
    for (pass = 0; pass < 2; pass++) {
        if (pass == 0) f = File.open(outFull);
        else f = File.open(outSlim);
        head = "SourceFile,ResultsFolder";
        if (pass == 0) {
            for (k = 0; k < refCols.length; k++) head = head + "," + refCols[k];
        } else {
            for (k = 0; k < slimIdx.length; k++) head = head + "," + refCols[slimIdx[k]];
        }
        print(f, head);
        for (i = 0; i < files.length; i++) {
            srcName = substring(files[i], 0, lengthOf(files[i]) - lengthOf(cellsCsvSuffix));
            lines = combReadLines(resultsDir + files[i]);
            if (lines.length < 2) continue;
            rawHeader = split(lines[0], ",");
            cols = combCleanHeader(rawHeader);
            map = newArray(refCols.length);
            for (k = 0; k < refCols.length; k++) map[k] = combIndexOf(cols, refCols[k]);
            for (j = 1; j < lines.length; j++) {
                line = lines[j];
                if (lengthOf(line) == 0) continue;
                vals = split(line, ",");
                outLine = srcName + "," + folderName;
                if (pass == 0) {
                    for (k = 0; k < refCols.length; k++) {
                        if (map[k] >= 0 && map[k] < vals.length) outLine = outLine + "," + vals[map[k]];
                        else outLine = outLine + ",NaN";
                    }
                } else {
                    for (k = 0; k < slimIdx.length; k++) {
                        q = map[slimIdx[k]];
                        if (q >= 0 && q < vals.length) outLine = outLine + "," + vals[q];
                        else outLine = outLine + ",NaN";
                    }
                }
                print(f, outLine);
                if (pass == 0) totalRows++;
            }
        }
        File.close(f);
    }
    print("Combined tables: " + totalRows + " rows from " + files.length + " image(s) -> " +
          combinedTableName + ".csv / " + combinedTableName + "_slim.csv");
}

// Reads a text file and returns its lines without Windows "\r" characters.
function combReadLines(path) {
    txt = File.openAsString(path);
    if (lengthOf(txt) == 0) {
        emptyLines = newArray(0);
        return emptyLines;
    }
    txt = replace(txt, "\r", "");
    lines = split(txt, "\n");
    return lines;
}

// Trims header names; the unnamed first column (row number) becomes RowInFile.
function combCleanHeader(cols) {
    for (k = 0; k < cols.length; k++) {
        c = String.trim(cols[k]);
        if (c == "") c = "RowInFile";
        cols[k] = c;
    }
    return cols;
}

// Position of 'name' in 'arr', or -1.
function combIndexOf(arr, name) {
    for (k = 0; k < arr.length; k++) {
        if (arr[k] == name) return k;
    }
    return -1;
}

// Turns the dialog text (e.g. ".tif, .tiff" or "lif") into an array of lower-case
// extensions with a leading dot. Falls back to ".lif" if nothing usable is given.
function parseExtensions(txt) {
    exts = newArray(0);
    txt = "" + txt;
    if (txt == "0") txt = "";      // an empty text field can arrive as "0"
    txt = replace(txt, ";", ",");
    txt = replace(txt, " ", ",");
    rawParts = split(txt, ",");
    for (k = 0; k < rawParts.length; k++) {
        e = toLowerCase(String.trim(rawParts[k]));
        if (e == "") continue;
        if (!startsWith(e, ".")) e = "." + e;
        exts = Array.concat(exts, e);
    }
    if (exts.length == 0) exts = newArray(".lif");
    return exts;
}

// =============================================================================
//  7) PARAMETERS FILE
// =============================================================================

function writeParametersFile(resultsDir) {
    getDateAndTime(y, mo, dow, d, h, mi, s, ms);
    timestamp = "" + y + "-" + pad2(mo + 1) + "-" + pad2(d) + " " + pad2(h) + ":" + pad2(mi) + ":" + pad2(s);

    f = File.open(resultsDir + parametersTxtName);
    print(f, scriptName + " - run parameters");
    print(f, "Timestamp: " + timestamp);
    print(f, "Run mode (Segment/Update): " + processingMode);
    print(f, "Image selection mode: " + processMode);
    print(f, "Image file chosen in dialog: " + file_name);
    print(f, "-----------------------------------------------------------");
    print(f, "resultsSubfolderName = " + resultsSubfolderName);
    print(f, "fileSuffix = " + fileSuffix);
    print(f, "fileExtensions = " + String.join(fileExtensions, ", "));
    print(f, "saveCombinedTables = " + saveCombinedTables + " (" + combinedTableName + ".csv / _slim.csv)");
    print(f, "");
    print(f, "CH_HES = " + CH_HES);
    print(f, "CH_DAPI = " + CH_DAPI);
    print(f, "");
    print(f, "hesColor = " + hesColor);
    print(f, "hesMin = " + hesMin);
    print(f, "hesMax = " + hesMax);
    print(f, "dapiColor = " + dapiColor);
    print(f, "dapiMin = " + dapiMin);
    print(f, "dapiMax = " + dapiMax);
    print(f, "");
    print(f, "cellposeEnvPath = " + cellposeEnvPath);
    print(f, "cellposeEnvType = " + cellposeEnvType);
    print(f, "cellposeModel = " + cellposeModel);
    print(f, "cellposeModelPath = " + cellposeModelPath);
    print(f, "cellposeDiameter = " + cellposeDiameter);
    print(f, "cellposeExtraArgs = " + cellposeExtraArgs);
    print(f, "");
    print(f, "measurementsString = " + measurementsString);
    print(f, "");
    print(f, "Hes background: GLOBAL = median (and mean) of all Hes pixels outside the nuclei grown by bgGapPx;");
    print(f, "                RING   = median (and mean) per nucleus, from bgGapPx to bgGapPx+bgRingWidthPx");
    print(f, "                         outside its edge, all nuclei (+gap) excluded.");
    print(f, "                Corrections subtract the MEDIAN; negative values are kept.");
    print(f, "bgGapPx = " + bgGapPx);
    print(f, "bgRingWidthPx = " + bgRingWidthPx);
    print(f, "bgRingMinPixels = " + bgRingMinPixels);
    print(f, "saveBgQcImage = " + saveBgQcImage);
    print(f, "");
    summaryColumnsStr = String.join(summaryColumns, ", ");
    print(f, "summaryColumns = " + summaryColumnsStr);
    print(f, "");
    print(f, "intensityColorColumn = " + intensityColorColumn);
    print(f, "hesColorMin = " + hesColorMin);
    print(f, "hesColorMax = " + hesColorMax);
    print(f, "calBarLocation = " + calBarLocation);
    print(f, "calBarFill = " + calBarFill);
    print(f, "calBarLabelColor = " + calBarLabelColor);
    print(f, "calBarNumberLabels = " + calBarNumberLabels);
    print(f, "calBarDecimal = " + calBarDecimal);
    print(f, "calBarFontSize = " + calBarFontSize);
    print(f, "calBarZoom = " + calBarZoom);
    print(f, "outlineWidth = " + outlineWidth);
    print(f, "scaleBar = " + scaleBarLength + " (calibrated units), thickness " + scaleBarThickness +
              " px, font " + scaleBarFontSize + ", " + scaleBarColor + ", " + scaleBarLocation);
    print(f, "");
    print(f, "Output suffixes: " + cellsCsvSuffix + ", " + compositeTifSuffix + ", " +
              compositeFlatTifSuffix + ", " + overlayTifSuffix + ", " + hesOnlyOverlayTifSuffix +
              ", " + hesColorTifSuffix + ", " + hesColorFlatTifSuffix + ", " + bgQcTifSuffix);
    print(f, "roiZipSuffix = " + roiZipSuffix);
    print(f, "roiZipManualSuffix = " + roiZipManualSuffix);
    print(f, "Summary file: " + summaryCsvName);
    File.close(f);

    print("Parameters saved to: " + resultsDir + parametersTxtName);
}

function pad2(n) {
    if (n < 10) {
        padded = "0" + n;
        return padded;
    }
    padded = "" + n;
    return padded;
}


// =============================================================================
//  8) HELPER FUNCTIONS
// =============================================================================

function ensureResultsDir(dir) {
    rd = dir + resultsSubfolderName + File.separator;
    if (!File.isDirectory(rd)) {
        print("Results folder does not exist yet - creating: " + rd);
        File.makeDirectory(rd);
        if (!File.isDirectory(rd)) {
            exit("Could not create the Results folder: " + rd +
                 " (check that the parent folder is writable).");
        }
    } else {
        print("Using existing Results folder: " + rd);
    }
    return rd;
}

// Builds a filesystem-safe "<lifBaseName>_<series>" prefix from a Bio-Formats
// window title such as "48h_notch_control_cd31AF555_Hes1AF647.lif - Image1-1".
// Falls back to sanitizing the raw title if no " - " separator is found.
function getBaseName(imgTitle, lifBaseName) {
    // drop a trailing file extension (e.g. an opened "<lif>.lif - Image 2.tif")
    t0 = imgTitle;
    lt = toLowerCase(t0);
    knownExt = newArray(".lif", ".tiff", ".tif");
    for (k = 0; k < knownExt.length; k++) {
        if (endsWith(lt, knownExt[k])) {
            t0 = substring(t0, 0, lengthOf(t0) - lengthOf(knownExt[k]));
            k = knownExt.length;
        }
    }
    // a lif base name taken from such a TIF's file name still ends in ".lif - <series>"
    lifSep = indexOf(lifBaseName, ".lif - ");
    if (lifSep > 0) lifBaseName = substring(lifBaseName, 0, lifSep);
    sepIndex = lastIndexOf(t0, " - ");
    if (lifBaseName != "" && sepIndex >= 0) {
        seriesPart = substring(t0, sepIndex + 3);
        base = lifBaseName + "_" + seriesPart;
    } else {
        base = t0;
    }
    cleanBase = sanitizeName(base);
    return cleanBase;
}

function sanitizeName(s) {
    cleaned = replace(s, "[^A-Za-z0-9_-]+", "_");
    return cleaned;
}

// Returns elements present in 'after' but not in 'before' (order preserved).
function arrayDiff(after, before) {
    result = newArray(0);
    for (i = 0; i < after.length; i++) {
        if (!arrayContains(before, after[i])) {
            result = Array.concat(result, after[i]);
        }
    }
    return result;
}

// Returns the first title in 'after' that is not present in 'before', or "".
function firstNewTitle(after, before) {
    for (i = 0; i < after.length; i++) {
        if (!arrayContains(before, after[i])) {
            newTitle = after[i];
            return newTitle;
        }
    }
    return "";
}

function arrayContains(arr, val) {
    for (i = 0; i < arr.length; i++) {
        if (arr[i] == val) return true;
    }
    return false;
}

// @File dialog fields that are left blank can be harvested as the literal
// string "null" (or occasionally "None") rather than an empty string.
function isUnset(v) {
    trimmed = String.trim(v);
    lower = toLowerCase(trimmed);
    unset = (trimmed == "" || lower == "null" || lower == "none");
    return unset;
}
