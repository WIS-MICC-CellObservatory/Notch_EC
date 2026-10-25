// @String(label="Run mode", choices=("Segment", "Update"), style="list", value="Segment") runMode
// @String(label="Image selection", choices=("Single open image", "Single file (.lif or .tif)", "Whole folder"), style="list") processMode
// @String(label="File extension(s), comma-separated (for Single file / Whole folder), e.g. .lif or .tif, .tiff", value=".lif") fileTypeChoice
// @File(label="Image file, .lif or .tif (ignored for 'Single open image'; for 'Whole folder', pick any file inside the target folder)", style="open", required=false) file_name
// @File(label="Cellpose conda env folder (required for Segment mode)", style="directory", required=false, value="C:\\ProgramData\\anaconda3\\envs\\cellpose4") cellposeEnvPath
// @File(label="Cellpose model_path folder (required for Segment mode)", style="directory", required=false, value="C:\\Fiji\\fiji-latest-win64-jdk\\Fiji") cellposeModelPath
// @String(label="Cellpose model", choices=("cpsam", "cpsam_v2", "cpdino"), style="list", value="cpsam") cellposeModelChoice
// @String(label="Custom Cellpose model name (optional - overrides the list above if filled)", required=false, value="") cellposeModelCustom
// @String(label="Results subfolder name (leave blank for automatic 'Results_<model>')", required=false, value="Results") resultsSubfolderName
//
// =============================================================================
//  Workflow 1 : CDH5-AF647 images  -  Cell segmentation + elongation (AR) analysis
// =============================================================================
//
//  WHAT THIS SCRIPT DOES (per image):
//    1. Duplicates the membrane + nuclei channels from a 4-channel .lif or .tif image
//    2. Builds a 2-channel composite with fixed, comparable display ranges
//    3. Segmentation - one of two modes (chosen from a dialog at start):
//         - Segment: runs Cellpose (model chosen in the dialog - e.g. cpsam,
//           cpsam_v2, cpdino - using both channels) and converts
//           the resulting label image to ROIs, then SAVES those ROIs as
//           <prefix>_RoiSet.zip so they can be manually corrected later.
//         - Update: skips Cellpose entirely. Looks for
//           <prefix>_RoiSet_manual.zip first (a manually-corrected copy you
//           made yourself); if it doesn't exist, falls back to the original
//           <prefix>_RoiSet.zip. A label image is then reconstructed from
//           those ROIs so the AR color-coding step still works.
//           MANUAL CORRECTION STEP (done by you, in Fiji, between runs):
//             a) Run in Segment mode once to produce <prefix>_RoiSet.zip
//             b) Open that zip into the ROI Manager (More >> Open...)
//             c) Correct/add/delete ROIs as needed
//             d) ROI Manager > More >> Save... and save it as EXACTLY
//                <prefix>_RoiSet_manual.zip in the same Results folder
//             e) Re-run this script in Update mode - it will pick up your
//                corrected ROIs automatically for that image.
//    4. Measures shape descriptors for every cell using MorphoLibJ's "Analyze
//       Regions", run directly on the label image (NOT ROI Manager "Measure" -
//       the ROI-polygon Fit-Ellipse AR was found to be unreliable for Cellpose's
//       irregular/concave shapes; region-based Ellipse.Elong tracks visual
//       elongation correctly instead).
//    5. Saves:
//         A) <prefix>_Cells.csv        - per-cell measurements (Analyze Regions table)
//         B) <prefix>_Composite.tif    - membrane+nuclei composite (fixed LUT/range)
//            <prefix>_Composite_flat.tif - flattened (RGB) version of the same
//         C) <prefix>_Overlay.tif      - ROIs (no labels) flattened on composite
//         C2) Nuclei-only segments (Cellpose sometimes segments only the nucleus
//             of a cell): segments with Area <= nucleiOnlyMaxArea are
//               - drawn in a second colour in <prefix>_Overlay_marked.tif
//                 (saveMarkedOverlay), and/or
//               - left out of <prefix>_Overlay_filtered.tif and
//                 <prefix>_Cells_filtered.csv, and shown as background in
//                 <prefix>_AR_filtered.tif / _AR_filtered_flat.tif (saveFilteredOutputs).
//             <prefix>_Cells.csv and <prefix>_Overlay.tif always keep ALL segments.
//         D) <prefix>_AR.tif + <prefix>_AR_flat.tif - Ellipse.Elong color-coded
//              label image (with Fire LUT + calibration bar), flattened version included
//         (Segment mode only) <prefix>_RoiSet.zip - the ROI set itself
//    6. Appends one row (mean/min/max of every Analyze Regions descriptor) to a
//       single run-wide summary table
//
//  At the end of the run, the script also writes:
//    - Summary_Cells.csv          (one row per image, whole run)
//    - All_Cells_Shape.csv + All_Cells_Shape_slim.csv  (every <image>_Cells.csv in the
//      Results folder combined, with SourceFile + ResultsFolder columns;
//      ready for the Python scripts in analysis/)
//    - Workflow1_CDH5AF647_Parameters.txt  (all parameters used + timestamp)
//
//  RUN MODES (all chosen from ONE dialog when you click Run):
//    A) Single already-open image   (leave "Image file" blank)
//    B) Single file (.lif or .tif)   (pick the file in the dialog; all series
//                                     inside it are opened and processed; a .tif
//                                     is one image, named after its file)
//    C) Whole folder                 (pick ANY *CDH5AF647*.lif or .tif file that lives
//                                     in the target folder - the script then
//                                     processes every matching file in that
//                                     file's parent folder)
//
//  The SAME dialog also lets you choose Segment vs Update mode, (for
//  Segment mode) the Cellpose conda-env and model_path folders and the
//  Cellpose model, and the name of the Results subfolder, so nothing needs to
//  be edited in the code just to run the script day to day. Fields not needed
//  for your choice (e.g. the image file in "Single open image" mode, or the
//  Cellpose folders/model in "Update" mode) can be left blank / ignored.
//
//  CELLPOSE MODEL:
//    Pick cpsam, cpsam_v2 or cpdino from the list. To use any other model
//    (e.g. a custom-trained one stored in the model_path folder), type its
//    name in "Custom Cellpose model name" - when filled, it overrides the list.
//    Whatever name is used is passed as-is to the BIOP Cellpose wrapper
//    ("model=..."), so it must be a name your Cellpose installation can
//    resolve (built-in model, or a model file inside model_path).
//
//  RESULTS LOCATION:
//    A subfolder named from the dialog ("Results subfolder name") is created
//    next to the source file(s) (mode A falls back to a folder you choose if
//    the open image has no file on disk). If left blank, it defaults to
//    "Results_<model>" (e.g. Results_cpsam_v2), so runs with different models
//    never overwrite each other. In Update mode, give the SAME subfolder name
//    as the Segment run whose ROIs you want to reuse/correct.
//
//  REQUIRED DEPENDENCIES (enable via Help > Update... > Manage update sites):
//    - BIOP update site        -> Cellpose wrapper ("Cellpose ...") and
//                                  "Label image to ROIs"
//    - IJPB-plugins (MorphoLibJ) update site -> LabelToValuePlugin
//                                  ("Assign Measure to Label")
//    - Bio-Formats              -> bundled with Fiji, used to read .lif files (.tif files open directly)
//    - A local Cellpose Python/conda environment with cellpose (v4) and the
//      model(s) you want to use (cpsam / cpsam_v2 / cpdino / custom) installed and reachable at CELLPOSE_ENV_PATH below. Test it once
//      via Plugins > BIOP > Cellpose before running this script.
//
//  IMPLEMENTATION NOTES:
//    - Interactive "Color Balance..." / "Enhance Contrast" were replaced with
//      explicit setMinAndMax() calls driven by variables below, so that every
//      image in the dataset gets the SAME display range and can be compared
//      side by side.
//    - The duplicated "call(...LabelToValuePlugin...)" + "run(Assign Measure
//      to Label)" pair is redundant (both do the same
//      thing); only the direct call() is kept here for a fully scriptable,
//      dialog-free run.
//    - Cellpose IS called with a diameter parameter (cellposeDiameter, default 40.0):
//      cpsam's underlying code divides by this value, so omitting it / setting it to
//      0 crashes with a ZeroDivisionError. It does not need to be tuned per-image,
//      but the key must be present.
//
// =============================================================================


// =============================================================================
//  1) USER-EDITABLE PARAMETERS
// =============================================================================

scriptName            = "Workflow1_CDH5AF647";   // used for the parameters file name
// resultsSubfolderName is set from the dialog (see the @String parameters at the
// very top of this file), not here. If left blank it becomes "Results_<model>".
fileSuffix            = "CDH5AF647";             // used to find files in folder mode (must appear in filename)
// Input file types come from the "File extension(s)" field in the dialog
// (fileTypeChoice), e.g. ".lif", ".tif, .tiff" or ".lif, .tif". Case-insensitive;
// a missing leading dot is added. A .lif is processed series by series; any other
// file (e.g. .tif) is opened as one image.
fileExtensions = parseExtensions(fileTypeChoice);

// Combined per-cell tables written into the Results folder at the end of every run
// (all <image>_Cells.csv files in that folder, input for the Python scripts in analysis/)
saveCombinedTables    = true;
combinedTableName     = "All_Cells_Shape";
combinedSlimCols      = newArray("Label", "Area", "Perimeter", "Circularity",
    "Ellipse.Radius1", "Ellipse.Radius2", "Ellipse.Orientation", "Ellipse.Elong",
    "Convexity", "MaxFeretDiam", "OBox.Elongation", "GeodesicElongation",
    "AverageThickness", "Centroid.X", "Centroid.Y");

// ---- Channel roles (1-based, as in the ORIGINAL 4-channel image) -----------
CH_MEMBRANE = 3;
CH_NUCLEI   = 4;

// ---- Display settings for the 2-channel composite ---------------------------
membraneColor = "Magenta";
membraneMin   = 0;
// Membrane display MAX per time point (the signal gets brighter over time, so one
// fixed max saturates the later time points). The time point is recognised by
// these labels in the .lif / image name (capitals ignored); the first match wins.
// Images whose name contains none of the labels use membraneMax (default below).
// Display only: pixel values, segmentation and measurements are not affected.
timePointLabels       = newArray("12h", "48h", "day5");   // text that identifies each time point in the image name
membraneMaxByTime     = newArray(5000,  8000,  12000);  // membrane display max for each of these time points, same order
membraneMax   = 5000; // used when no time point label is found in the name

nucleiColor   = "Blue";
nucleiMin     = 0;
nucleiMax     = 20000;

// ---- Cellpose settings -------------------------------------------------------
// cellposeEnvPath, cellposeModelPath and the model name (cellposeModelChoice /
// cellposeModelCustom -> cellposeModel) are set from the dialog (see the
// parameters at the very top of this file), not here.
cellposeEnvType      = "conda";
cellposeDiameter     = 40.0; // cpsam still divides by this value internally (0 causes a
                              // ZeroDivisionError) - 40.0 is the default used;
                              // adjust to your typical cell diameter in pixels if segmentation
                              // quality is off, but do not set it to 0.
cellposeExtraArgs    = "";   // extra "key=value " flags appended verbatim if needed

// ---- Per-cell shape measurements (MorphoLibJ "Analyze Regions", run on the label
// image itself - more robust than ROI-polygon-based "Measure" for irregular shapes) --
analyzeRegionsOptions = "area perimeter circularity centroid " +
    "equivalent_ellipse ellipse_elong. convexity max._feret_diameter " +
    "oriented_box_elong. max._inscribed_disc geodesic_elong. average_thickness";

// ---- Columns included in the run-wide SUMMARY table (mean/min/max of each) -----
// The per-cell CSV always keeps every Analyze Regions column (incl. Centroid,
// Ellipse.Center, Ellipse.Orientation, Feret angle); only THIS list controls what
// also gets aggregated into Summary_Cells.csv. Edit freely - names must match the
// per-cell CSV's column headers exactly.
summaryColumns = newArray("Area", "Perimeter", "Circularity", "Ellipse.Radius1",
    "Ellipse.Radius2", "Ellipse.Elong", "ConvexArea", "Convexity", "MaxFeretDiam",
    "OBox.Elongation", "InscrDisc.Radius", "GeodesicElongation", "AverageThickness");

// ---- Elongation color-coding settings -----------------------------------------
arColorColumn = "Ellipse.Elong"; // column from the Analyze Regions table used for
                                  // the color-coded image below (elongation is >= 1)
arColorMin = 0;   // display range starts at 0 although elongation is always >= 1:
                   // with a minimum of 1, round cells (~1) would be as black as the
                   // background (0) in the Fire LUT and disappear
arColorMax = 13.0;

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

// ---- Nuclei-only segments ---------------------------------------------------------
// Cellpose sometimes segments only the nucleus of a cell (sometimes with a little
// cytoplasm). Segments with an area at or below nucleiOnlyMaxArea (calibrated units,
// um2) are treated as nuclei-only, where
//   nucleiOnlyMaxArea = nucleiAreaReference x nucleiOnlyFactor
// nucleiAreaReference = 362 um2, the 95th percentile of nuclear area measured by
// Workflow2 (DAPI nuclei); the factor (> 1) also catches a nucleus plus a little
// cytoplasm. 1.2 -> 434 um2. Set nucleiOnlyFactor = 1 for the plain 95th percentile.
// The original _Cells.csv and _Overlay.tif are always saved with ALL segments.
nucleiAreaReference = 362;                    // um2
nucleiOnlyFactor    = 1.2;
nucleiOnlyMaxArea   = nucleiAreaReference * nucleiOnlyFactor;   // 434.4 um2
saveMarkedOverlay   = true;   // _Overlay_marked.tif: kept cells yellow, nuclei-only in nucleiOnlyColor
saveFilteredOutputs = true;   // _Overlay_filtered.tif and _Cells_filtered.csv without nuclei-only segments
nucleiOnlyColor     = newArray(0, 255, 255);  // RGB, cyan

// ---- Output file naming -------------------------------------------------------
cellsCsvSuffix       = "_Cells.csv";
compositeTifSuffix   = "_Composite.tif";
compositeFlatTifSuffix = "_Composite_flat.tif";
overlayTifSuffix     = "_Overlay.tif";
overlayMarkedTifSuffix   = "_Overlay_marked.tif";
overlayFilteredTifSuffix = "_Overlay_filtered.tif";
cellsFilteredCsvSuffix   = "_Cells_filtered.csv";
arFilteredTifSuffix      = "_AR_filtered.tif";
arFilteredFlatTifSuffix  = "_AR_filtered_flat.tif";
arTifSuffix          = "_AR.tif";
arFlatTifSuffix      = "_AR_flat.tif";
roiZipSuffix         = "_RoiSet.zip";        // saved automatically after every Segment run
roiZipManualSuffix   = "_RoiSet_manual.zip"; // YOU create this: open roiZipSuffix in the ROI
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

// --- resolve the Cellpose model name (custom text field wins over the list) ------
if (isUnset(cellposeModelCustom)) {
    cellposeModel = cellposeModelChoice;
} else {
    cellposeModel = String.trim(cellposeModelCustom);
}
if (processingMode == "Update") {
    print("Cellpose model: " + cellposeModel + "   (ignored in Update mode - ROIs are loaded from disk)");
} else {
    print("Cellpose model: " + cellposeModel);
}

// --- resolve the Results subfolder name (blank -> Results_<model>) ----------------
if (isUnset(resultsSubfolderName)) {
    resultsSubfolderName = "Results_" + cellposeModel;
} else {
    resultsSubfolderName = String.trim(resultsSubfolderName);
}
// keep it a single, filesystem-safe folder name (no slashes, colons, etc.)
resultsSubfolderName = sanitizeName(resultsSubfolderName);
if (resultsSubfolderName == "" || resultsSubfolderName == "_") {
    exit("The Results subfolder name is empty or invalid after cleanup - please enter a valid name.");
}
print("Results subfolder: " + resultsSubfolderName);

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
        exit("No image is open. Open the raw CDH5AF647 image first, or choose a different Image selection mode.");
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

    // --- determine duplicated-stack channel order -----------------------------
    chLow  = minOf(CH_MEMBRANE, CH_NUCLEI);
    chHigh = maxOf(CH_MEMBRANE, CH_NUCLEI);
    if (CH_MEMBRANE == chLow) {
        membraneIdx = 1; nucleiIdx = 2;
    } else {
        membraneIdx = 2; nucleiIdx = 1;
    }

    // --- 1) duplicate the two relevant channels --------------------------------
    selectWindow(imgTitle);
    cleanBaseName = sanitizeName(baseName);
    dupTitle = cleanBaseName + "_dup";
    if (isOpen(dupTitle)) close(dupTitle);
    run("Duplicate...", "duplicate title=" + dupTitle + " channels=" + chLow + "-" + chHigh);

    // --- 2) composite with fixed display ranges --------------------------------
    selectWindow(dupTitle);
    run("Make Composite");
    Stack.setChannel(membraneIdx);
    run(membraneColor);
    imgMembraneMax = getMembraneMaxForImage(lifBaseName + " " + imgTitle,
                                            timePointLabels, membraneMaxByTime, membraneMax);
    print("     Membrane display range: " + membraneMin + " - " + imgMembraneMax);
    setMinAndMax(membraneMin, imgMembraneMax);
    Stack.setChannel(nucleiIdx);
    run(nucleiColor);
    setMinAndMax(nucleiMin, nucleiMax);

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
        cellposeArgs = "env_path=" + cellposeEnvPath + " env_type=" + cellposeEnvType +
                       " model=" + cellposeModel + " model_path=" + cellposeModelPath +
                       " diameter=" + cellposeDiameter +
                       " additional_flags=" + cellposeExtraArgs +
                       " ch1=" + membraneIdx + " ch2=" + nucleiIdx;
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

        // Check whether Cellpose actually found any objects BEFORE calling "Label
        // image to ROIs" - on an all-background (max=0) label image, that command
        // pops up a blocking "ROI Manager: The list is empty" dialog that would
        // otherwise stall/halt an unattended batch run.
        selectWindow(labelTitle);
        getStatistics(dummyArea, dummyMean, labelMin, labelMax);
        if (labelMax == 0) {
            print("WARNING: Cellpose found 0 cells for " + imgTitle + " - recording 0 cells and moving to next image.");
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
        // same order ROIs are stored/measured) so the AR color-coding step downstream has
        // a label image to work with, matching Cellpose's own output convention.
        selectWindow(dupTitle);
        width = getWidth(); height = getHeight();
        getVoxelSize(lbPw, lbPh, lbPd, lbUnit);
        labelTitle = cleanBaseName + "_labels";
        if (isOpen(labelTitle)) close(labelTitle);
        newImage(labelTitle, "16-bit black", width, height, 1);
        // copy the calibration, so Analyze Regions reports um / um2 as in Segment mode
        // (a new image is uncalibrated, which would give areas in pixels)
        setVoxelSize(lbPw, lbPh, lbPd, lbUnit);
        for (r = 0; r < nRoisLoaded; r++) {
            roiManager("select", r);
            setColor(r + 1);
            fill();
        }
        run("Select None");
    }
    nCells = roiManager("count");


    // --- 5) measure per-cell shape descriptors -----------------------------------
    // NOTE: we intentionally do NOT use ROI Manager "Measure" here. Its Fit-Ellipse
    // AR is computed from the ROI's traced polygon and can be misleading for the
    // concave, irregular shapes Cellpose produces (a fairly round cell can get a
    // spuriously high AR). MorphoLibJ's "Analyze Regions", run directly on the
    // label image, computes region-based descriptors (incl. Ellipse.Elong) from
    // the actual pixel mask / inertia matrix, which tracks visual elongation much
    // more reliably (checked visually against the ROI-based fit).
    selectWindow(labelTitle);
    tableTitlesBefore = getList("window.titles");
    run("Analyze Regions", analyzeRegionsOptions);
    tableTitlesAfter = getList("window.titles");
    morphTableName = firstNewTitle(tableTitlesAfter, tableTitlesBefore);
    if (morphTableName == "") {
        print("ERROR: Analyze Regions did not produce a results table for " + imgTitle + " - skipping this image.");
        if (isOpen(dupTitle)) close(dupTitle);
        if (isOpen(labelTitle)) close(labelTitle);
        roiManager("reset");
        run("Clear Results");
        return;
    }
    nCells = Table.size(morphTableName);
    print("     Analyze Regions measured " + nCells + " cell(s).");

    // --- A) save per-cell results -------------------------------------------------
    Table.save(resultsDir + baseName + cellsCsvSuffix, morphTableName);
    print("     Saved: " + baseName + cellsCsvSuffix);

    // --- append to run-wide summary table -----------------------------------------
    appendSummaryRow(lifBaseName, imgTitle, nCells, morphTableName);
    if (saveMarkedOverlay || saveFilteredOutputs) {
        nNucOnlyTable = countAtOrBelow(morphTableName, "Area", nucleiOnlyMaxArea);
        if (nNucOnlyTable >= 0) {
            Table.set("N_NucleiOnly", summaryRowIndex - 1, nNucOnlyTable, "Summary");
            Table.set("N_Cells_Kept", summaryRowIndex - 1, nCells - nNucOnlyTable, "Summary");
            print("     Nuclei-only segments (Area <= " + nucleiOnlyMaxArea + "): " + nNucOnlyTable + " of " + nCells);
        }
    }
    finalizeSummary(resultsDir); // checkpoint: save Summary_Cells.csv after every image, so
                                  // results are never lost even if a later image/file fails

    // --- C) overlay (ROIs, no labels) flattened on composite -----------------------
    saveFlattenedWithOutlines(dupTitle, resultsDir + baseName + overlayTifSuffix);
    print("     Saved: " + baseName + overlayTifSuffix);

    // --- C2) nuclei-only segments: marked overlay and/or filtered outputs ---------
    if (saveMarkedOverlay || saveFilteredOutputs) {
        selectWindow(dupTitle);
        getPixelSize(nuUnit, nuPw, nuPh);
        if (nuUnit == "pixel" || nuUnit == "pixels" || nuUnit == "") {
            print("     WARNING: image has no spatial calibration - nucleiOnlyMaxArea (um2) cannot be applied;");
            print("              marked/filtered outputs skipped for this image.");
        } else {
            nuFlags = nucleiOnlyFlags(dupTitle, nucleiOnlyMaxArea);
            keptIdx = newArray(0);
            nucIdx = newArray(0);
            for (k = 0; k < nuFlags.length; k++) {
                if (nuFlags[k] == 1) nucIdx = Array.concat(nucIdx, k);
                else keptIdx = Array.concat(keptIdx, k);
            }
            if (saveMarkedOverlay) {
                saveFlattenedWithRoiSubsets(dupTitle, resultsDir + baseName + overlayMarkedTifSuffix, keptIdx, nucIdx);
                print("     Saved: " + baseName + overlayMarkedTifSuffix + "  (" + nucIdx.length + " nuclei-only marked)");
            }
            if (saveFilteredOutputs) {
                noIdx = newArray(0);
                saveFlattenedWithRoiSubsets(dupTitle, resultsDir + baseName + overlayFilteredTifSuffix, keptIdx, noIdx);
                print("     Saved: " + baseName + overlayFilteredTifSuffix);
                nDel = saveFilteredCellsCsv(resultsDir + baseName + cellsCsvSuffix,
                                            resultsDir + baseName + cellsFilteredCsvSuffix, nucleiOnlyMaxArea);
                if (nDel >= 0) print("     Saved: " + baseName + cellsFilteredCsvSuffix + "  (" + nDel + " nuclei-only rows removed)");
            }
        }
    }

    // --- D) elongation (Ellipse.Elong) color-coded label image ---------------------
    selectWindow(labelTitle);
    titlesBeforeAR = getList("image.titles");
    call("inra.ijpb.plugins.LabelToValuePlugin.process", "Table=" + morphTableName,
         "Column=" + arColorColumn, "Min=" + arColorMin, "Max=" + arColorMax);
    titlesAfterAR = getList("image.titles");
    arTitle = firstNewTitle(titlesAfterAR, titlesBeforeAR);
    if (arTitle == "") {
        print("WARNING: " + arColorColumn + " color-coded image was not created for " +
              imgTitle + " - Cells.csv/Overlay were still saved, but this image is missing.");
    } else {
        selectWindow(arTitle);
        setMinAndMax(arColorMin, arColorMax);
        run("Fire");
        run("Calibration Bar...", "location=[" + calBarLocation + "] fill=" + calBarFill +
            " label=" + calBarLabelColor + " number=" + calBarNumberLabels +
            " decimal=" + calBarDecimal + " font=" + calBarFontSize +
            " zoom=" + calBarZoom + " overlay");
        saveAs("Tiff", resultsDir + baseName + arTifSuffix);
        rename(arTitle);

        run("Select None"); // never burn a stray selection into the flattened copy
        run("Flatten");
        saveAs("Tiff", resultsDir + baseName + arFlatTifSuffix);
        close(); // close flattened AR copy
        print("     Saved: " + baseName + arTifSuffix + " / " + baseName + arFlatTifSuffix);
    }

    // --- D2) elongation color-coded image WITHOUT nuclei-only segments --------------
    if (saveFilteredOutputs) {
        saveFilteredColorCodedFromAR(arTitle, labelTitle, nucleiOnlyMaxArea,
            resultsDir + baseName + arFilteredTifSuffix, resultsDir + baseName + arFilteredFlatTifSuffix);
    }

    // --- cleanup ------------------------------------------------------------------
    if (isOpen(dupTitle)) close(dupTitle);
    if (isOpen(labelTitle)) close(labelTitle);
    if (isOpen(arTitle)) close(arTitle);
    if (isOpen(morphTableName)) close(morphTableName);
    roiManager("reset");
    run("Clear Results");
    print("  <- Done: " + imgTitle);
}


// =============================================================================
//  6) SUMMARY TABLE (one row per image, whole run)
// =============================================================================

function appendSummaryRow(lifBaseName, imgTitle, nCells, tableName) {
    Table.set("File", summaryRowIndex, lifBaseName, "Summary");
    Table.set("Image", summaryRowIndex, imgTitle, "Summary");
    Table.set("N_Cells", summaryRowIndex, nCells, "Summary");

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
            values = Table.getColumn(colName, tableName);
            Array.getStatistics(values, minV, maxV, meanV, stdV);
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


// Returns the membrane display max for an image: the value whose time point label
// (e.g. "48h") appears in 'name' (case-insensitive), or 'fallback' if none does.
function getMembraneMaxForImage(name, labels, values, fallback) {
    lname = toLowerCase(name);
    for (k = 0; k < labels.length; k++) {
        if (indexOf(lname, toLowerCase(labels[k])) >= 0) return values[k];
    }
    print("     WARNING: no time point label (" + String.join(labels, ", ") + ") in '" + name +
          "' - using the fallback membrane max " + fallback);
    return fallback;
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

// ---- nuclei-only helpers ---------------------------------------------------------

// Returns an array with one entry per ROI in the ROI Manager: 1 if the ROI's area
// (calibrated) is <= maxArea (nuclei-only), else 0.
function nucleiOnlyFlags(srcTitle, maxArea) {
    selectWindow(srcTitle);
    roiManager("Show None");
    n = roiManager("count");
    flags = newArray(n);
    for (i = 0; i < n; i++) {
        roiManager("select", i);
        if (getValue("Area") <= maxArea) flags[i] = 1;
    }
    roiManager("Deselect");
    run("Select None");
    return flags;
}

// Like saveFlattenedWithOutlines(), but draws two sets of ROIs: 'yellowIdx' in yellow
// and 'markIdx' in nucleiOnlyColor. ROIs in neither set are not drawn.
function saveFlattenedWithRoiSubsets(srcTitle, outPath, yellowIdx, markIdx) {
    flatID = makeFlatCopyWithScaleBar(srcTitle);
    run("Line Width...", "line=" + outlineWidth);
    drawRoiSubset(yellowIdx, 255, 255, 0);
    drawRoiSubset(markIdx, nucleiOnlyColor[0], nucleiOnlyColor[1], nucleiOnlyColor[2]);
    saveAndCloseFlatCopy(flatID, outPath, srcTitle);
}

// Draws the outlines of the ROIs with the given indices into the current image.
function drawRoiSubset(idx, r, g, b) {
    if (idx.length == 0) return;
    setForegroundColor(r, g, b);
    roiManager("select", idx);
    roiManager("Draw");             // draws the selected ROIs only
    roiManager("Deselect");
    run("Select None");
}

// Number of rows of 'tableName' whose 'col' value is <= maxValue (-1 if no such column).
function countAtOrBelow(tableName, col, maxValue) {
    heads = split(Table.headings(tableName), "\t");
    found = false;
    for (k = 0; k < heads.length; k++) if (heads[k] == col) found = true;
    if (!found) {
        print("     WARNING: no '" + col + "' column - nuclei-only segments not counted.");
        return -1;
    }
    vals = Table.getColumn(col, tableName);
    c = 0;
    for (k = 0; k < vals.length; k++) if (vals[k] <= maxValue) c++;
    return c;
}

// Saves a copy of the per-cell CSV without the rows whose Area is <= maxArea.
// Returns the number of rows removed (-1 on failure).
function saveFilteredCellsCsv(rawPath, outPath, maxArea) {
    before = getList("window.titles");
    Table.open(rawPath);
    after = getList("window.titles");
    t = firstNewTitle(after, before);
    if (t == "") t = File.getName(rawPath);
    if (countAtOrBelow(t, "Area", maxArea) < 0) {
        if (isOpen(t)) { selectWindow(t); run("Close"); }
        return -1;
    }
    area = Table.getColumn("Area", t);
    nDel = 0;
    for (i = area.length - 1; i >= 0; i--) {
        if (area[i] <= maxArea) {
            Table.deleteRows(i, i, t);
            nDel++;
        }
    }
    Table.update(t);
    Table.save(outPath, t);
    if (isOpen(t)) { selectWindow(t); run("Close"); }
    return nDel;
}

// Saves the elongation color-coded image with the nuclei-only segments shown as background.
function saveFilteredColorCodedFromAR(arSrcTitle, labTitle, maxArea, outPath, outFlatPath) {
    // Copy of the unfiltered elongation image with the nuclei-only segments set to 0
    // (background). The values of all kept cells are copied pixel by pixel, so they
    // keep exactly the same value and color as in the unfiltered image. (Re-running
    // the label-to-value mapping on a size-filtered label image is NOT safe: once
    // labels are removed, table rows and labels can be matched up wrongly.)
    if (arSrcTitle == "" || !isOpen(arSrcTitle)) {
        print("     WARNING: no unfiltered color-coded image - filtered version skipped.");
        return;
    }
    selectWindow(labTitle);
    getPixelSize(lfUnit, lfPw, lfPh);
    if (lfUnit == "pixel" || lfUnit == "pixels" || lfUnit == "") {
        print("     WARNING: label image has no spatial calibration - filtered color-coded image skipped.");
        return;
    }
    maxPixels = floor(maxArea / (lfPw * lfPh));
    titlesBefore = getList("image.titles");
    run("Label Size Filtering", "operation=Greater_Than size=" + maxPixels);
    titlesAfter = getList("image.titles");
    maskTitle = firstNewTitle(titlesAfter, titlesBefore);
    if (maskTitle == "") {
        print("     WARNING: Label Size Filtering failed - filtered color-coded image skipped.");
        return;
    }
    // turn the kept labels into a 0/1 mask (label values themselves are not used)
    selectWindow(maskTitle);
    run("Min...", "value=0");
    run("Max...", "value=1");
    run("32-bit");
    imageCalculator("Multiply create 32-bit", arSrcTitle, maskTitle);
    arfID = getImageID();
    run("Remove Overlay");
    setMinAndMax(arColorMin, arColorMax);
    run("Fire");
    run("Calibration Bar...", "location=[" + calBarLocation + "] fill=" + calBarFill +
        " label=" + calBarLabelColor + " number=" + calBarNumberLabels +
        " decimal=" + calBarDecimal + " font=" + calBarFontSize +
        " zoom=" + calBarZoom + " overlay");
    saveAs("Tiff", outPath);
    selectImage(arfID);
    run("Select None");
    run("Flatten");
    if (getImageID() != arfID) {
        saveAs("Tiff", outFlatPath);
        close();
    }
    selectImage(arfID);
    close();
    if (isOpen(maskTitle)) close(maskTitle);
    print("     Saved: " + File.getName(outPath) + " / " + File.getName(outFlatPath) +
          "  (segments <= " + maxPixels + " pixels set to background)");
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
    print(f, "CH_MEMBRANE = " + CH_MEMBRANE);
    print(f, "CH_NUCLEI = " + CH_NUCLEI);
    print(f, "");
    print(f, "membraneColor = " + membraneColor);
    print(f, "membraneMin = " + membraneMin);
    tpPairs = "";
    for (k = 0; k < timePointLabels.length; k++) {
        tpPairs = tpPairs + timePointLabels[k] + "=" + membraneMaxByTime[k];
        if (k < timePointLabels.length - 1) tpPairs = tpPairs + ", ";
    }
    print(f, "membraneMax per time point = " + tpPairs + "  (fallback membraneMax = " + membraneMax + ")");
    print(f, "nucleiColor = " + nucleiColor);
    print(f, "nucleiMin = " + nucleiMin);
    print(f, "nucleiMax = " + nucleiMax);
    print(f, "");
    print(f, "cellposeEnvPath = " + cellposeEnvPath);
    print(f, "cellposeEnvType = " + cellposeEnvType);
    if (processingMode == "Update") {
        print(f, "cellposeModel = " + cellposeModel + "   (not used in Update mode)");
    } else {
        print(f, "cellposeModel = " + cellposeModel);
    }
    print(f, "cellposeModelPath = " + cellposeModelPath);
    print(f, "cellposeDiameter = " + cellposeDiameter);
    print(f, "cellposeExtraArgs = " + cellposeExtraArgs);
    print(f, "");
    print(f, "analyzeRegionsOptions = " + analyzeRegionsOptions);
    print(f, "");
    summaryColumnsStr = String.join(summaryColumns, ", ");
    print(f, "summaryColumns = " + summaryColumnsStr);
    print(f, "");
    print(f, "arColorColumn = " + arColorColumn);
    print(f, "arColorMin = " + arColorMin);
    print(f, "arColorMax = " + arColorMax);
    print(f, "calBarLocation = " + calBarLocation);
    print(f, "calBarFill = " + calBarFill);
    print(f, "calBarLabelColor = " + calBarLabelColor);
    print(f, "calBarNumberLabels = " + calBarNumberLabels);
    print(f, "calBarDecimal = " + calBarDecimal);
    print(f, "calBarFontSize = " + calBarFontSize);
    print(f, "calBarZoom = " + calBarZoom);
    print(f, "outlineWidth = " + outlineWidth);
    print(f, "nucleiOnlyMaxArea = " + nucleiOnlyMaxArea + " = nucleiAreaReference " + nucleiAreaReference +
              " x nucleiOnlyFactor " + nucleiOnlyFactor + "  (segments with Area <= this are nuclei-only)");
    print(f, "saveMarkedOverlay = " + saveMarkedOverlay + ", saveFilteredOutputs = " + saveFilteredOutputs);
    print(f, "nucleiOnlyColor (RGB) = " + nucleiOnlyColor[0] + "," + nucleiOnlyColor[1] + "," + nucleiOnlyColor[2]);
    print(f, "scaleBar = " + scaleBarLength + " (calibrated units), thickness " + scaleBarThickness +
              " px, font " + scaleBarFontSize + ", " + scaleBarColor + ", " + scaleBarLocation);
    print(f, "");
    print(f, "Output suffixes: " + cellsCsvSuffix + ", " + compositeTifSuffix + ", " +
              compositeFlatTifSuffix + ", " + overlayTifSuffix + ", " + arTifSuffix +
              ", " + arFlatTifSuffix);
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
// window title such as "12h_notch_CONTROL_CDH5AF647.lif - Image1-1".
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
