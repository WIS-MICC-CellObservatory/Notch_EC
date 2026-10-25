# Sample data

## Input images

Raw .lif files are usually too large for GitHub (100 MB per file). The sample images are therefore single multi-channel TIFFs, exported from the .lif files in Fiji (open the .lif series, then File › Save As › Tiff), with the original channel order and calibration:

- Workflow 1 (cell elongation): 4 channels, membrane (CDH5) = channel 3, nuclei = channel 4; file names contain `CDH5AF647`
- Workflow 2 (nuclear Hes1): 5 channels, Hes1 = channel 3, DAPI = channel 5; file names contain `Hes`

File names keep the time point and condition (e.g. `day5_notch_act_CDH5AF647.lif - Image 2.tif`), so both the macros and the Python analysis recognise them. To run the examples, choose *Whole folder* or *Single file* and set *File extension(s)* to `.tif`.

The full raw data are available at [Zenodo / BioImage Archive link].

## Example outputs (`Out/`)

The main README shows six example outputs from day 5, as PNG images (GitHub cannot display TIFFs):

| File | Content |
|---|---|
| `day5_notch_control_CDH5AF647_Image_3_AR_filtered_flat.png` | Workflow 1, day 5 control: cells color-coded by elongation, nuclei-only segments as background |
| `day5_notch_act_CDH5AF647_Image_2_AR_filtered_flat.png` | Workflow 1, day 5 Notch activation: as above |
| `day5_notch_inh_CDH5AF647_Image_5_AR_filtered_flat.png` | Workflow 1, day 5 Notch inhibition: as above |
| `day5_notch_control_cd31AF555_Hes1AF647_Image_3_HesIntensity_flat.png` | Workflow 2, day 5 control: nuclei color-coded by background-corrected Hes1 intensity |
| `day5_notch_act_cd31AF555_Hes1AF647_Image_2_HesIntensity_flat.png` | Workflow 2, day 5 Notch activation: as above |
| `day5_notch_inh_cd31AF555_Hes1AF647_Image_2_HesIntensity_flat.png` | Workflow 2, day 5 Notch inhibition: as above |

Each PNG is the macro's `_flat.tif` output for that image, saved from Fiji with File › Save As › PNG (optionally reduced first with Image › Scale, e.g. to 1000 × 1000 pixels). To show other images, add the PNG here and change the file name in the image table at the top of the main README.
