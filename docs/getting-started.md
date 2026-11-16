# Getting started

This page walks you through a complete run of the EcoTaxa Tools Box (R version), from the EcoTaxa export to the figures. The pipeline is interactive: each step is a dialog window or an editor in RStudio. The tables below list the windows in the order in which they appear, and what to answer.

Instrument-specific metadata are described in the instrument pages: [ZooScan](instruments/zooscan.md) · [FlowCam](instruments/flowcam.md) · [PlanktoScope](instruments/planktoscope.md) .

## 1. Before you start

**Export your data from EcoTaxa**

1. In your EcoTaxa project, go to **Export → General export**.
2. Keep the default settings, except *Separate (in ZIP sub-directories) output by*, which must be set to **Sample**.
3. Unzip the archive and gather all the tsv files (one per sample) in a single folder.

**Check your metadata**

All mandatory metadata must be filled in, in the expected units. The list depends on the instrument: see the instrument pages. Missing metadata are replaced by 1: the pipeline does not stop, but the values are wrong. A warning listing the missing metadata is printed in the R console for each sample.

**First run?** Use a test dataset (see the [README](../README.md#quick-start-with-the-test-datasets)).

## 2. Launch the pipeline

In RStudio:

```r
library(EcotaxaTools)
EcoTaxa_Tool_Box()
```

## 3. Processing the tsv files

| Window | What to do |
|---|---|
| *Instrument* | Choose PlanktoScope, FlowCam or ZooScan. |
| *Please select the data directory* | Click OK, then select the folder with your tsv files (one file per sample). All tsv files in this folder are processed. |
| Warning on acquisitions (only if a sample has several acquisitions) | Check that all acquisitions belong to the same sample. If not, create a separate tsv file for the additional acquisition and restart. |
| *You can now edit metadata* | Click OK. The metadata of all samples open in an editor in RStudio (see below). |
| Metadata editor | Check and correct the metadata, then click **Synchronise** and **Done**. The original tsv files are not modified. Do not change `unique_id`. |
| *Do you want to use the edited metadata?* | **Yes** to use the metadata you have just checked; **No** to select another metadata table (`.csv`), e.g. one saved during a previous run. |
| *These taxa do not exist in the database. Would you like to create them?* (only if new taxa are found) | See section 4. |

A new folder `EcoTaxa_data_analysis_<date>` is created inside the data folder: all tables and figures are saved there (see section 7).

During processing, the pipeline computes for each object the biovolume (three methods) and its size class, then for each sample the abundance, the biovolume and the size spectra, per taxonomic annotation. The content of the tables is described in [Outputs](outputs.md).

**Metadata editor.** The editor shows one line per acquisition (`unique_id`). The original metadata are saved in `metadata/original_metadata.csv`, and the edited metadata in `metadata/edited_metadata_<date>.csv`, which you can reload in a later run (answer **No** to *Do you want to use the edited metadata?* and select this file). Samples are numbered (`sample_num`) in chronological order, from the edited dates and times.

## 4. Taxonomic grouping

Each taxonomic annotation is assigned to three grouping levels, defined in the table `trophic_affiliation_of_organisms` provided with the package:

- **living / non-living** (`Type`);
- **plankton group** (`Sub_type`, e.g. bacillariophyta, copepoda, detritus);
- **trophic group** (`Value`), coded as follows:

| Code | Trophic group |
|---|---|
| -1 | none (not feeding) |
| 1 | phototrophs |
| 1.5 | mixotrophs |
| 2 | grazers |
| 2.5 | omnivorous |
| 3 | predators |
| 3.5 | unknown |

**New taxa.** If your dataset contains taxa that are not yet in the table, a window lists them and asks *Would you like to create them?*

| Choice | What happens |
|---|---|
| Yes | An editor opens: for each new taxon, choose its category (living status and plankton group) and its trophic group, then click **Synchronise** and **Done**. The completed table is saved in `summary/trophic_affiliation_of_organisms.csv`. |
| No, then *import an existing database?* Yes | Select an existing affiliation table (`.csv`), e.g. one saved during a previous run. |
| No, then No | The new taxa are classified as *temporary*, with an unknown trophic group (3.5). |

The table of the package itself is not modified: to reuse your additions, reload `summary/trophic_affiliation_of_organisms.csv` in the next runs.

**Temporary categories (t00X).** Objects that cannot be identified are usually annotated in broad categories such as "other<plankton", which may end up containing a large and heterogeneous part of the plankton. As an alternative, used at the Quantitative Imaging Platform of Villefranche-sur-Mer (PIQv), visually consistent groups of unidentified objects can be annotated in EcoTaxa as temporary morphotypes (`t001`, `t002`…). This keeps them isolated: they can be identified later, or analysed separately.

Using these categories is optional. In the R version, objects annotated in the EcoTaxa category *temporary* are automatically classified as *temporary*, with an unknown trophic group (3.5).

## 5. Analysis

The figures are produced automatically, without additional windows:

- **metadata**: sampling map (with sample numbers and track) and number of acquisitions per sample;
- **global plots**: abundance and biovolume per sample (total, living, non-living), relative composition of the living per plankton group, maps of abundance, biovolume and Shannon index, NBSS per sample and for the whole project, biovolume composition per size class, trophic pyramid of all samples;
- **per-sample plots**: one synoptic sheet per sample.

The global plots and per-sample plots use the ellipsoidal biovolume and focus on living organisms (detritus and non-living objects excluded).

## 6. Optional multivariate analyses

At the end, the pipeline asks *Do you want to continue the pipeline with multivariate analysis?* This module is an example tool for exploring the data: normalization, ordination and clustering methods should be selected according to your dataset and scientific objectives.

| Window | What to do |
|---|---|
| *What data should be used for analysis?* | Abundance or biovolume, per taxonomic annotation, plankton group or trophic group. |
| *Do you want to delete some groups or taxa before the analysis?* | **Yes**: an editor opens, set each group to *Keep* or *Remove*, then click **Synchronise** and **Done**. **No**: you can import a selection saved during a previous run, or keep all groups. |
| *Do you want to normalise your data?* | None, log(x+1), double cube root, relative abundance, Hellinger, or centred and scaled. |
| Distance for PCoA | Euclidean (recommended after Hellinger transformation) or Bray-Curtis (recommended otherwise). |
| *Do you want to keep N clusters?* | Keep the number of clusters proposed automatically for the hierarchical clustering, or enter your own. |

The module runs a PCA, a PCoA and a hierarchical clustering, and projects the results on maps.

## 7. Where are the results?

All files are saved in the folder `EcoTaxa_data_analysis_<date>`, created inside your data folder:

| File or folder | Content |
|---|---|
| `<sample_id>.csv` | One table per sample: each object with its biovolumes and size classes (see [Outputs](outputs.md)) |
| `metadata/` | Original and edited metadata |
| `summary/` | Processed database, size spectra and summary tables (see [Outputs](outputs.md)) |
| `graph/` | Sampling map and number of acquisitions per sample |
| `graph/global raw analysis/` | Figures for the whole dataset |
| `graph/raw analysis per sample/` | One synoptic sheet per sample |
| `multivariates analysis/` | Outputs of the optional multivariate module |

These figures are direct, unedited outputs, intended for data exploration. Check them and interpret them critically before using them in a publication.

## 8. Common problems

| Problem | Likely cause |
|---|---|
| No file is processed | The export was not separated by sample, or the wrong folder was selected. |
| Warning *These metadata do not have value. Default is 1* | A mandatory metadata field is missing or empty: correct it in the metadata editor or in EcoTaxa. |
| Abundances or biovolumes are implausible | A metadata field is in the wrong unit (e.g. mL instead of L), or was replaced by 1. Check the instrument page. |
| The folder selection or the editor does not open | Run the pipeline from RStudio. |
