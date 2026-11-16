# EcoTaxa Tools Box Pipeline (R)

The EcoTaxa Tools Box is an open-source pipeline designed to process and explore taxonomically annotated images hosted on [EcoTaxa](https://ecotaxa.obs-vlfr.fr/). An equivalent version is available in MATLAB: [zoemgt/EcoTaxa-Tools-Box-Pipeline](https://github.com/zoemgt/EcoTaxa-Tools-Box-Pipeline).

It reads the tsv files exported from EcoTaxa for one instrument and computes key ecological descriptors for each sample. Each instrument has specific variable names and units: the pipeline accounts for them, so that the descriptors are computed in the same way and expressed in the same units whatever the instrument (FlowCam, PlanktoScope, ZooScan). For each sample, it computes:

- abundance (ind m⁻³) and biovolume (mm³ m⁻³), per taxonomic annotation and per grouping level (living vs non-living, plankton groups, trophic groups);
- biovolume size spectra (BSS) and normalized biovolume size spectra (NBSS);
- standard graphical analyses at three levels: per-sample (or station) plots, global plots summarizing the entire dataset, and NBSS size spectra analyses.

## How it works

Run a single function: **`EcotaxaTools()`**. It guides you through successive steps:

1. **Instrument**: selection of the instrument (FlowCam, PlanktoScope, ZooScan)
2. **Processing**: the metadata of all tsv files are checked and can be edited; biovolumes, abundances and size spectra are computed for each object and sample; the results are saved as csv tables in the 'summary' folder.
3. **Analysis**: figures are produced for the metadata, the whole project, for each sample, and are saved as images in the 'graph' folder.

At the end, `EcotaxaTools()` offers an optional module for multivariate analyses (PCA, PCoA, hierarchical clustering). See [Getting started](docs/getting-started.md#6-optional-multivariate-analyses).

## Supported instruments

| Instrument | Size range | Instrument-specific notes |
|---|---|---|
| [ZooScan](docs/instruments/zooscan.md) | > 200 µm | Several size fractions per sample, sub-sampling via `acq_sub_part` |
| [FlowCam](docs/instruments/flowcam.md) | 20–200 µm | One fraction per sample |
| [PlanktoScope](docs/instruments/planktoscope.md) | 20–200 µm | Concentration and dilution factors |

Each instrument page lists the mandatory EcoTaxa metadata, their expected units and the volume conversion formula.

## Requirements

- R (≥ 4.1) and **RStudio** (recommended: the folder selection and the metadata editor open in RStudio).
- The R packages listed in `DESCRIPTION` (installed automatically with the package): dplyr, ggplot2, DataEditR, svDialogs, sf, rnaturalearth, vegan, FactoMineR…

## Installation

```r
install.packages("devtools")
devtools::install_github("Adeelaiide/EcotaxaTools", force = TRUE)
```

## Quick start (with the test datasets)

Small test datasets (3 samples per instrument) are available in the `test_data/` folder of the [MATLAB repository](https://github.com/zoemgt/EcoTaxa-Tools-Box-Pipeline/tree/main/test_data):

| Instrument | Folder | Source |
|---|---|---|
| ZooScan | `test_data/zooscan/` | Tara Pacific expedition (Mériguet et al., 2025), [EcoTaxa project #11292](https://ecotaxa.obs-vlfr.fr/prj/11292) |
| FlowCam | `test_data/flowcam/single_fraction/` | Tonga, [EcoTaxa project #2841](https://ecotaxa.obs-vlfr.fr/prj/2841) |
| PlanktoScope | `test_data/planktoscope/` | Ada-II 2023 expedition (Mériguet et al., 2026, under review), [EcoTaxa project #10056](https://ecotaxa.obs-vlfr.fr/prj/10056) |

1. In RStudio, run:
   ```r
   library(EcotaxaTools)
   EcoTaxa_Tool_Box()
   ```
2. Choose the instrument in the list (e.g. **PlanktoScope**).
3. Select the test folder of your instrument (e.g. `test_data/planktoscope/`).
4. Check the metadata in the editor; edit them if necessary and save changes with "Synchronise" and "Done"; then answer the dialogs (edited metadata, new taxa).
5. The tables and figures are written in a new folder `EcoTaxa_data_analysis_<date>` created inside the data folder.

Each dialog is explained step by step in [Getting started](docs/getting-started.md), and instrument-specific steps in the instrument pages.

## Before using your own data

- Export from EcoTaxa with **Export → General export**, and set *Separate (in ZIP sub-directories) output by* to **Sample**.
- Check that all mandatory metadata are filled in and in the expected units (see the instrument pages). Missing metadata are replaced by 1, which produces wrong values without stopping the pipeline.

## Documentation

- [Getting started](docs/getting-started.md): step-by-step walkthrough of the pipeline
- [Outputs](docs/outputs.md): content of the csv tables and figures, variables and units
- Instrument pages: [ZooScan](docs/instruments/zooscan.md) · [FlowCam](docs/instruments/flowcam.md) · [PlanktoScope](docs/instruments/planktoscope.md) · [IFCB](docs/instruments/ifcb.md)

## Citation

If you use this pipeline, please cite:

Mériguet, Z., Perruchon, A., Zaccomer, H., et al. (2026) Studying plankton communities from annotated quantitative images on EcoTaxa: an automated plankton community processing and analysis pipeline. *Journal of Plankton Research* (under review).

## License and contact

GPL (≥ 3). See [LICENSE](LICENSE.md).

Adélaïde Perruchon (adelaide.perruchon@senckenberg.de), Hugo Zaccomer (hugo.zaccomer@imev-mer.fr), Laboratoire d'Océanographie de Villefranche-sur-Mer.
