# Outputs

The pipeline saves its results as csv tables (semicolon-separated, decimal comma: read them with `readr::read_csv2()`), in the folder `EcoTaxa_data_analysis_<date>` created inside your data folder. The function `EcotaxaTools()` also returns the processed database as an R data frame.

```r
library(readr)
db  <- read_csv2("EcoTaxa_data_analysis_<date>/summary/Processed database.csv")
bss <- read_csv2("EcoTaxa_data_analysis_<date>/summary/BSS.csv")
```

| File | Content |
|---|---|
| `<sample_id>.csv` | One table per sample, one line per object: measurements, biovolumes, size classes, conversion factor |
| `summary/BSS.csv` | Abundance and biovolume per sample, taxonomic annotation, size class and biovolume method |
| `summary/Processed database.csv` | `BSS.csv` merged with the metadata used and the taxonomic grouping |
| `summary/Abundance.csv` | Abundance per taxonomic annotation (rows) and sample (columns), with totals |
| `summary/Biovolume_elli.csv`, `Biovolume_plain.csv`, `Biovolume_riddled.csv` | Biovolume per taxonomic annotation and sample, for each biovolume method |
| `summary/summary_all.csv` | Abundance and biovolume per sample and taxonomic annotation, all methods, without size classes |
| `summary/taxo.csv` | Taxonomic levels and grouping of each annotation |
| `summary/metadata_used.csv` | Metadata used for the calculations |

## 1. Size classes

All descriptors are computed on a common grid of **213 biovolume size classes**, identical for all instruments. The classes follow a geometric progression: the upper limit of each class is 2^0.25 times its lower limit (Platt & Denman, 1978). The grid covers biovolumes from 10⁻¹² to about 10⁴ mm³, i.e. equivalent spherical diameters (ESD) from about 0.1 to 27 000 µm. Each instrument only fills the classes corresponding to its own size range.

| Variable | Content | Unit |
|---|---|---|
| `class` | Number of the size class (1 to 213) | |
| `min`, `max` | Lower and upper biovolume limits of the class | mm³ |
| `norm` | Width of the class (`max − min`), used to normalize the spectra | mm³ |
| `ESD` | Upper limit of the class expressed as ESD (in `Processed database.csv`) | µm |

## 2. Per-object tables (`<sample_id>.csv`)

Each line is an object, with its EcoTaxa fields and the following variables:

| Variable | Content | Unit |
|---|---|---|
| `major`, `minor` | Major and minor axes | mm |
| `area`, `area_exc` | Area, area excluding holes | mm² |
| `BV_plain` | Biovolume from object area | mm³ |
| `BV_riddled` | Biovolume from excluded area | mm³ |
| `BV_elli` | Ellipsoidal biovolume (from major and minor axes) | mm³ |
| `ESD_plain`, `ESD_riddled` | ESD from area and from excluded area | mm |
| `size_class.id_*`, `size_class.min_*`, `size_class.max_*`, `size_class.norm_*` | Size class of the object for each biovolume method | |
| `conver.uniqueID`, `conver.sample` | Conversion factor from imaged objects to concentration, per acquisition and per sample | m⁻³ |

The three biovolume methods are described in Vandromme et al. (2012). The ellipsoidal method is recommended: it is used for all figures.

## 3. Size spectra and summary tables (`BSS.csv`, `Processed database.csv`)

Each line corresponds to one sample, one taxonomic annotation, one size class and one biovolume method:

| Variable | Content | Unit |
|---|---|---|
| `sample_id` | Sample identifier | |
| `object_annotation_category`, `object_annotation_hierarchy` | Taxonomic annotation | |
| `type` | Biovolume method: `elli`, `plain` or `riddled` | |
| `AB` | Abundance in the size class | ind m⁻³ |
| `BV` | Biovolume in the size class (biovolume size spectrum, BSS) | mm³ m⁻³ |

The NBSS is obtained by dividing the biovolume of each class by its width: **NBSS = `BV` / `norm`** (mm³ m⁻³ mm⁻³). It is not stored as a column; the function `BSS_table(x, NBSS = TRUE)` returns it directly.

`Processed database.csv` also contains the metadata used (one set per sample), the sample number (`sample_num`, chronological order) and the taxonomic grouping:

| Variable | Content |
|---|---|
| `Type` | Living status (living, non_living, temporary) |
| `Sub_type` | Plankton group |
| `Value` | Trophic group code (see [Getting started](getting-started.md#4-taxonomic-grouping)) |
| `Trophic_lvl` | Trophic group name |

## 4. Interpreting the outputs

- Abundances and biovolumes depend on the quality of the annotation and of the metadata: check volumes and units before any interpretation. Missing metadata are replaced by 1.
- Size classes at the edges of the instrument's size range are under-sampled.
- The figures produced by the pipeline are exploratory. Choose the descriptors and methods adapted to your scientific question before any publication.
