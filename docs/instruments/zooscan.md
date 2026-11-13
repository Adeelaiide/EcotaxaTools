# ZooScan

The ZooScan is an automated imaging system for preserved mesozooplankton (> 200 µm) (Gorsky et al., 2010). Images are segmented with ZooProcess and uploaded to EcoTaxa. In the pipeline, choose **ZooScan** in the instrument list.

Back to: [Getting started](../getting-started.md) · [Outputs](../outputs.md)

## 1. Required metadata

These fields are used in calculations: they must be filled in, in the expected units. Missing values are replaced by 1.

| EcoTaxa field | Content | Expected unit |
|---|---|---|
| `sample_id` | Sample identifier | |
| `sample_tot_vol` | Volume of water filtered by the net | m³ |
| `acq_id` | Size fraction identifier | |
| `acq_sub_part` | Sub-sampling factor of the fraction (e.g. 4 if 1/4 of the sample was scanned) | |
| `process_particle_pixel_size_mm` | Pixel size | mm/pixel |
| `object_lat`, `object_lon` | Sampling position | decimal degrees |
| `object_date`, `object_time` | Sampling date and time | |
| `object_major`, `object_minor` | Major and minor axes | pixels |
| `object_area`, `object_area_exc` | Area, area excluding holes | pixels |
| `object_annotation_category`, `object_annotation_hierarchy`, `object_annotation_status` | Taxonomic annotation and validation status | |

These fields are read but not used in calculations. They must be present in the export and can be empty: `sample_scan_operator`, `sample_barcode`, `acq_min_mesh`, `acq_max_mesh`.

## 2. Volume conversion

For each fraction, the conversion factor from imaged objects to concentration is:

$$conver = \frac{S}{V_{filtered}}$$

- $S$: sub-sampling factor of the fraction (`acq_sub_part`; e.g. 4 if 1/4 of the sample was scanned);
- $V_{filtered}$: volume of water filtered by the net (`sample_tot_vol`, m³).

Abundance (ind m⁻³) = number of objects × `conver`.

## 3. Instrument-specific steps

No additional window. Each size fraction (`acq_id`) is an acquisition: its metadata appear on a separate line in the metadata editor. Check that the sub-sampling factor `acq_sub_part` is correct for each fraction.

## 4. Specificities

- All three biovolume methods are available (plain area, riddled area, ellipsoid).
- Touching objects must be separated in ZooProcess before upload: they bias abundance and size estimates (Vandromme et al., 2012).

## 5. Test dataset

Folder `test_data/zooscan/` of the [MATLAB repository](https://github.com/zoemgt/EcoTaxa-Tools-Box-Pipeline/tree/main/test_data): 3 samples from the Tara Pacific expedition (Mériguet et al., 2025), [EcoTaxa project #11292](https://ecotaxa.obs-vlfr.fr/prj/11292).

1. Run `EcotaxaTools()` and choose **ZooScan**.
2. Select `test_data/zooscan/`.
3. Follow the steps described in [Getting started](../getting-started.md).
