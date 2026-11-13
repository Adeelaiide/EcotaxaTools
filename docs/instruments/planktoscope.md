# PlanktoScope

The PlanktoScope is a low-cost, open-source imaging platform for microplankton (Pollina et al., 2022). Images are segmented by the PlanktoScope software and uploaded to EcoTaxa. In the pipeline, choose **PlanktoScope** in the instrument list.

Back to: [Getting started](../getting-started.md) · [Outputs](../outputs.md)

## 1. Required metadata

These fields are used in calculations: they must be filled in, in the expected units. Missing values are replaced by 1.

| EcoTaxa field | Content | Expected unit |
|---|---|---|
| `sample_id` | Sample identifier | |
| `sample_total_volume` | Volume of water collected (e.g. filtered by the net) | L |
| `sample_concentrated_sample_volume` | Volume of the concentrated sample | mL |
| `sample_dilution_factor` | Dilution factor (1 if absent or empty, see section 2) | no unit |
| `acq_imaged_volume` | Volume imaged by the PlanktoScope | mL |
| `process_pixel` | Pixel size | µm/pixel |
| `acq_id` | Acquisition identifier | |
| `object_lat`, `object_lon` | Sampling position | decimal degrees |
| `object_date`, `object_time` | Sampling date and time | |
| `object_major`, `object_minor` | Major and minor axes | pixels |
| `object_area`, `object_area_exc` | Area, area excluding holes | pixels |
| `object_annotation_category`, `object_annotation_hierarchy`, `object_annotation_status` | Taxonomic annotation and validation status | |

These fields are read but not used in calculations. They must be present in the export and can be empty: `sample_operator`, `acq_minimum_mesh`, `acq_maximum_mesh`. The field `acq_celltype` is optional.

## 2. Volume conversion

For each sample:

$$conver = \frac{V_{concentrated}}{V_{imaged} \times V_{total} \times D}$$

- $V_{concentrated}$: volume of the concentrated sample (`sample_concentrated_sample_volume`, converted from mL to m³);
- $V_{imaged}$: volume imaged (`acq_imaged_volume`, converted from mL to m³), summed over the acquisitions of the sample;
- $V_{total}$: volume of water collected (`sample_total_volume`, converted from L to m³);
- $D$: dilution factor (`sample_dilution_factor`, 1 if absent or empty).

Abundance (ind m⁻³) = number of objects × `conver`.

## 3. Instrument-specific steps

No additional window. If a sample was imaged in several acquisitions, they appear on separate lines in the metadata editor, and their imaged volumes are summed.

## 4. Specificities

- All three biovolume methods are available (plain area, riddled area, ellipsoid).
- Small, round objects close to the resolution limit of the PlanktoScope are difficult to identify: they are often annotated as unidentified.

## 5. Test dataset

Folder `test_data/planktoscope/` of the [MATLAB repository](https://github.com/zoemgt/EcoTaxa-Tools-Box-Pipeline/tree/main/test_data): 3 samples from the Ada-II 2023 expedition (Mériguet et al., 2026, under review), [EcoTaxa project #10056](https://ecotaxa.obs-vlfr.fr/prj/10056).

1. Run `EcotaxaTools()` and choose **PlanktoScope**.
2. Select `test_data/planktoscope/`.
3. Follow the steps described in [Getting started](../getting-started.md).
