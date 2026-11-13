# FlowCam

The FlowCam is an imaging flow cytometer for microplankton (Sieracki et al., 1998). Images are processed with VisualSpreadsheet and uploaded to EcoTaxa. In the pipeline, choose **FlowCam** in the instrument list.

Back to: [Getting started](../getting-started.md) · [Outputs](../outputs.md)

## 1. Required metadata

These fields are used in calculations: they must be filled in, in the expected units. Missing values are replaced by 1.

| EcoTaxa field | Content | Expected unit |
|---|---|---|
| `sample_id` | Sample identifier | |
| `sample_initial_col_vol_m3` | Volume of water collected (e.g. filtered by the net) | m³ |
| `sample_conc_vol_ml` | Volume of the concentrated sample | mL |
| `acq_fluid_volume_imaged` | Volume of fluid imaged by the FlowCam | mL |
| `process_pixel` | Pixel size | µm/pixel |
| `acq_id` | Acquisition identifier | |
| `object_lat`, `object_lon` | Sampling position | decimal degrees |
| `object_date`, `object_time` | Sampling date and time | |
| `object_major`, `object_minor` | Major and minor axes | pixels |
| `object_area`, `object_area_exc` | Area, area excluding holes | pixels |
| `object_annotation_category`, `object_annotation_hierarchy`, `object_annotation_status` | Taxonomic annotation and validation status | |

These fields are read but not used in calculations. They must be present in the export: `object_id`, `acq_raw_image_total`, `acq_min_esd`, `acq_max_esd`. The fields `sample_volconc` and `acq_celltype` are optional.

## 2. Volume conversion

For each sample:

$$conver = \frac{V_{concentrated}}{V_{imaged} \times V_{collected}}$$

- $V_{concentrated}$: volume of the concentrated sample (`sample_conc_vol_ml`, converted from mL to m³);
- $V_{imaged}$: volume of fluid imaged (`acq_fluid_volume_imaged`, converted from mL to m³), summed over the acquisitions of the sample;
- $V_{collected}$: volume of water collected (`sample_initial_col_vol_m3`, m³).

Abundance (ind m⁻³) = number of objects × `conver`.

## 3. Instrument-specific steps

No additional window.

The R version processes one fraction per sample, with all images validated. Datasets in which only a subset of the images was validated for one size fraction (e.g. 30 % of the objects < 45 µm in the Tara Pacific FlowCam datasets) are currently handled by the MATLAB version only (see its [FlowCam page](https://github.com/zoemgt/EcoTaxa-Tools-Box-Pipeline/blob/main/docs/instruments/flowcam.md)).

## 4. Specificities

- All three biovolume methods are available (plain area, riddled area, ellipsoid).

## 5. Test dataset

Folder `test_data/flowcam/single_fraction/` of the [MATLAB repository](https://github.com/zoemgt/EcoTaxa-Tools-Box-Pipeline/tree/main/test_data): 3 samples from Tonga, [EcoTaxa project #2841](https://ecotaxa.obs-vlfr.fr/prj/2841).

1. Run `EcotaxaTools()` and choose **FlowCam**.
2. Select `test_data/flowcam/single_fraction/`.
3. Follow the steps described in [Getting started](../getting-started.md).
