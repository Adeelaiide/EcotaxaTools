# IFCB

The Imaging FlowCytobot (IFCB) is a submersible imaging flow cytometer for nano- and microplankton (Sosik & Olson, 2007). Images are processed with the IFCB software and uploaded to EcoTaxa.

Back to: [Getting started](../getting-started.md) · [Outputs](../outputs.md)

> **The IFCB is listed in the instrument menu but is not yet implemented in the R version.** To process IFCB data, use the [MATLAB version](https://github.com/zoemgt/EcoTaxa-Tools-Box-Pipeline/blob/main/docs/instruments/ifcb.md).

## 1. Required metadata

| EcoTaxa field | Content | Expected unit |
|---|---|---|
| `sample_id` | Sample identifier | |
| `acq_volume_sampled` | Volume of water analysed | mL |
| `acq_resolution_pixel_per_micron` | Image resolution | pixels/µm |
| `acq_id` | Acquisition identifier | |
| `object_lat`, `object_lon` | Sampling position | decimal degrees |
| `object_date`, `object_time` | Sampling date and time | |
| `object_major_axis_length`, `object_minor_axis_length` | Major and minor axes | pixels |
| `object_surface_area` | Area | pixels |
| `object_summed_biovolume` | Biovolume computed by the IFCB software | pixels³ |
| `object_summed_surface_area` | Summed surface area computed by the IFCB software | pixels |
| `object_annotation_category`, `object_annotation_hierarchy`, `object_annotation_status` | Taxonomic annotation and validation status | |

## 2. Volume conversion

The IFCB analyses a known volume of water directly, without concentration:

$$conver = \frac{1}{V_{sampled}}$$

- $V_{sampled}$: volume of water analysed (`acq_volume_sampled`, converted from mL to m³).

## 3. Specificities

- The IFCB export has no `area_exc` field: the riddled-area biovolume is not available.
