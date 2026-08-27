# Example World-Camera Images

This directory contains manually selected example frames from the following recordings:

- `indoor_1.tiff` and `indoor_2.tiff`: `FLIC_2001` `walkIndoor`
- `outdoor_1.tiff` and `outdoor_2.tiff`: `FLIC_2001` `walkOutdoor`
- `planetarium_1.tiff`: the planetarium recording

These are completely raw world-camera frames. No processing, correction, linearization, debayering, or other transformation has been applied.

## Raw-frame indices

The zero-based global indices below count frames physically present in the naturally ordered raw world-camera chunks. They do not include synthetic frames for timestamp gaps. Run the example-frame discovery section in `code/defineWorldCameraCalibration/populateData.ipynb` while the raw recordings are available to replace this provisional block.

<!-- populateData:example-frame-indices:start -->
| Frame | Source recording | Zero-based global raw-frame index |
| --- | --- | ---: |
| `indoor_1.tiff` | FLIC_2001 walkIndoor | 4777 |
| `indoor_2.tiff` | FLIC_2001 walkIndoor | 42993 |
| `indoor_3.tiff` | FLIC_2001 walkIndoor | 21713 |
| `outdoor_1.tiff` | FLIC_2001 walkOutdoor | 4512 |
| `outdoor_2.tiff` | FLIC_2001 walkOutdoor | 40615 |
| `planetarium_1.tiff` | Fels Planetarium | 96978 |
<!-- populateData:example-frame-indices:end -->

## Digital-gain values

Record the digital-gain (`DGain`) value associated with each frame below.

| Frame | Source recording | `DGain` |
| --- | --- | --- |
| `indoor_1.tiff` | `FLIC_2001` `walkIndoor` | 4.830781965284408 |
| `indoor_2.tiff` | `FLIC_2001` `walkIndoor` | 5.982647368887896 |
| `outdoor_1.tiff` | `FLIC_2001` `walkOutdoor` | 1.0 |
| `outdoor_2.tiff` | `FLIC_2001` `walkOutdoor` | 1.0 |
| `planetarium_1.tiff` | Planetarium recording | 4.2608509465799225 |
