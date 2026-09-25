# Example World-Camera Images

This directory contains manually selected example frames from the following recordings:

- `indoor1`, `indoor2`, and `indoor3`: `FLIC_2001` `walkIndoor`
- `outdoor1` and `outdoor2`: `FLIC_2001` `walkOutdoor`
- `outdoor3`: `FLIC_1054` `sitBiopond`
- `planetarium1`: the planetarium recording

These are completely raw world-camera frames. No processing, correction, linearization, debayering, or other transformation has been applied. The `outdoor3` example uses frame 5000 from subject/recording `FLIC_1054 sitBiopond`.

## Raw-frame indices

The zero-based global indices below count frames physically present in the naturally ordered raw world-camera chunks. They do not include synthetic frames for timestamp gaps. These indices are the durable selection rule: `code/defineWorldCameraCalibration/empiricalDataPrep/populateData.ipynb` reads this block to regenerate the paired `*_AGCandMS_*.mat` files from the raw recordings.

<!-- populateData:example-frame-indices:start -->
| Frame | Source recording | Zero-based global raw-frame index |
| --- | --- | ---: |
| `indoor1` | FLIC_2001 walkIndoor | 4777 |
| `indoor2` | FLIC_2001 walkIndoor | 42993 |
| `indoor3` | FLIC_2001 walkIndoor | 21713 |
| `outdoor1` | FLIC_2001 walkOutdoor | 4512 |
| `outdoor2` | FLIC_2001 walkOutdoor | 40615 |
| `outdoor3` | FLIC_1054 sitBiopond | 5000 |
| `planetarium1` | Fels Planetarium | 96978 |
<!-- populateData:example-frame-indices:end -->

## Digital-gain values

Record the digital-gain (`DGain`) value associated with each frame below.

| Frame | Source recording | `DGain` |
| --- | --- | --- |
| `indoor1` | `FLIC_2001` `walkIndoor` | 4.830781965284408 |
| `indoor2` | `FLIC_2001` `walkIndoor` | 5.982647368887896 |
| `outdoor1` | `FLIC_2001` `walkOutdoor` | 1.0 |
| `outdoor2` | `FLIC_2001` `walkOutdoor` | 1.0 |
| `planetarium1` | Planetarium recording | 4.2608509465799225 |
