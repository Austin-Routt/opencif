# The paper's image in this folder

## Source

H. Xia, B. C. Strachan, S. C. Gifford and S. S. Shevkoplyas, "A high-throughput microfluidic approach for
1000-fold leukocyte reduction of platelet-rich plasma", *Scientific Reports*, 2016, 6, 35943,
DOI 10.1038/srep35943, with its Author Correction (2022, DOI 10.1038/s41598-022-18449-5). Equation 1's w_c(i) and
w_c(ref) were inverted in the article, and the correction gives the corrected form f_gap(i) = (w_c(ref)/w_c(i))²
f_gap(ref).

## Extraction

The image is the publisher's embedded bitmap, extracted without resampling:

```
pdfimages -png -f 5 -l 5 "<the article PDF>" xia5     # xia5-000.png: Fig. 3 (journal page 5)
```

We re-extracted it on 2026-10-18 from the article PDF in the project with poppler `pdfimages` 24.02.0. The sha256
below reproduces, and it is the raster that the trace records of rc6, rc7, and RC8 name.

| file | figure | size (px) | sha256 |
|---|---|---|---|
| `xia5-000.png` | Fig. 3, with (a) the RS-7 device's schematic at a 5 mm bar and (b–i) fluorescence micrographs along the device | 2025 × 3138 (gray) | `e9fb7fbb348567d3ab7606479c83342bcf2f5ee3d6b4c162f7c345e5ef339991` |

The publisher embedded Fig. 3 as a JPEG (one gray channel, 450 ppi), which `pdfimages -png` decodes.

## What the package uses

- Panels c, d, e, f, and g (the micrographs with post rows) are traced at rc6's scale of 1.798 µm per pixel, which is
  the mean of the printed widths, the post columns' separations, and the loop curvature
  (`tree/opencif/data/r6/rc6_values.json`). Panel f is traced at its own 1.774 µm per pixel (its post period). Their
  boxes in the raster are set in `tree/provenance/rc7_traces/code/panels_xia_more.py` (`CROPS`) and `panels_xia.py`
  (`PANELS['g']`).
- Panel a, the whole-device schematic (columns 765–1203, rows 40–2920), is used at its 5 mm bar, read with its
  anti-aliased ends (178.32 px), for the whole device's registration (`xia_fig3a_device.py`).
- Panel b (the inlet's serpentine head, from which the presets draw the head "as depicted") and panels h and i (the two
  collection ports) are not traced, as in rc6.

## Caption, verbatim

The caption of Fig. 3 (p. 5) reads as follows.

> Figure 3. Images of an RS-CIF device performing leukocyte reduction of PRP. (a) Schematic illustration of the RS-7
> device. (b-i) Fluorescent images of DNA-stained leukocytes within a PRP sample as they are progressively concentrated in
> the centre channel, while flowing through the device: (b) device inlet; (c) side channel architecture transition;
> (d) end of the second (right) and fourth (left) loop of the device; (e) the third (right) and fifth (left) legs of the
> device; (f) the seventh leg of the device; (g) device outlet; (h) central channel (retentate) collection port for
> highly-concentrated leukocytes; (i) side channel (filtrate) collection port for leukocyte reduced PRP.
> Filtrate-to-retentate volumetric flow ratio is approximately 15:1. Scale Bars: (a), 5 mm; (b–i), 250 μm.

## The taper, verbatim

The text on p. 4 reads "The degree of centre channel narrowing relative to the degree of side channel widening at
each filtration gap, i, in an RS-CIF device is determined numerically (using formulas for rectangular channels based on
Yang et al.) in order to maintain a relatively constant amount of (maximum) shear experienced along each leg of the
device, while still satisfying the recursive framework constraints of OR-CIF." Fig. 1d(i) adds that "the width of the
central channel, wc(i) gradually decreases until it reaches 150 μm, and remains constant thereafter".

## Use

Xia et al. 2016 is published under a Creative Commons Attribution 4.0 International License (© The Author(s) 2016), and
the figure is reproduced here unchanged, with this attribution. opencif distributes it in `notebooks/figures/` with the
other papers' figures, outside its MIT License (`NOTICE.md`).
