# The paper's image in this folder

## Source

S. C. Gifford, A. M. Spillane, S. M. Vignes and S. S. Shevkoplyas, "Controlled incremental filtration: a
simplified approach to design and fabrication of high-throughput microfluidic devices for selective enrichment of
particles", *Lab Chip*, 2014, 14, 4496–4505, DOI 10.1039/c4lc00785a.

## Extraction

The image is the publisher's embedded bitmap, extracted without resampling:

```
pdfimages -png -f 5 -l 5 "<the article PDF>" gif5     # gif5-001.png: Fig. 4 (journal page 4500); gif5-000.png is Fig. 3, not used
```

We re-extracted it on 2026-10-18 from the article PDF in the project with poppler `pdfimages` 24.02.0, and the sha256
below reproduces.

| file | figure | size (px) | sha256 |
|---|---|---|---|
| `gif5-001.png` | Fig. 4 (the device schematic on the left, and the inset micrograph and its 300 µm bar on the right) | 829 × 607 | `e70a8b38b87f96393df21e404802f51159d435f163f5e55f3a3b48bc3b523a0b` |

The publisher embedded Fig. 4 as a JPEG (one gray channel, 300 ppi), which `pdfimages -png` decodes. Another JPEG
decoder gives slightly different pixels. For example, Pillow's decoding of the same stream (`pdfimages -j`, sha256
`8e8f7f46dbc32f8cb1368d299a1e6f8e8d9ae09a7154aaf616ca90dd5c936b0b`) differs by up to 9 gray levels in 255. Therefore,
re-extract the figure with `pdfimages -png` to reproduce the sha256 above. The trace and the figures use this PNG.

## The micrograph

The micrograph is the box `[198, 54, 829, 553]` (columns, rows) of the raster, which the trace code finds
(`trace_lib.micrograph_box`), and the figures show exactly that box. The schematic on the left is not used, because
Gifford says that Figs. 4 and 6 "are not drawn to scale, for clarity" (p. 4502).

## Caption, verbatim

The caption of Fig. 4 (p. 4500) reads as follows.

> Fig. 4 CIF-based device for separation of particles with different size. A value of fgap = 3.34 × 10−4 was chosen
> in-between those observed to retain beads with diameter 4.7 μm and diameter 8.3 μm and used to pattern a full-length
> device. Overall design of the device is shown schematically on the left. The inset shows an image of the second and
> fourth legs of the device during operation. Selective retention and concentration of only the larger-size particles
> in the central channel is clearly visible. The depth of the device was 80 μm ± 3 μm; the width of the central
> channel was 100 μm and the final width of the side channels was 300 μm. Arrows indicate direction of flow. Scale bar
> is 300 μm.

## Text, verbatim

The text on p. 4500 reads "The depth of the device shown in Fig. 4 was 80 μm; the width of the central channel was
100 μm and the final width of the side channels was 300 μm, resulting in ~10× enrichment by the end of the device for
particle populations that are maintained in the central channel."

## Use

The article is copyright © The Royal Society of Chemistry 2014, and it carries no Creative Commons license. The image
is reproduced with the approval of Dr. S. S. Shevkoplyas (corresponding author) and the permission of the Royal
Society of Chemistry, and opencif distributes it in `notebooks/figures/`, outside its MIT License (`NOTICE.md`).
