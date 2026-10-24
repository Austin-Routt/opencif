# The paper's image in this folder

## Source

S. C. Gifford, A. M. Spillane, S. M. Vignes and S. S. Shevkoplyas, "Controlled incremental filtration: a
simplified approach to design and fabrication of high-throughput microfluidic devices for selective enrichment of
particles", *Lab Chip*, 2014, 14, 4496–4505, DOI 10.1039/c4lc00785a.

## Extraction

The image is the publisher's embedded bitmap, extracted without resampling:

```
pdfimages -png -f 7 -l 7 "<the article PDF>" gif7     # gif7-000.png: Fig. 6 (journal page 4502)
```

We re-extracted it on 2026-10-18 from the article PDF in the project with poppler `pdfimages` 24.02.0. The sha256
below reproduces, and it is the raster that the trace records of rc6, rc7, and RC8 name.

| file | figure | size (px) | sha256 |
|---|---|---|---|
| `gif7-000.png` | Fig. 6 (the device schematic on the left, and the close-up micrograph of the ends of four legs on the right) | 829 × 615 | `182826e7aaccf14ce9ad72f1c6411ad662bf96c19b5814dec242ebce81a00379` |

The publisher embedded Fig. 6 as a JPEG, like Fig. 4, and `pdfimages -png` decodes it, whereas another JPEG decoder
(Pillow on the `pdfimages -j` stream) gives slightly different pixels. Therefore, re-extract the figure with
`pdfimages -png` to reproduce the sha256 above. The trace and the figures use this PNG.

## The micrograph

The micrograph is the box `[199, 59, 828, 541]` (columns, rows) of the raster. The trace code finds it
(`trace_lib.micrograph_box`) and uses the same box as its crop, and the figures show exactly that box. The schematic
on the left is not used, because Gifford says that Figs. 4 and 6 "are not drawn to scale, for clarity" (p. 4502).

## The scale

The figure's 250 µm bar (rows 548–553, below the box) runs off the raster's right edge, so its visible 54.6 px only
bound the scale from above (at most 4.5807 µm per pixel). Therefore, the trace of rc6, rc7, and RC8 takes the scale
from the printed 125 µm central channel, between the inner tip lines of each leg's two post rows. That scale is
4.5248 µm per pixel (4.5143–4.5302 over the legs), inside the bar's bound
(`tree/provenance/rc6_scripts/gifford_fig6_measurement.json`).

## Caption, verbatim

The caption of Fig. 6 (p. 4502) reads as follows.

> Fig. 6 CIF-based device for enrichment of platelets in platelet-rich plasma (PRP). Overall design of the device is
> shown schematically on the left. The value of f gap for this device was 1.04 × 10−4. A close-up view of a
> representative functional device in operation, showing incremental concentration of platelets in successive legs of
> the CIF array, is shown on the right. The depth of the device was 150 μm ±5 μm; the width of the central channel was
> 125 μm and the final width of the side channels was 140 μm such that the ratio of total flow into the side channel
> collection port was ~2.5× that into the central collection port. The device consistently retained >80% of platelets
> in the central channel, representing a ~3× enrichment in platelet concentration overall.

## Text, verbatim

The text on p. 4502 reads as follows.

> We patterned a full-length (~40 cm, in total) CIF device with f gap =1.04 × 10−4 (w c = 125 μm, d = 150 μm) in order
> to reliably concentrate the vast majority of platelets in a flowing suspension of PRP into the central channel. We
> selected the final width of the side-channels (w s(i) = 144 μm) to achieve our ultimate desired level of platelet
> enrichment after a single pass through the device, however other levels of enrichment are easily selected and
> implemented by setting the final value of w s(i) during the recursive design process. Fig. 6 illustrates the overall
> design and operation of this device, which consistently produced ~3× enrichment of platelets with 80–85% yield. In a
> typical experiment (shown in Fig. 6), the input platelet count (PLT) was 384 × 10³ μL, and the output PLT for the
> central channel was 1039 × 10³ μL (enriched platelets) and for the side channels PLT was 93 × 10³ μL (purified
> plasma). The yield of platelets in the central (enrichment) channel was 83.2%, with a total throughput of 55 μL min−1,
> generated simply from hanging the PRP reservoir at a height of 5 ft above the CIF device (i.e. no pump required).

(The PDF's text layer drops the asterisk of f*gap and renders w s(i) with a stray glyph, and both are restored above as
the printed symbols. The two printed final widths differ, with 140 µm in the caption and 144 µm in the text. The design
preset uses the caption's width and the aligned preset the text's, and register GIF6-2 weighs them against the
micrograph.)

## Use

The article is copyright © The Royal Society of Chemistry 2014, and it carries no Creative Commons license. The image
is reproduced here with the approval of Dr. S. S. Shevkoplyas (corresponding author) for use in A. H. Routt's
dissertation. It belongs in the dissertation's evidence package rather than in opencif's public repository.
