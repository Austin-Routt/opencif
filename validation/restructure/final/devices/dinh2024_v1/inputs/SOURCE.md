# The paper's images in this folder

## Source

M. T. P. Dinh, A. Mukhamedshin, K. Abhishek, F. W. Lam, S. C. Gifford and S. S. Shevkoplyas, "Separation of
platelets by size in a microfluidic device based on controlled incremental filtration", *Lab Chip*, 2024, 24,
913–923, DOI 10.1039/d3lc00842h.

## Extraction

The images are the publisher's embedded bitmaps, extracted without resampling:

```
pdfimages -png -f 5 -l 5 "<the article PDF>" dinh5     # dinh5-000.png: Fig. 2 (journal page 917)
pdfimages -png -f 6 -l 6 "<the article PDF>" dinh6     # dinh6-000.png: Fig. 3 (journal page 918)
```

We re-extracted them on 2026-10-18 from the PDF in the project
(the copy "Downloaded by University of Houston on 3/11/2024"), and both sha256 below reproduce.

| file | figure | size (px) | sha256 |
|---|---|---|---|
| `dinh5-000.png` | Fig. 2 (panel A, the layout; panel B, four micrographs) | 982 × 1940 | `a21e2e3ce0710e32410ba4941f348c6a1ef561b748f569e24decce1dea9fbbe6` |
| `dinh6-000.png` | Fig. 3 (panel A, the assembled device; panel B, histograms; panel C, two micrographs) | 1835 × 1316 | `14c01a4ac50adf606a89a7a4292e67e8e70bf5babe9021dfa3d6edea86c90ed5` |

## Captions, verbatim

The captions of Fig. 2 (p. 917) and Fig. 3 (p. 918) read as follows.

> Fig. 2 Design of the microfluidic device. (A) The overall layout of the CIF device. Inset magnifies the aggregate
> filter and the initial part of the meandering side channel sections. (B) Microscopic images of the key sections of the
> device. Width of the meandering side channel segments was 28 μm. Width of each side channel at the start of the
> linear section was wS(in) = 21 μm (for all values of i > f\*gap), widening near the outlet to wS(out) = 44 μm, 48 μm, and
> 57 μm for f\*gap values of 1.04 × 10−4, 1.12 × 10−4, and 1.28 × 10−4, respectively. All filtration gaps in the linear
> section were 19 μm-wide. The middle channel narrowed from 75 μm near the inlet to 52 μm near the outlet in all
> designs. Depth of all channels was 140 μm. The scale bars are 50 μm.

> Fig. 3 CIF device operation. (A) An assembled CIF device (∼3.2 mm × 75.8 mm) filled with orange dye for contrast.
> Dashed rectangles indicate the approximate locations of images shown in panel C. (B) Typical PLT histograms for the
> input sample (PRP), middle channel output (CIFmid), and side channel output (CIFside) of the device. Dashed line
> indicates cell volume of 12 fL. (C) PLTs retained in the middle channel accumulate gradually near the filtration gaps:
> (i) virtually no PLT accumulation in the beginning of the linear section of the device; (ii) significant PLT
> accumulation is visible further along the device. Arrows indicate flow direction. Scale bars are 100 μm.

(The phrase "for all values of i > f\*gap" is as printed, which we checked on the rendered page as well as in the text
layer.)

## Which images are traced

The trace uses the four micrographs of Fig. 2B, each at its own 50 µm bar. Fig. 3C's two micrographs are included but
not traced, because at their printed 100 µm bars the post pitch reads 51.2 to 51.4 µm, against 45.1 µm on all three
post-bearing panels of Fig. 2B at theirs. The authors' intensity-profile axes in Fig. 3C confirm the 45.1 µm pitch
scale, so Fig. 3C's 100 µm bars are about 12% to 14% short and Fig. 2B's bars are right (register INV-5, rc7). Fig. 3C
could be traced at the pitch scale, although this package does not trace it.

## Use

The article is copyright © The Royal Society of Chemistry 2024, and it carries no Creative Commons license. The images
are reproduced with the approval of Dr. S. S. Shevkoplyas (corresponding author) and the permission of the Royal
Society of Chemistry, and opencif distributes them in `notebooks/figures/`, outside its MIT License (`NOTICE.md`).
