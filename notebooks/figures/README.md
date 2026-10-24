# The published papers' figures

These six images are figures of the three published papers that opencif reproduces. They are the publishers' embedded
bitmaps, extracted from the articles' PDFs with poppler's `pdfimages -png` without resampling. The table gives each
one's sha256, and the trace suites and the device kits refuse an image whose sha256 differs from the one they pin.
The traces notebook reads this
folder by default; the suites read it with `OPENCIF_PAPER_RASTERS=notebooks/figures`, and the device kits with
`tools/device_package.py DEVICE --rasters notebooks/figures`.

**These figures are not covered by opencif's MIT License** (`NOTICE.md`). Gifford et al. (2014) and Dinh et al. (2024)
are © the Royal Society of Chemistry, and their figures are reproduced here with its permission; any other use needs
the Royal Society of Chemistry's permission. Xia et al. (2016) is published under a Creative Commons Attribution 4.0
International License (© The Author(s) 2016), and its figures are reproduced unchanged with the attribution below.

| file | article | figure | extraction | sha256 |
|---|---|---|---|---|
| `dinh5-000.png` | Dinh et al. 2024 | Fig. 2 (A, the layout; B, four micrographs) | `pdfimages -png -f 5 -l 5 <PDF> dinh5` | `a21e2e3ce0710e32410ba4941f348c6a1ef561b748f569e24decce1dea9fbbe6` |
| `dinh6-000.png` | Dinh et al. 2024 | Fig. 3 (A, the assembled device; B, histograms; C, two micrographs) | `pdfimages -png -f 6 -l 6 <PDF> dinh6` | `14c01a4ac50adf606a89a7a4292e67e8e70bf5babe9021dfa3d6edea86c90ed5` |
| `gif5-001.png` | Gifford et al. 2014 | Fig. 4 (the schematic, and the inset micrograph with its 300 µm bar) | `pdfimages -png -f 5 -l 5 <PDF> gif5` | `e70a8b38b87f96393df21e404802f51159d435f163f5e55f3a3b48bc3b523a0b` |
| `gif7-000.png` | Gifford et al. 2014 | Fig. 6 (the schematic, and the micrograph of the ends of four legs) | `pdfimages -png -f 7 -l 7 <PDF> gif7` | `182826e7aaccf14ce9ad72f1c6411ad662bf96c19b5814dec242ebce81a00379` |
| `xia3-000.png` | Xia et al. 2016 | Fig. 1 (read only for its scale bars, `provenance/rc7_traces/code/measure_scales.py`) | `pdfimages -png -f 3 -l 3 <PDF> xia3` | `5f4132c81c92399b4211daab3ca001e12084b43cb7d7fb0591d037824468f67f` |
| `xia5-000.png` | Xia et al. 2016 | Fig. 3 (a, the RS-7 schematic; b–i, micrographs along the device) | `pdfimages -png -f 5 -l 5 <PDF> xia5` | `e9fb7fbb348567d3ab7606479c83342bcf2f5ee3d6b4c162f7c345e5ef339991` |

The articles:

- S. C. Gifford, A. M. Spillane, S. M. Vignes and S. S. Shevkoplyas, "Controlled incremental filtration: a simplified
  approach to design and fabrication of high-throughput microfluidic devices for selective enrichment of particles",
  *Lab Chip*, 2014, 14, 4496–4505, DOI 10.1039/c4lc00785a. © The Royal Society of Chemistry 2014.
- M. T. P. Dinh, A. Mukhamedshin, K. Abhishek, F. W. Lam, S. C. Gifford and S. S. Shevkoplyas, "Separation of platelets
  by size in a microfluidic device based on controlled incremental filtration", *Lab Chip*, 2024, 24, 913–923, DOI
  10.1039/d3lc00842h. © The Royal Society of Chemistry 2024.
- H. Xia, B. C. Strachan, S. C. Gifford and S. S. Shevkoplyas, "A high-throughput microfluidic approach for 1000-fold
  leukocyte reduction of platelet-rich plasma", *Scientific Reports*, 2016, 6, 35943, DOI 10.1038/srep35943, published
  under a Creative Commons Attribution 4.0 International License (https://creativecommons.org/licenses/by/4.0/); the
  figures are unchanged.

Some figures were embedded as JPEG streams, so another JPEG decoder gives slightly different pixels: re-extract with
`pdfimages -png` to reproduce the sha256 values above (each kit's `inputs/SOURCE.md` gives the details).
