# What the license covers

The MIT License in `LICENSE` covers opencif: the package `opencif/`, its tests, tools, specifications, notebooks,
documents, device kits, and the scripts in `provenance/` and `validation/`. Four parts of this repository are
distributed with opencif under other terms.

- `reference/CIF_v2.py` (version 2.1.0) is the laboratory's design code, on which opencif's kernels and physics
  auditor are based. It is included byte for byte, so that the tests can compare opencif with it, and the package does
  not import it. The MIT License does not extend to it. Its rights remain with its authors, and any use beyond running
  this repository's tests needs their permission.
- `foundation/` holds Appendix A, Appendix B and the glossary of A. H. Routt's dissertation, frozen on 2026-09-29, with
  the laboratory's design specification v5 and the list of errata for the next unfreeze. They are included so that
  the tests can trace the code's constants to them. They are © their authors, all rights reserved, and they may be
  read with this repository but not republished without permission.
- `reference_glossary/` holds the glossary's LaTeX preamble, which the tests read, under the same terms as
  `foundation/`.
- `notebooks/figures/` holds six figures of the published papers, as the publishers embedded them in the articles'
  PDFs, which the trace suites, the device kits and the traces notebook read (`notebooks/figures/README.md`). Gifford
  et al. (2014) and Dinh et al. (2024) are © the Royal Society of Chemistry, and their figures are reproduced with its
  permission. Xia et al. (2016) is published under CC BY 4.0, and its figures are reproduced unchanged with attribution.
  The MIT License does not extend to them: any other use of Gifford et al.'s and Dinh et al.'s figures needs the Royal
  Society of Chemistry's permission, and any use of Xia et al.'s needs the attribution that CC BY 4.0 requires. The figures that the kits and
  the traces notebook draw contain these images and carry the same terms.
