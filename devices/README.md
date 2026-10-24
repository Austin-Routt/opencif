# The device kits

The four kits in this folder reproduce the published devices with opencif and register each photomask on the paper's
own micrographs. They are Gifford et al. (2014, Fig. 4 and Fig. 6), Xia et al. (2016, RS-7), and Dinh et al. (2024,
V1). Each kit builds the device's photomasks from its presets, checks that every photomask is byte for byte the
release's (which is RC8's), traces the photomask onto the micrographs with the trace code in
`provenance/rc7_traces/code`, compares every value of the new registration with RC8's record, and draws the figures.

## The publishers' figures

The micrographs are the publishers'. Gifford et al. (2014) and Dinh et al. (2024) are © the Royal Society of Chemistry,
and Xia et al. (2016) is published under CC BY 4.0. The figures are in `notebooks/figures/`, reproduced with the
Royal Society of Chemistry's permission and with Xia et al.'s attribution, outside the MIT License (`NOTICE.md`). Each kit's
`MANIFEST.json` names the figures that it needs and gives the sha256 of each, and `inputs/SOURCE.md` says where each
figure comes from and how it was extracted from the paper's PDF. A kit refuses a figure whose sha256 differs. The
figures that a kit draws contain the publishers' images and carry the same terms.

## Running a kit

From the repository's folder, run

```
python3 tools/device_package.py gifford2014_fig4 --rasters notebooks/figures --out OUT_DIR
```

where `notebooks/figures` holds the figures (or any folder of figures extracted as `inputs/SOURCE.md` says). The command assembles the package in `OUT_DIR/opencif_gifford2014_fig4/`
as a clean room, with the repository's package and the files that the kit needs under `tree/`, the kit's own files
beside it, and the figures under `inputs/`, each checked against its sha256. It then runs the device's script there in
its own process. The script runs its steps in order, and it stops with exit code 2 at the first departure from the
release's records. A step can be run alone by naming it after the options (for example `check build`).

The steps are the same for every device. The check step confirms that the package's opencif is the release (its code
stamp equals the one in `tree/validation/RELEASE.json`) and that every file is the one that `MANIFEST.json` names. The
build step builds each photomask with `pipeline.run` and compares its sha256 with the release's photomask
(`reference/<name>_dxf_sha256.txt`). The trace step registers the photomask on each micrograph and compares every value
of the registration with RC8's trace record (`reference/`). The render step draws the zoomed micrograph panels with the
mask's edges over them, the whole device with the micrographs placed by their registrations, and that drawing cut into
rows, each as a PDF and a PNG. Some devices add steps of their own, which their scripts (`make_<device>.py`) describe;
Gifford's Fig. 4 adds the taper-direction test (V13).

## What the kits show and what they do not

A kit shows that the release builds the same photomasks as RC8 and that the photomasks register on the papers'
micrographs as RC8's did, with every number of the registration equal. It does not measure anything new, because the
build is the generator's and the registrations are the trace code's. `docs/REPRODUCTIONS.md` gives each paper's printed
quantities beside the reproduced ones, `docs/LIMITATIONS.md` (item 28) and `provenance/RC8_DISCREPANCIES.md` record what
the registrations do not match (for example, the neighboring leg's side channels on Gifford's Fig. 4 inset, GIF4-1), and
`docs/PRESET_PROVENANCE.md` gives the source of every value of the presets.

## Files

| path | contents |
|---|---|
| `_shared/devicekit.py` | the steps that every kit runs |
| `<device>/make_<device>.py` | the device's script, with everything that is specific to the device |
| `<device>/MANIFEST.json` | every file of the package with its source and sha256, and the figures that `--rasters` must supply |
| `<device>/reference/` | the release's photomask sha256 sums and exit codes (`EXPECT.json`) |
| `<device>/inputs/SOURCE.md` | the source of each figure and how it was extracted |
