# tests/test_figure_anchors.py on simulated trees (after the rc5 reviews and the checks of their verification)

Each tree links the final package (`opencif/`), presets and provenance and holds every run record; they differ only in
what is missing.  `exit_codes.txt` gives each run's exit status (0 ALL OK, 1 FAILED, 4 INCOMPLETE).

| log | tree | expected |
|---|---|---|
| `no_pickles.log` | records without their run pickles (as the evaluation archive ships them) | INCOMPLETE, 18 ok, 17 not run, exit 4 |
| `source_anchors_only.log` | the same, `--source-anchors-only` | ALL OK, 17 checks, exit 0 |
| `source_anchors_only_no_provenance.log` | the same without `provenance/`, `--source-anchors-only` | INCOMPLETE (the length-note checks cannot run), exit 4 |
| `one_pickle_missing.log` | every pickle but m02's | INCOMPLETE (the record check ran on 45, not on m02), exit 4 |
| `one_preset_value_changed.log` | a copy of the presets with one `length_bound` value changed by 0.1 mm, `--source-anchors-only` | FAILED, exit 1 |

The full tree's own run is `logs/test_figure_anchors.log` (ALL OK, 35 checks).
