# rc6 patch (1) tried on a scratch copy (after the rc5 reviews)

The verification report's §8.2 patch (1), a table-domain status for the outlet stretch-end term with a refusal in
strict mode, is deferred to rc6 because it changes hashed code.  These files try it on a scratch copy of the rc5
package; nothing here changes `opencif/`.

- `apply_patch1.py SRC_REPO DST_TREE` copies `opencif/` into DST_TREE (with `presets/` and `examples/` linked) and
  applies the patch as the report writes it, each edit by exact string replacement (it stops if a target is not found
  exactly once).  The status is recorded where `_plan_layout` and `_plan_paperclip` make their final outlet plan, and
  the strict refusal is raised in `plan_layout` only on the plan its loop-radius search accepts.
- `run_forced.py TREE R_END|none ARGS…` runs the opencif CLI from that tree.  With a number, it replaces
  `cifsched.outlet_exp_domain` by a test version whose table ends at that width ratio, so the shipped designs (whose
  ratios lie inside the real tables) exercise the refusal path; with `none` the patch runs as written.  It prints how
  often the status was evaluated and at which ratios.
- `cases3.sh` is the run list (written for the author's scratch paths); `outputs/status.txt` gives each case's exit
  code and wall time, `outputs/out_*.txt` the CLI output, and `outputs/build_evidence.txt` the two m02 builds' DXF
  sha256s and record fields.  `specs/` holds the two matrix specs with `gates.mode: evaluate`.

| case | command | result |
|---|---|---|
| m02 (linear, teardrop, strict), table end forced to 1.5 | `check` | exit 5: "the outlet stretch's end lies outside the solved 2-D tables (forced test: width ratio 2.898 beyond 1.50): build it in evaluate or reproduce mode, where the end values are held and reported, or change the outlet"; the status evaluated once, at the settled ratio |
| m04 (linear, ring, strict), forced | `check` | exit 5; evaluated once, ratio 2.194 |
| m12 (paperclip, ring, strict), forced (the `_plan_paperclip` placement) | `check` | exit 5; evaluated once, ratio 2.194 |
| m02 and m12 in evaluate mode, forced | `check` | exit 0, no refusal |
| m02, m04, m12, not forced | `check` | exit 0 |
| m02, not forced | `build` | exit 0; DXF byte-identical to the shipped record's; `design.outlet.stretch_end.table_domain` "inside the solved tables" |
| m02 in evaluate mode, forced | `build` | exit 0; `table_domain` carries the forced note; DXF byte-identical to the shipped one; evaluated four times, all at the settled ratio 2.898 |

Not part of this trial: the verdict clause for evaluate and reproduce mode (`[1-D model: …]`), the unit tests the
report lists, and a design whose loop-radius search rejects a radius (none of the matrix designs does: each strict
case above evaluated the status once).
