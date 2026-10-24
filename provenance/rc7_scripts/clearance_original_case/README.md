# The case that started the clearance guarantee, at the final code (register BUILD-2)

`g4B_fold20.yaml` is the spec A. H. Routt's figure came from (2026-10-19 16:08): Gifford Fig. 4's printed design at a
20 µm floor with a packing head (W 2400 µm), the legs in path order at the inset's 818 µm
(`layout.leg_spacing_um: 818`), leg 4 ending 8.2 mm short.  At rc7's first freeze (code 4a624626f570918c) it was drawn
and written: the head lay across leg 2, and the build FAILED seven checks after a 14-minute build.

At the final code 52ba91eddb3f07b8:
- `opencif check g4B_fold20.yaml` → exit 5 in 50 s (`check_52ba91eddb3f07b8.out`): the drawn parts collide, 'inlet and
  head' and 'leg 2' 0.0 µm apart, with the remedies (the generator's own spacing for this head, 1668 µm);
- `opencif build g4B_fold20.yaml` → exit 5 in 66 s, REFUSED (layout), no photomask written (`g4B_fold20.json`, the
  refusal record; `g4B_fold20_refused_layout.png`, the two parts outlined in red; `build_52ba91eddb3f07b8_tail.out`).

The same parts at the generator's own spacing are `provenance/rc7_traces/variants/gifford2014_fig4_fold20.yaml`
(legs 1-2 1672 µm apart): built, exit 0, every check passing.
