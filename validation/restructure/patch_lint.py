"""opencif 1.0.0: what pyflakes finds after the restructure (imports and assignments left unused; comments that still
describe the old layout).  Pure removals: no expression with an effect is touched (the full rebuild checks every
photomask byte for byte).  Applied after patch_names.py.   usage: python -I patch_lint.py TREE"""
import os
import sys

T = sys.argv[1]
n_done = 0


def edit(rel, old, new, count=1):
    global n_done
    p = os.path.join(T, rel)
    s = open(p, encoding='utf-8').read()
    assert s.count(old) == count, (rel, s.count(old), old[:100])
    s = s.replace(old, new)
    open(p, 'w', encoding='utf-8').write(s)
    n_done += 1


edit('opencif/published.py', '\nimport numpy as np\n', '\n')
edit('opencif/compose.py', 'from shapely.geometry import LineString, Point, box\n', 'from shapely.geometry import Point, box\n')
edit('opencif/checks/drc.py', 'import numpy as np\nimport ezdxf\n', 'import ezdxf\n')
edit('opencif/checks/drc.py', '    # connected components of the whole fluid (cheap: union of all)\n    t1 = time.time()\n',
     '    # connected components of the whole fluid (cheap: union of all)\n')
edit('opencif/physics/schedule.py', 'import math, sys\n', 'import math\n')
edit('opencif/physics/schedule.py',
     'Frozen physics used (v114 record; the official kernels are imported from CIF_v2.py so the numbers are the engine\'s):',
     'Frozen physics used (v114 record; the kernels are opencif.physics.kernels, which return the design code\'s numbers bit\n'
     'for bit: tests/test_kernels.py):')
edit('opencif/physics/schedule.py', '  * duct law r(w, h): Gifford eqn 5, 50 odd terms, sorted (official CIF_v2.r);',
     '  * duct law r(w, h): Gifford eqn 5, 50 odd terms, sorted (kernels.r = CIF_v2.r);')
edit('opencif/figures/proof.py', 'from opencif.physics import schedule as S  # noqa: E402\n', '')
edit('opencif/figures/proof.py', "    rc_ = rep['validation']['realised_cutoff']\n    k0 = max(0, NM - 5)\n",
     "    k0 = max(0, NM - 5)\n")
edit('opencif/figures/drawing.py', "matplotlib.use('Agg')\nimport matplotlib.pyplot as plt\n", "matplotlib.use('Agg')\n")
edit('opencif/geometry/depicted.py', '    s, sa, w = spec.s, spec.s_a, spec.w_arm\n', '    sa, w = spec.s_a, spec.w_arm\n')
edit('opencif/geometry/layout.py', '''    from scipy.optimize import brentq
    spec = sc.spec
    R = S.R_of(sc, k)
    Lam = spec.Lam_L
    PL = Lam + math.pi * R
    pt = Lam / 2.0 + math.pi * R / 4.0
''', '''    spec = sc.spec
    R = S.R_of(sc, k)
    Lam = spec.Lam_L
    PL = Lam + math.pi * R
''')
edit('opencif/geometry/layout.py', "        (xa, ya), (xb, yb_) = p.p[-2], p.p[-1]\n", "        xb, yb_ = p.p[-1]\n")
edit('opencif/geometry/shapes.py', "            R2, t2 = layout_arms(n2)\n            xb = sc.x[k] - sa / 2.0\n",
     "            R2, t2 = layout_arms(n2)\n")
edit('opencif/geometry/dinh_inlet.py', '\nimport math\n', '\n')
edit('tests/test_bend_fd.py', 'from opencif import build  # noqa: F401  (puts the engine on the path)\n', '')
edit('tests/test_composition.py', 'import copy\nimport os\n', 'import os\n')
edit('tests/test_options_sheet.py', '# the defaults the sheet states that live outside spec.DEFAULTS, against the code\nfrom opencif.physics import schedule as S\n',
     '# the defaults the sheet states that live outside spec.DEFAULTS, against the code\n')
edit('tests/test_units.py', 'from opencif.geometry import layout as L\nfrom opencif import published as D\n',
     'from opencif.geometry import layout as L\n')
edit('tools/cfd_kit.py', 'import json\nimport math\nimport os\n', 'import json\nimport os\n')
edit('tools/run_batch.py', 'import os, subprocess, sys, time, json\n', 'import os, subprocess, sys, time\n')
edit('tools/dissertation_figures.py', 'from matplotlib.patches import Rectangle, Circle      # noqa: E402\n',
     'from matplotlib.patches import Rectangle              # noqa: E402\n')
edit('tools/dissertation_figures.py',
     "import opencif                                        # noqa: E402,F401  (puts the engine, the pickles' classes, on the path)",
     "import opencif                                        # noqa: E402,F401  (the package the pickles' classes belong to)")
print('patch_lint: %d edits' % n_done)
