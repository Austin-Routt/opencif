"""test_options_sheet.py: the options sheet against the code and the generator call.

The suite checks that the options sheet (docs/OPTIONS.md) documents every option that the code accepts, and that the
single generator call reaches every part. It is fast because it draws nothing. Run it with
`python3 tests/test_options_sheet.py`."""
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'tools'))
from opencif import spec as SP, api
from opencif import pipeline as _PL_CODE                # noqa: E402  (every suite's log names the code it ran on)
print('opencif %s (code %s)' % (_PL_CODE.CODE['version'], _PL_CODE.CODE['sha256']))
import options_sheet as OS

fails = []
def check(name, cond, info=''):
    print('%-74s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)

missing = [(s, k) for s, k, t in OS.leaves() if k is not None and (s, k) not in OS.KEYS]
check('every schema key has a row in the sheet', not missing, str(missing))
check('every schema section is described', set(SP.SCHEMA) <= set(OS.SECTIONS), str(set(SP.SCHEMA) - set(OS.SECTIONS)))
vals = [(s, k, v) for (s, k), ch in SP.CHOICES.items() for v in ch if (s, k, v) not in OS.VALUES]
check('every allowed value of every enumerated option is described', not vals, str(vals))
extra = [k for k in OS.KEYS if k[0] not in SP.SCHEMA or not isinstance(SP.SCHEMA[k[0]], dict) or k[1] not in SP.SCHEMA[k[0]]]
check('the sheet documents no option the code does not accept', not extra, str(extra))
for (s, k, v), d in OS.VALUES.items():
    if d[2] and not d[2].startswith('tools/') and '*' not in d[2]:
        if not os.path.exists(os.path.join(ROOT, d[2])):
            fails.append('example %s' % d[2])
            print('missing example file', d[2])
        else:
            dd = SP.load(os.path.join(ROOT, d[2]))
            got = dd.get(s, {}).get(k, SP.DEFAULTS.get((s, k)))
            if s == 'meander' and not dd.get('meander', {}).get('enabled', False):
                got = None
            if got != v:
                fails.append('example %s does not use %s.%s = %s' % (d[2], s, k, v))
                print('example %s has %s.%s = %r, not %r' % (d[2], s, k, got, v))
check('every example the sheet cites exists and uses the value it illustrates', not [f for f in fails if f.startswith('example')])
# The defaults that the sheet states and that live outside spec.DEFAULTS are checked against the code.
spA, opA = SP.to_engine(SP.load(api.build_spec(4.0, 0.5, meander='arms')))
spL, opL = SP.to_engine(SP.load(api.build_spec(4.0, 0.5, layout='paperclip')))
check('sheet defaults: arms 3, legs 1 (linear) / auto, port gap 150 (layout.inlet_local), outlet stretch 200, '
      'tolerances 0.01 / 1e-3, DRC on',
      spA.n_arms_fixed == 3 and opA['legs'] == 1 and opL['legs'] == 'auto' and spA.L_out == 200.0
      and opA['cut_tol'] == 0.01 and opA['F_tol'] == 1e-3 and opA['drc'] is True
      and 'clear=150.0' in open(os.path.join(ROOT, 'opencif', 'geometry', 'layout.py')).read()
      and OS.EXTRA_DEFAULTS[('meander', 'arms')] == '`3`' and OS.EXTRA_DEFAULTS[('outlet', 'stretch_um')] == '`200`')
cur = open(OS.OUT).read() if os.path.exists(OS.OUT) else ''
check('docs/OPTIONS.md is current (python3 tools/options_sheet.py)', cur == OS.render())
# Every part argument of generate() reaches its spec section, and the spec it builds passes the schema.
d = api.build_spec(3.0, 0.6, process='dinh2024', chip=(25, 75, 800), layout='paperclip', legs='auto',
                   taper=('constant_shear', 90, 60), meander={'style': 'arms', 'arms': 3, 'channel_um': 28}, seal='require',
                   inlet='filter', outlet={'style': 'ring', 'port_radius_um': 800}, flow_uL_min=20, mode='evaluate', array=2,
                   overrides={'checks': {'cutoff_tolerance': 0.02}, 'physics': {'trim_F': True}})
SP.load(d)
check('generate(): targets', d['targets'] == dict(cutoff_um=3.0, filtrate_fraction=0.6))
check('generate(): named process', d['process'] == api.PROCESSES['dinh2024'])
check('generate(): chip tuple', d['chip'] == dict(width_mm=25.0, length_mm=75.0, margin_um=800.0))
check('generate(): layout and legs', d['layout'] == dict(style='paperclip', legs='auto'))
check('generate(): taper tuple', d['centre_channel'] == dict(taper='constant_shear', width_in_um=90.0, width_min_um=60.0))
check('generate(): meander dict and seal', d['meander'] == dict(style='arms', arms=3, channel_um=28, enabled=True, seal='require'))
check('generate(): inlet, outlet, flow, mode, array', d['inlet'] == dict(style='filter') and d['outlet']['style'] == 'ring'
      and d['operation'] == dict(flow_uL_min=20.0) and d['gates'] == dict(mode='evaluate') and d['array'] == dict(count=2))
check('generate(): overrides merged', d['checks'] == dict(cutoff_tolerance=0.02) and d['physics'] == dict(trim_F=True))
d0 = api.build_spec(4.0, 0.5)
SP.load(d0)
sp0, op0 = SP.to_engine(d0)
check('generate(4, 0.5) defaults: generic19, w_c = D/0.05, pack head sealed, plain inlet, separate outlet, strict',
      d0['process'] == api.PROCESSES['generic19'] and d0['centre_channel']['width_in_um'] == 80.0 and sp0.head == 'pack'
      and op0['seal'] == 'require' and op0['inlet'] == 'plain' and op0['outlet'] == 'separate' and op0['mode'] == 'strict'
      and 'width_in_um = D/0.05' in d0['description'])
dn = api.build_spec(4.0, 0.5, meander=None, start='floor')
check('generate(): meander=None gives a no-head device', dn['meander'] == dict(enabled=False) and dn['start'] == dict(style='floor'))
dr = api.build_spec(4.0, None, filtrate_to_retentate_ratio=1.0)
check('generate(): the flow-ratio form of F', abs(SP.load(dr)['targets']['filtrate_fraction'] - 0.5) < 1e-15)
try:
    api.build_spec(4.0, 0.5, process='nonesuch'); ok = False
except ValueError:
    ok = True
check('generate(): an unknown process name is an error', ok)
try:
    SP.load(api.build_spec(4.0, 0.5, inlet='sieve')); ok = False
except SP.SpecError:
    ok = True
check('generate(): an unknown part value is refused by the schema', ok)
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
