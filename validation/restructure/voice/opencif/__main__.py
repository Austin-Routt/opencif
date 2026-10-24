"""python -m opencif build SPEC.yaml [--out DIR] [--no-render] [--no-drc]
   python -m opencif check SPEC.yaml          (validate the spec, plan the design, lay out its outline in memory, and
                                              evaluate the gates, with no photomask written; exit 5 refused, 3 outside
                                              the envelope, 0 inside, 6 spec error)
   python -m opencif audit RUN_DIR [--name NAME] [--out FILE] [--no-drc]
                                             (run the physics auditor on a finished build's photomask file from
                                              scratch, which plans and draws again, reads the DXF back, measures it,
                                              and checks P0-P10, R0, M14, and R6b; exit 0 pass, 3 pass with M14
                                              outside the envelope, 2 fail, 4 incomplete)
   python -m opencif generate D F [--process P] [--chip WxL[xMARGIN_UM]] [--layout L] [--legs N|auto]
                              [--taper constant|constant_shear] [--width-in UM] [--width-min UM]
                              [--meander pack|arms|hairpins|narrow|none] [--seal require|auto|allow]
                              [--inlet plain|filter|dinh_filter] [--outlet separate|teardrop|ring]
                              [--start sealed|floor] [--flow UL_MIN] [--mode strict|evaluate|reproduce] [--array N]
                              [--name NAME] [--out DIR] [--set section.key=value ...] [--check] [--no-render]
                              [--no-drc]
        (opencif.generate from the command line, with D in um and F as the targets; docs/OPTIONS.md lists every
        option)
Exit codes: 0 pass, 2 fail, 3 outside the validity envelope (evaluate/reproduce), 4 incomplete, 5 refused,
6 spec error."""
import argparse
import json
import sys

from . import pipeline as PL
from . import spec as SP


def main(argv=None):
    """Parse the command line, run the subcommand, and return its exit code (6 for a spec error)."""
    from . import __version__
    from . import api
    C = SP.CHOICES
    ap = argparse.ArgumentParser(prog='opencif', description=__doc__.split('Exit codes')[0].strip(),
                                 epilog='Exit codes' + __doc__.split('Exit codes')[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', action='version', version='opencif %s (code %s)' % (__version__, PL.CODE['sha256']))
    ap.add_argument('--list-processes', action='store_true',
                    help='list the named processes (--process NAME) with their parameters and sources, and exit')
    sub = ap.add_subparsers(dest='cmd')
    b = sub.add_parser('build', help='build a spec file: plan, draw, measure and assess (writes the DXF, record and figures)')
    b.add_argument('spec', help='a spec YAML or JSON file (docs/OPTIONS.md)')
    b.add_argument('--out', help='output folder (default runs/<name>)')
    b.add_argument('--no-render', action='store_true', help='skip the three figures')
    b.add_argument('--no-drc', action='store_true', help='skip the DRC (the verdict is then INCOMPLETE, exit 4)')
    c = sub.add_parser('check', help='plan a spec, lay out its outline in memory (no photomask written) and evaluate its gates '
                                     'and its parts\' clearance (exit 5 refused, 3 outside, 0 inside)')
    c.add_argument('spec', help='a spec YAML or JSON file')
    au = sub.add_parser('audit', help='run the physics auditor (P0-P10, R0, M14, R6b) on a finished build\'s photomask '
                                      'file, from scratch (exit 0 pass, 3 pass with M14 outside the envelope, 2 fail, '
                                      '4 incomplete)')
    au.add_argument('run_dir', help='the build\'s folder (its record <name>.json and photomask <name>.dxf)')
    au.add_argument('--name', help='the build\'s name, when the folder holds more than one record')
    au.add_argument('--out', help='also write the audit as JSON to this file')
    au.add_argument('--no-drc', action='store_true', help='skip the DRC (P9 is then incomplete, exit 4)')
    gp = sub.add_parser('generate', help='design from the two targets D and F (opencif.generate on the command line)',
                        description='Design a CIF device from its cutoff D and filtrate fraction F; every part has a '
                                    'default (docs/OPTIONS.md lists every option, its default and its source).')
    gp.add_argument('cutoff_um', type=float, help='target cutoff diameter D (um)')
    gp.add_argument('filtrate_fraction', type=float, help='target filtrate fraction F (0 < F < 1)')
    gp.add_argument('--process', metavar='NAME|JSON',
                    help='named process (%s; default generic19) or a JSON process block, e.g. \'{"depth_um": 80, ...}\''
                    % ', '.join(api.PROCESSES))
    gp.add_argument('--chip', metavar='WxL[xMARGIN]', help='chip width x length in mm, optional margin in um (e.g. 25x75)')
    gp.add_argument('--layout', choices=C[('layout', 'style')], help='device layout (default linear)')
    gp.add_argument('--legs', metavar='N|auto', help='number of legs, or auto (the fewest that fit the chip)')
    gp.add_argument('--taper', choices=C[('centre_channel', 'taper')], help='centre-channel taper (default constant)')
    gp.add_argument('--meander', choices=C[('meander', 'style')] + ('none',),
                    help='meander head style, or none (default pack)')
    gp.add_argument('--seal', choices=C[('meander', 'seal')], help='gap-1 sealing policy of the head (default require)')
    gp.add_argument('--inlet', choices=C[('inlet', 'style')], help='inlet style (default plain)')
    gp.add_argument('--outlet', choices=C[('outlet', 'style')], help='outlet style (default separate)')
    gp.add_argument('--start', choices=C[('start', 'style')], help='side-channel start without a head (default sealed)')
    gp.add_argument('--mode', choices=C[('gates', 'mode')], help='gates mode (default strict)')
    gp.add_argument('--name', help='device name (default cif_D<D>_F<F>)')
    gp.add_argument('--out', help='output folder')
    gp.add_argument('--width-in', type=float, metavar='UM', help='centre width at gap 1 in um (default D/0.05)')
    gp.add_argument('--width-min', type=float, metavar='UM', help='narrowest centre width of a constant-shear taper (um)')
    gp.add_argument('--flow', type=float, metavar='UL_MIN', help='operating flow in uL/min (default 30; sets the Reynolds gates)')
    gp.add_argument('--array', type=int, metavar='N', help='also compose N parallel copies')
    gp.add_argument('--set', action='append', default=[], metavar='SECTION.KEY=VALUE',
                    help='any spec option, e.g. --set process.aspect_ratio_max=25 (repeatable)')
    gp.add_argument('--check', action='store_true', help='plan, gate and draw the outline in memory (its collisions), '
                                                          'nothing written')
    gp.add_argument('--no-render', action='store_true', help='skip the three figures')
    gp.add_argument('--no-drc', action='store_true', help='skip the DRC (verdict INCOMPLETE)')
    a = ap.parse_args(argv)
    if a.list_processes:
        for k, v in api.PROCESSES.items():
            print('%-12s %s\n             source: %s' % (k, json.dumps(v), api.PROCESS_SOURCES.get(k, '')))
        return 0
    if not a.cmd:
        ap.print_usage()
        return 6
    try:
        return _main(a)
    except SP.SpecError as e:
        print('SPEC ERROR: %s' % e)
        return 6
    except (ValueError, TypeError) as e:
        if a.cmd == 'generate' and getattr(e, '_opencif_args', False):
            print('SPEC ERROR: %s' % e)
            return 6
        raise


def _audit(a):
    """Run `opencif audit RUN_DIR`, the physics auditor on a finished build's photomask file, and return its exit
    code."""
    from .checks import auditor as AU
    r = AU.audit_run(a.run_dir, a.name, drc=not a.no_drc)
    for c in r['checks']:
        print('%-6s %-9s %-8s %s' % (c['id'], 'ok' if c['ok'] else ('reported' if c['role'] == 'report' else 'FAIL'),
                                     c['role'], c['check']))
    if a.out:
        with open(a.out, 'w') as fh:
            json.dump(r, fh, indent=1, default=float)
        print('audit: %s' % a.out)
    status = 'INCOMPLETE' if (a.no_drc and r['hard_failures'] == ['P9']) else r['status']
    print(AU.summary_line(r) + ('' if status == r['status'] else ' (INCOMPLETE: the DRC was not run)'))
    return {'PASS': 0, 'PASS_OUTSIDE_ENVELOPE': 3, 'FAIL': 2, 'INCOMPLETE': 4}[status]


def _main(a):
    """Run the subcommand that the parsed arguments name and return its exit code."""
    if a.cmd == 'audit':
        return _audit(a)
    if a.cmd == 'generate':
        try:
            kw = _generate_args(a)
        except (ValueError, TypeError) as e:              # a bad argument is a spec error (exit 6) rather than a crash
            e._opencif_args = True
            raise
        from . import api
        r = api.generate(a.cutoff_um, a.filtrate_fraction, out=a.out, render=not a.no_render, drc=not a.no_drc,
                         check_only=a.check, verbose=True, **kw)
        if r.record:
            _print_status(r.record)
        print(r.verdict)
        print('spec: %s' % r.spec_path)
        if r.dxf:
            print('mask: %s' % r.dxf)
            print('record: %s' % r.record_path)
            for f in r.figures:
                print('figure: %s' % f)
        return r.exit_code if r.exit_code is not None else 2
    return _build_or_check(a)


def _print_status(rec):
    """Print the status codes and the source trace in the command-line summary (ADJ-C6)."""
    a = (rec or {}).get('assessment', {}) or {}
    if a.get('status_codes'):
        print('status: ' + ', '.join('%s %s' % (k, v) for k, v in a['status_codes'].items()))
    if a.get('status', {}).get('source_trace'):
        print('source trace: ' + a['status']['source_trace'])


def _generate_args(a):
    """Turn the `generate` arguments into keyword arguments for api.generate, and check them by loading the spec that
    they build."""
    from . import api
    kw = {}
    if a.process:
        kw['process'] = json.loads(a.process) if a.process.strip().startswith('{') else a.process
    if a.chip:
        kw['chip'] = tuple(float(v) for v in a.chip.lower().split('x'))
    for opt in ('layout', 'taper', 'seal', 'inlet', 'outlet', 'start', 'name', 'mode'):
        if getattr(a, opt):
            kw[opt] = getattr(a, opt)
    if a.legs:
        kw['legs'] = a.legs if a.legs == 'auto' else int(a.legs)
    if a.meander:
        kw['meander'] = None if a.meander == 'none' else a.meander
    if a.width_in is not None:
        kw['width_in_um'] = a.width_in
    if a.width_min is not None:
        kw['width_min_um'] = a.width_min
    if a.flow is not None:
        kw['flow_uL_min'] = a.flow
    # The test is against None, so that --array 0 is refused rather than dropped (RC8 code review, n1).
    if a.array is not None:
        kw['array'] = a.array
    ov = {}
    for item in a.set:
        key, _, val = item.partition('=')
        sec, _, k = key.partition('.')
        try:
            import yaml
            v = yaml.safe_load(val)
        except ImportError:
            v = json.loads(val)
        ov.setdefault(sec, {})[k] = v
    if ov:
        kw['overrides'] = ov
    SP.load(api.build_spec(a.cutoff_um, a.filtrate_fraction, **kw))      # argument and schema errors surface here
    return kw


def _build_or_check(a):
    """Build a spec file, or check it without writing a photomask, and return the exit code."""
    if a.cmd == 'build':
        d = SP.load(a.spec)
        if a.no_drc:
            d.setdefault('checks', {})['drc'] = False
        rec = PL.run(d, a.out, render=not a.no_render)
        _print_status(rec)                                # before the verdict, which stays the last line printed
        print(rec['assessment']['verdict'])
        return rec['assessment'].get('exit_code', 2)
    d = SP.load(a.spec)
    from opencif.physics import schedule as S
    S.clear_caches()                                      # as in pipeline.run
    spec, opts, sc, g, env_bad, notes = PL.prepare(d, print)
    from . import build as BU
    # The layout and the outline that a build would draw are made in memory, with nothing written, and the collision
    # checks that a build applies before it writes the photomask are applied to them, so that a check refuses what a
    # build would refuse (the clearance guarantee).
    lay = PL.layout_check(spec, sc, opts, notes, print) if not PL.refusal(opts, env_bad)[0] else None
    print(json.dumps(dict(N=sc.N if sc is not None else None, NM=sc.NM if sc is not None else None,
                          F_design=sc.F if sc is not None else None, sealed=PL.is_sealed(sc) if sc is not None else False,
                          D_um=S.KAPPA_G * spec.delta_star, kernel_pair=BU.KERNEL_PAIR,
                          seal_refusal=opts.get('seal_refusal'), outlet_refusal=opts.get('outlet_refusal'),
                          layout_refusal=opts.get('layout_refusal'), preflight_refusal=opts.get('preflight_refusal'),
                          envelope_failures=env_bad, layout=lay, notes=notes, gates=g),
                     indent=1, default=float))
    # The exit code is the one that a build would give on the planned and outlined design, as in
    # generate(check_only=True), which is 5 for a refused design, 3 for one outside the envelope (evaluate and
    # reproduce modes), and 0 for one inside it. The drawn mask's measurement, network, and DRC are not run.
    return 5 if PL.refusal(opts, env_bad)[0] else (3 if env_bad else 0)


if __name__ == '__main__':
    sys.exit(main())
