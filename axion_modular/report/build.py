"""report: one self-contained HTML page from the module manifests.   python -m report.build [--out FILE]

Reads : <module>/outputs/manifest.json of every module that has been run.
Writes: report/index.html   status per module (checks, runtime, artifacts with links and hashes), the dependency graph,
                            and a STALE flag when an upstream artifact recorded as input has changed since the run.

This replaces, for the modules refactored so far, the bespoke builders of the legacy report/ folder (build.py,
build_radial.py, ...), which are hard-wired to the outputs of the steps not yet ported (N3, N4, A). New modules plug in by
adding one line to MODULES.
"""
import argparse
import html
import json
from pathlib import Path

from common.artifacts import fingerprint

ROOT = Path(__file__).resolve().parent.parent
OUT_DEFAULT = Path(__file__).resolve().parent / 'index.html'

# module -> (title, upstream modules)
MODULES = {
    'R_paper_refs': ('Paper equations and tables as data', []),
    'D0_derive_action': ('Action -> field and background equations', []),
    'D1_derive_radial': ('Radial linearization and toy model (Paper I)', ['D0_derive_action', 'R_paper_refs']),
    'D2_derive_nonradial': ('Cowling and beyond-Cowling equations (Paper II)', ['D0_derive_action', 'R_paper_refs']),
    'S1_star_background': ('BSk26 star with condensate: profile tables', []),
    'N1_modes_cowling': ('Cowling modes and gravitational-atom levels', ['S1_star_background', 'R_paper_refs']),
    'N2_modes_fullgr_noaxion': ('Polytrope modes without axion vs the literature (validation of the extended system)', ['D2_derive_nonradial', 'R_paper_refs']),
    'T_toy_model': ('Toy-model numerics (Paper I Fig. 3, branches)', ['D1_derive_radial']),
}


def _find_input(rec):
    """Locate a recorded input: path relative to the project root, else by file name in the module outputs."""
    p = Path(rec['path'])
    for cand in ([p] if p.is_absolute() else [ROOT / p]):
        if cand.exists():
            return cand
    for m in MODULES:
        cand = ROOT / m / 'outputs' / rec['name']
        if cand.exists():
            return cand
    return None


def load(mod):
    p = ROOT / mod / 'outputs' / 'manifest.json'
    return json.loads(p.read_text(encoding='utf-8')) if p.exists() else None


def stale_inputs(man):
    bad = []
    for rec in man.get('inputs', []):
        f = _find_input(rec)
        if f is None:
            bad.append((rec['name'], 'missing'))
        elif fingerprint(f) != rec['sha256']:
            bad.append((rec['name'], 'changed'))
    return bad


def _row(mod, title, ups):
    man = load(mod)
    if man is None:
        return f'<tr class="none"><td><b>{mod}</b><br><small>{html.escape(title)}</small></td><td colspan="5">not run</td></tr>', None
    s = man.get('summary', {})
    stale = stale_inputs(man)
    cls = 'fail' if s.get('n_fail') else ('stale' if stale else 'ok')
    status = f"{s.get('n_ok', 0)} ok / {s.get('n_fail', 0)} fail / {s.get('n_descriptive', 0)} descr."
    arts = ''.join(f'<li><a href="../{mod}/outputs/{html.escape(a["path"])}">{html.escape(a["path"])}</a> '
                   f'<small>{a["bytes"] // 1024} kB, {a["sha256"][:10]}</small></li>' for a in man.get('artifacts', []))
    stale_txt = ('<br><b>STALE:</b> ' + ', '.join(f'{n} ({why})' for n, why in stale)) if stale else ''
    prod = man.get('producer', {})
    row = (f'<tr class="{cls}"><td><b>{mod}</b><br><small>{html.escape(title)}</small></td>'
           f'<td>{status}{stale_txt}</td><td>{man.get("seconds")} s</td>'
           f'<td><small>{", ".join(ups) or "-"}</small></td>'
           f'<td><small>{prod.get("created", "")}<br>git {html.escape(str((prod.get("git") or {}).get("commit")))}</small></td>'
           f'<td><details><summary>{len(man.get("artifacts", []))} artifacts</summary><ul>{arts}</ul></details></td></tr>')
    return row, man


def checks_section(mod, man):
    items = []
    for c in man.get('checks', []):
        tag = 'RES' if c['ok'] is None else ('OK' if c['ok'] else 'FAIL')
        note = f' <small>-- {html.escape(c["note"])}</small>' if c.get('note') else ''
        items.append(f'<li class="{tag}"><code>{tag}</code> {html.escape(c["name"])}{note}</li>')
    return f'<details><summary><b>{mod}</b>: {len(items)} checks</summary><ul class="checks">{"".join(items)}</ul></details>'


def build(out=OUT_DEFAULT):
    rows, sections = [], []
    for mod, (title, ups) in MODULES.items():
        row, man = _row(mod, title, ups)
        rows.append(row)
        if man:
            sections.append(checks_section(mod, man))
    graph = ' &nbsp; '.join(f'{u} &rarr; <b>{m}</b>' for m, (_, ups) in MODULES.items() for u in ups)
    page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><title>axion_modular status</title>
<style>
body{{font:14px/1.45 system-ui,sans-serif;margin:2rem;max-width:1200px;color:#1f2328;background:#fff}}
@media (prefers-color-scheme: dark){{body{{color:#e6edf3;background:#0d1117}}a{{color:#79c0ff}}}}
table{{border-collapse:collapse;width:100%}}td,th{{border:1px solid #8883;padding:6px 8px;vertical-align:top;text-align:left}}
tr.ok td:nth-child(2){{border-left:5px solid #2da44e}}tr.fail td:nth-child(2){{border-left:5px solid #cf222e}}
tr.stale td:nth-child(2){{border-left:5px solid #d4a72c}}tr.none{{opacity:.55}}
ul.checks{{list-style:none;padding-left:.5rem}}li.FAIL{{color:#cf222e}}li.RES{{opacity:.8}}code{{font-size:12px}}
</style></head><body>
<h1>axion_modular: module status</h1>
<p>Generated from the manifests in <code>&lt;module&gt;/outputs/manifest.json</code>. Green: all checks pass;
red: at least one check fails; amber: an input artifact changed after the run (re-run the module).
Modules N3, N4 and A are not part of this layout yet.</p>
<p><small>Dependencies: {graph}</small></p>
<table><tr><th>module</th><th>checks</th><th>time</th><th>reads</th><th>produced</th><th>artifacts</th></tr>
{''.join(rows)}</table>
<h2>Checks</h2>
{''.join(sections)}
</body></html>'''
    Path(out).write_text(page, encoding='utf-8')
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(OUT_DEFAULT))
    print('wrote', build(ap.parse_args().out))
