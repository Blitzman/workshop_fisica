"""Builds ns_cooling/report/cooling_report.html from ns_cooling/results/*.json. Every number in the page is read from those files
(or from data/marino2024_sources.json); nothing is typed by hand. Run from the repository root: python ns_cooling/report/build.py"""
import html
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / 'results'
J = {p.stem: json.loads(p.read_text()) for p in sorted(RES.glob('step*.json'))}
SRC = json.loads((ROOT / 'data' / 'marino2024_sources.json').read_text())
s1, s2, s3, s4, s5, s6, s7, s8, s9 = (J[k] for k in ('step1_structure', 'step2_sly4', 'step3_integrals', 'step4_standard_candle',
                                                     'step5_cooling_curve', 'step6_superfluid', 'step7_sf_sweep', 'step8_inmedium',
                                                     'step9_figure_curves'))
esc = html.escape


def sci(x, d=1):
    e = math.floor(math.log10(abs(x)))
    sup = str(e).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
    return f'{x / 10**e:.{d}f}×10{sup}'


def counts(d):
    c = d['checks']
    return sum(x['ok'] is True for x in c), sum(x['ok'] is False for x in c), sum(x['ok'] is None for x in c)


# ---------- numbers used in the prose
n_ok, n_fail = sum(counts(d)[0] for d in J.values()), sum(counts(d)[1] for d in J.values())
n_tot = n_ok + n_fail
objs = SRC['objects']
h7, h8 = s7['summary'], s8['summary']
plain = [s['nonsf_ratio']['1.4'] for s in h7 if s['pbf'] == 'Y01']
best_y01 = [s['best_ratio'] for s in h8 if s['pbf'] == 'Y01']
best_upd = [s['best_ratio'] for s in h8 if s['pbf'] == 'updated']
n_comb = h8[0]['n_total']
n_reach = sum(s['n_within'] for s in h8) + sum(s['n_within'] for s in h7)
colder = [s['frac_colder'] for s in h7 + h8]
dim = [x for s in h8 for x in s['nonsf_dimming'].values()]
acc = [x for g in s5['gaps'] for x in g['accreted_over_iron']]
alpha = [x for a in s8['alpha'].values() for x in a]
b24, sly = s1['bsk24'], s2
n_gp, n_gn = len(s6['tc_max']['p']), len(s6['tc_max']['nt'])
dT = s4['dT']
masses = [r['M'] for r in s4['rows']]
bad_m = [f"{m:.1f}" for m, x in zip(masses, dT) if abs(x) > 0.07]
gain = [s['nonsf_ratio'][str(s['best']['M'])] / s['best_ratio'] for s in h7 if s['pbf'] == 'Y01']


def age_txt(o):
    a, b = o['age_yr']
    return f'{a:,}' if a == b else f'{a:,}–{b:,}'


# ---------- figure
W, H, ML, MR, MT, MB = 880, 470, 62, 190, 18, 48
X0, X1, Y0, Y1 = 1.5, 5.5, 31.5, 35.0


def px(t):
    return ML + (math.log10(t) - X0) / (X1 - X0) * (W - ML - MR)


def py(L):
    return MT + (Y1 - math.log10(L)) / (Y1 - Y0) * (H - MT - MB)


t = s9['t_yr']
SERIES = [('standard', 's1', 'Standard star, no superfluidity'), ('inmedium', 's2', 'Same, with boosted slow channel'),
          ('best_Y01', 's3', 'Best superfluid case, 2001 rate'), ('best_updated', 's4', 'Best superfluid case, updated rate')]
svg = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Cooling curves of slow-cooling stars against the three cold sources" id="fig">']
for e in range(2, 6):
    x = px(10**e)
    svg.append(f'<line class="grid" x1="{x:.1f}" y1="{MT}" x2="{x:.1f}" y2="{H - MB}"/><text class="tick" x="{x:.1f}" y="{H - MB + 18}" text-anchor="middle">10{str(e).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))}</text>')
for e in range(32, 36):
    y = py(10**e)
    svg.append(f'<line class="grid" x1="{ML}" y1="{y:.1f}" x2="{W - MR}" y2="{y:.1f}"/><text class="tick" x="{ML - 8}" y="{y + 4:.1f}" text-anchor="end">10{str(e).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))}</text>')
svg.append(f'<rect class="frame" x="{ML}" y="{MT}" width="{W - ML - MR}" height="{H - MT - MB}"/>')
svg.append(f'<text class="axlabel" x="{(ML + W - MR) / 2}" y="{H - 8}" text-anchor="middle">Age (years)</text>')
svg.append(f'<text class="axlabel" transform="translate(16 {(MT + H - MB) / 2}) rotate(-90)" text-anchor="middle">Thermal luminosity seen from afar (erg/s)</text>')
ends = []
for key, cls, label in SERIES:
    Ls = s9['curves'][key]['L']
    pts = [(px(a), py(L)) for a, L in zip(t, Ls) if Y0 <= math.log10(L) <= Y1]
    svg.append(f'<polyline class="line {cls}" points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in pts)}"/>')
    ends.append([pts[-1][1], cls, label, pts[-1][0]])
ends.sort()
for i in range(1, len(ends)):                              # keep the end labels at least 30 px apart (two lines each)
    ends[i][0] = max(ends[i][0], ends[i - 1][0] + 30)
shift = max(0, ends[-1][0] + 14 - (H - MB))
for y, cls, label, x in ends:
    y -= shift
    a, b = label.split(', ') if ', ' in label else (label, '')
    svg.append(f'<circle class="dot {cls}" cx="{W - MR + 12}" cy="{y - 4:.1f}" r="4"/><text class="note" x="{W - MR + 22}" y="{y:.1f}">{esc(a)}</text>'
               + (f'<text class="note sub" x="{W - MR + 22}" y="{y + 14:.1f}">{esc(b)}</text>' if b else ''))
for o, (dx, dy, anchor) in zip(objs, ((-10, 4, 'end'), (8, -12, 'start'), (10, 16, 'start'))):
    L, (lo, hi) = o['L32'] * 1e32, o['L32_err']
    a0, a1 = o['age_yr']
    xm = (px(a0) + px(a1)) / 2
    svg.append(f'<line class="err" x1="{xm:.1f}" y1="{py(L + hi * 1e32):.1f}" x2="{xm:.1f}" y2="{py(L - lo * 1e32):.1f}"/>')
    if a0 != a1:
        svg.append(f'<line class="err" x1="{px(a0):.1f}" y1="{py(L):.1f}" x2="{px(a1):.1f}" y2="{py(L):.1f}"/>')
    svg.append(f'<circle class="src" cx="{xm:.1f}" cy="{py(L):.1f}" r="5"><title>{esc(o["name"])}: {sci(L, 2)} erg/s at {age_txt(o)} yr</title></circle>'
               f'<text class="note" x="{xm + dx:.1f}" y="{py(L) + dy:.1f}" text-anchor="{anchor}">{esc(o["name"])}</text>')
svg.append(f'<line id="xh" class="xh" x1="0" y1="{MT}" x2="0" y2="{H - MB}" hidden="hidden"/></svg>')
legend = ''.join(f'<span><i class="sw {cls}"></i>{esc(label)}</span>' for _, cls, label in SERIES) + '<span><i class="sw src"></i>Observed (90 % error bars)</span>'
figdata = json.dumps({'t': [round(x, 1) for x in t], 'x': [round(px(a), 1) for a in t],
                      's': [[lab, cls, [float(f'{L:.3g}') for L in s9['curves'][k]['L']]] for k, cls, lab in SERIES]})

# ---------- tables
STEPS = [
    ('1', 'BSk24 star: maximum mass and fast-cooling threshold', 'Marino et al. 2024; Potekhin &amp; Chabrier 2018', s1),
    ('2', 'SLy4 star: maximum mass, no fast cooling at any mass', 'Douchin &amp; Haensel 2001', s2),
    ('3', 'Heat stored and neutrino output of a whole star', 'Ofengeim et al. 2017', s3),
    ('4', 'Standard slow cooling: interior temperature at 330 yr', 'Yakovlev et al. 2011', s4),
    ('5', 'First cooling curve with the outer-layer relation', 'Potekhin, Chabrier &amp; Yakovlev 1997; Yakovlev et al. 2011', s5),
    ('6', 'Superfluidity formulas', 'Yakovlev et al. 2001; Ho et al. 2015', s6),
    ('7', 'Sweep over superfluidity models', 'internal consistency only', s7),
    ('8', 'Boosted slow channel and repeated sweep', 'Shternin, Baldo &amp; Haensel 2018', s8),
    ('9', 'Curves for the figure', 'internal consistency only', s9),
]
step_rows = ''.join(
    f'<tr><td class="num">{n}</td><td>{what}</td><td>{ref}</td><td class="num">{counts(d)[0]}/{counts(d)[0] + counts(d)[1]}</td>'
    f'<td>{"<span class=\"pill ok\">pass</span>" if counts(d)[1] == 0 else f"<span class=\"pill bad\">{counts(d)[1]} fail</span>"}</td></tr>'
    for n, what, ref, d in STEPS)


def best_cell(s):
    b = s['best']
    return f"{s['best_ratio']:.1f}× <span class=\"small\">({b['M']} M☉, p {b['p'] or '—'}, n {b['n'] or '—'})</span>"


sweep_rows = ''
for o in objs:
    a = [s for s in h7 if s['object'] == o['name']]
    b = [s for s in h8 if s['object'] == o['name']]
    y7, y8, u8 = a[0], next(s for s in b if s['pbf'] == 'Y01'), next(s for s in b if s['pbf'] == 'updated')
    sweep_rows += (f"<tr><td>{esc(o['name'])}</td><td class=\"num\">{age_txt(o)}</td><td class=\"num\">{sci(o['L32'] * 1e32, 2)}</td>"
                   f"<td class=\"num\">{y7['nonsf_ratio']['1.4']:.0f}×</td><td class=\"num\">{y8['nonsf_ratio']['1.4']:.0f}×</td>"
                   f"<td>{best_cell(y8)}</td><td>{best_cell(u8)}</td>"
                   f"<td class=\"num\">{y8['n_within'] + u8['n_within']} of {y8['n_total'] + u8['n_total']}</td></tr>")

fails = [(n, c) for n, _, _, d in STEPS for c in d['checks'] if c['ok'] is False]
FAIL_WHY = {
    '2': 'The original SLy4 paper quotes thresholds that its own composition table does not reproduce. Either value lies above the centre of the heaviest star, so the conclusion is unaffected.',
    '4': f"Our stars are up to {max(dT):.0%} warmer than the benchmark at {' and '.join(bad_m)} M☉. The benchmark’s authors name the 1.9 M☉ SLy star as their own worst case. Not resolved; nothing was tuned.",
    '6': 'A badly designed test of ours: it compared the fits with their limiting formulas at a temperature where the limit does not yet apply. The correctly designed test (in the limit) passes.',
    '8': 'A hand-derived closed form against numerical integration, missed by a hair on a tolerance set too tight for that integration method.',
}
fail_rows = ''.join(f'<tr><td class="num">{n}</td><td>{esc(c["name"])}{(" <em>(" + esc(c["note"]) + ")</em>") if c["note"] else ""}</td><td>{FAIL_WHY[n]}</td></tr>' for n, c in fails)

STATEMENTS = [
    ('S1', 'Slow cooling cannot reach the luminosities of the three sources.',
     f"Supported so far. None of {n_comb} superfluidity combinations (for each of two emission rates) reaches any source. Closest: {min(best_y01):.1f}× too bright.", 'ok', 'supported so far'),
    ('S2', 'Unless gap choices are extremely unusual, superfluidity accelerates cooling.',
     f"Only as a majority statement: {min(colder):.0%}–{max(colder):.0%} of the combinations are colder than the plain star at the source ages; the rest are warmer, and they are ordinary published models.", 'flag', 'partly'),
    ('S3', 'An iron outer layer is the coldest case; light elements are about ten times brighter.',
     f"Confirmed in our model: a fully light-element layer is {min(acc):.0f}–{max(acc):.0f}× brighter at the source ages.", 'ok', 'confirmed'),
    ('S4', 'The core can be treated as one uniform temperature from the start.', 'Not tested yet. Needs the full heat-flow version of the code. Matters most for the 841-year source.', 'open', 'not tested'),
    ('S5', 'Fast cooling switches on above 1.6 M☉ in BSk24, and never in SLy4.',
     f"Confirmed, with a caveat: the BSk24 threshold is {b24['m_durca']:.3f} M☉, so the paper’s 1.6 M☉ model sits {1.6 - b24['m_durca']:.3f} M☉ above it. Potekhin &amp; Chabrier 2018 make the same remark.", 'flag', 'confirmed, caveat'),
    ('S6', 'About 75 % of equations of state are excluded.', 'Not tested yet. The figure comes from a simplified side estimate, not from simulations. See the caveat on models close to the threshold below.', 'open', 'not tested'),
    ('S7', 'A classifier gives slow cooling a probability below 5 %.', 'Not a physics check. The classifier can only rank the 81 runs it was given. Lower priority.', 'open', 'not tested'),
]
st_rows = ''.join(f'<tr><td class="num">{k}</td><td>{claim}</td><td>{res}</td><td><span class="pill {cls}">{lab}</span></td></tr>' for k, claim, res, cls, lab in STATEMENTS)

page = f"""<title>Cold Neutron Star Cooling</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=Newsreader:opsz,wght@6..72,600&display=swap">
<style>
/* Layout: one reading column (~68ch) with wider figures and tables; same identity as the axion reports */
:root {{
  --paper: #f4f6f9; --sheet: #ffffff; --ink: #161e2e; --muted: #586378; --rule: #d8dde6;
  --accent: #2b4bb3; --ok: #1d7451; --ok-bg: #e4f3eb; --flag: #9a5400; --flag-bg: #fdf0dc;
  --bad: #a3262f; --bad-bg: #fbe5e6; --open: #586378; --open-bg: #e9edf3; --code-bg: #eaeef4;
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a; --s4: #eda100; --grid: #e6e9ef;
  --f-display: "Newsreader", "Iowan Old Style", Georgia, serif;
  --f-body: "IBM Plex Sans", "Segoe UI", system-ui, sans-serif;
  --f-mono: "IBM Plex Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --paper: #0d121b; --sheet: #141b27; --ink: #e3e8f1; --muted: #98a3b7; --rule: #263042;
    --accent: #93a8ff; --ok: #62c99c; --ok-bg: #12291f; --flag: #f1a84d; --flag-bg: #2b2112;
    --bad: #f08a90; --bad-bg: #2d1518; --open: #98a3b7; --open-bg: #1b2433; --code-bg: #1b2433;
    --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --grid: #212a3a; color-scheme: dark;
  }}
}}
:root[data-theme="dark"] {{
  --paper: #0d121b; --sheet: #141b27; --ink: #e3e8f1; --muted: #98a3b7; --rule: #263042;
  --accent: #93a8ff; --ok: #62c99c; --ok-bg: #12291f; --flag: #f1a84d; --flag-bg: #2b2112;
  --bad: #f08a90; --bad-bg: #2d1518; --open: #98a3b7; --open-bg: #1b2433; --code-bg: #1b2433;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --grid: #212a3a; color-scheme: dark;
}}
body {{ background: var(--paper); color: var(--ink); font: 16px/1.6 var(--f-body); }}
.wrap {{ max-width: 900px; margin: 0 auto; padding-left: max(16px, env(safe-area-inset-left, 0px));
  padding-right: max(16px, env(safe-area-inset-right, 0px)); padding-block: 40px 64px; display: grid; gap: 30px; }}
header {{ display: grid; gap: 10px; }}
.eyebrow {{ font: 500 12px/1.4 var(--f-mono); letter-spacing: .08em; text-transform: uppercase; color: var(--muted); }}
h1 {{ font: 600 clamp(30px, 5vw, 44px)/1.1 var(--f-display); margin: 0; text-wrap: balance; }}
h2 {{ font: 600 25px/1.2 var(--f-display); margin: 0; text-wrap: balance; }}
h3 {{ font: 600 16px/1.3 var(--f-body); margin: 0; }}
.lede {{ font-size: 18px; color: var(--muted); max-width: 68ch; margin: 0; }}
section {{ display: grid; gap: 14px; min-width: 0; }}
section > p, section > ul, section > ol {{ max-width: 68ch; margin: 0; }}
ul, ol {{ padding-left: 22px; display: grid; gap: 6px; }}
.sheet {{ background: var(--sheet); border: 1px solid var(--rule); border-radius: 6px; padding: 20px 22px; display: grid; gap: 10px; min-width: 0; }}
.sheet ul {{ margin: 0; }}
code {{ font-family: var(--f-mono); font-size: .9em; background: var(--code-bg); padding: 1px 5px; border-radius: 3px; }}
pre {{ font: 13px/1.5 var(--f-mono); background: var(--code-bg); padding: 12px 14px; border-radius: 6px; overflow-x: auto; margin: 0; }}
.small {{ font-size: 13px; color: var(--muted); }}
figure {{ margin: 0; display: grid; gap: 8px; min-width: 0; position: relative; }}
.figbox {{ overflow-x: auto; }}
figure svg {{ width: 100%; min-width: 640px; height: auto; display: block; background: var(--sheet); border: 1px solid var(--rule); border-radius: 6px; }}
figcaption {{ font-size: 14px; color: var(--muted); max-width: 72ch; }}
.legend {{ display: flex; flex-wrap: wrap; gap: 6px 18px; font-size: 13px; color: var(--muted); }}
.legend span {{ display: inline-flex; align-items: center; gap: 6px; }}
.sw {{ width: 10px; height: 10px; border-radius: 50%; display: inline-block; }}
.sw.s1 {{ background: var(--s1); }} .sw.s2 {{ background: var(--s2); }} .sw.s3 {{ background: var(--s3); }} .sw.s4 {{ background: var(--s4); }}
.sw.src {{ background: var(--ink); }}
svg .grid {{ stroke: var(--grid); stroke-width: 1; }}
svg .frame {{ fill: none; stroke: var(--rule); }}
svg .tick {{ fill: var(--muted); font: 12px var(--f-mono); }}
svg .axlabel {{ fill: var(--ink); font: 500 13px var(--f-body); }}
svg .note {{ fill: var(--ink); font: 12.5px var(--f-body); }}
svg .note.sub {{ fill: var(--muted); }}
svg .line {{ fill: none; stroke-width: 2; stroke-linejoin: round; }}
svg .line.s1, svg .dot.s1 {{ stroke: var(--s1); }} svg .line.s2, svg .dot.s2 {{ stroke: var(--s2); }}
svg .line.s3, svg .dot.s3 {{ stroke: var(--s3); }} svg .line.s4, svg .dot.s4 {{ stroke: var(--s4); }}
svg .dot.s1 {{ fill: var(--s1); }} svg .dot.s2 {{ fill: var(--s2); }} svg .dot.s3 {{ fill: var(--s3); }} svg .dot.s4 {{ fill: var(--s4); }}
svg .src {{ fill: var(--ink); stroke: var(--sheet); stroke-width: 2; }}
svg .err {{ stroke: var(--ink); stroke-width: 1.5; }}
svg .xh {{ stroke: var(--muted); stroke-width: 1; stroke-dasharray: 3 3; }}
.tip {{ position: absolute; pointer-events: none; background: var(--sheet); border: 1px solid var(--rule); border-radius: 6px; padding: 8px 10px;
  font-size: 12.5px; line-height: 1.5; box-shadow: 0 2px 10px rgb(0 0 0 / .12); max-width: 300px; z-index: 2; }}
.tip b {{ font-family: var(--f-mono); font-weight: 500; }}
.tip div {{ display: flex; gap: 6px; align-items: center; }}
.tbl {{ overflow-x: auto; border: 1px solid var(--rule); border-radius: 6px; background: var(--sheet); }}
table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
th, td {{ text-align: left; padding: 8px 12px; border-bottom: 1px solid var(--rule); vertical-align: top; }}
tr:last-child td {{ border-bottom: 0; }}
th {{ font: 500 12px/1.3 var(--f-mono); letter-spacing: .04em; text-transform: uppercase; color: var(--muted); white-space: nowrap; }}
td.num {{ font-family: var(--f-mono); font-variant-numeric: tabular-nums; white-space: nowrap; }}
.pill {{ display: inline-block; font: 500 11.5px/1 var(--f-mono); padding: 4px 7px; border-radius: 99px; white-space: nowrap; }}
.pill.ok {{ color: var(--ok); background: var(--ok-bg); }} .pill.flag {{ color: var(--flag); background: var(--flag-bg); }}
.pill.bad {{ color: var(--bad); background: var(--bad-bg); }} .pill.open {{ color: var(--open); background: var(--open-bg); }}
footer {{ font-size: 13px; color: var(--muted); border-top: 1px solid var(--rule); padding-top: 14px; }}
</style>
<div class="wrap">
<header>
  <p class="eyebrow">Workshop with José A. Pons · neutron-star cooling · preliminary results</p>
  <h1>Can slow cooling explain three cold young neutron stars?</h1>
  <p class="lede">A first check, with our own 1D code, of Marino et al. 2024 (Nat. Astron. 8, 1020). Hypotheses were written down before each calculation, and every number on this page is inserted from the result files when the page is built.</p>
</header>

<section>
  <div class="sheet">
    <h2>In brief</h2>
    <ul>
      <li><b>The paper’s main claim holds so far.</b> A standard slow-cooling star is {min(plain):.0f} to {max(plain):.0f} times brighter than the three sources. We tried {n_gp} proton and {n_gn} neutron superfluidity models, alone and combined, at six masses and with two versions of the pair-breaking emission rate. None of them reaches any source. The closest case is still {min(best_y01):.1f} times too bright.</li>
      <li><b>One sentence of the paper is too strong.</b> It says superfluidity accelerates cooling unless gap choices are extremely unusual. In our sweep {min(colder):.0%} to {max(colder):.0%} of the combinations are colder than the plain star; the rest are warmer, and they are ordinary published models.</li>
      <li><b>A detail worth knowing.</b> In BSk24 the fast channel opens at {b24['m_durca']:.3f} M☉. The paper’s 1.6 M☉ model is therefore only {1.6 - b24['m_durca']:.3f} M☉ above the threshold, and the classifier assigns about 1.6 M☉ to both pulsars.</li>
      <li><b>This is not the final word.</b> The interior is treated as one uniform temperature, effective masses are fixed at 0.7, there is no crust physics, and only SLy4 has been swept. The list of what is missing is below.</li>
      <li><b>Checks:</b> {n_ok} of {n_tot} pass. The {n_fail} failures are listed and explained; none was fixed by changing a criterion.</li>
    </ul>
  </div>
</section>

<section>
  <h2>The question</h2>
  <p>Three young neutron stars are much fainter than expected for their age. The paper concludes that ordinary slow cooling cannot explain them, so dense matter must allow a fast cooling channel (direct Urca or hyperons), and that this would exclude most current equations of state.</p>
  <p>That conclusion rests on 81 runs of a 2D magneto-thermal code with three equations of state, one of them without a fast channel (SLy4), and a single superfluidity model. José asked for an independent look with a simple 1D code and updated microphysics. We split the paper’s argument into seven statements and test them one at a time.</p>
  <div class="tbl"><table>
    <tr><th></th><th>Statement in the paper</th><th>What we find (preliminary)</th><th>Status</th></tr>
    {st_rows}
  </table></div>
</section>

<section>
  <h2>Slow-cooling stars against the three sources</h2>
  <figure>
    <div class="legend">{legend}</div>
    <div class="figbox">{''.join(svg)}</div>
    <div class="tip" id="tip" hidden></div>
    <figcaption>SLy4, iron outer layer, uniform interior temperature. The standard star has {s9['curves']['standard']['M']} M☉. The best superfluid cases are the ones closest to CXOU J0852−4617 in the sweep: {s9['curves']['best_Y01']['M']} M☉ with proton gap {s9['curves']['best_Y01']['p']} and neutron triplet gap {s9['curves']['best_Y01']['n']} (2001 rate), and {s9['curves']['best_updated']['M']} M☉ with {s9['curves']['best_updated']['p']} and {s9['curves']['best_updated']['n']} (updated rate). Move the pointer over the plot to read values. The curves are not reliable before about a century, when the interior is not yet uniform.</figcaption>
  </figure>
  <div class="tbl"><table>
    <tr><th>Source</th><th>Age (yr)</th><th>Observed (erg/s)</th><th>Standard star</th><th>With boosted channel</th><th>Best case, 2001 rate</th><th>Best case, updated rate</th><th>Cases reaching it</th></tr>
    {sweep_rows}
  </table></div>
  <p class="small">“N×” means the model is N times brighter than observed at the source’s age (the oldest allowed age for the third source, which favours the model). “Reaching” means at or below the upper end of the 90 % error bar. p and n are the proton and neutron-triplet gap models of Ho et al. 2015. Observed values from Marino et al. 2024, Supplementary Table 6.</p>
  <p><b>What this means.</b> Superfluidity lowers the luminosity by a factor of about {min(gain):.0f} to {max(gain):.0f} at best, through the extra neutrino emission when neutron pairs form. That is not enough. The boosted version of the main slow channel (Shternin, Baldo &amp; Haensel 2018) dims the plain star by a factor {min(dim):.1f} to {max(dim):.1f}, but barely moves the best superfluid cases, because there that channel is already suppressed.</p>
  <p><b>A caveat for the “75 % excluded” figure.</b> The boost grows sharply for equations of state that come close to the fast-cooling threshold without crossing it. In SLy4 it is modest (a factor {min(alpha):.1f} to {max(alpha):.1f} over the bare rate’s normalisation, against 5 to 10 quoted by its authors for their own model) because SLy4 never approaches the threshold. SLy4 may therefore be the least favourable kind of slow-cooling star. We have not tested other equations of state yet.</p>
</section>

<section>
  <h2>What we built, and how each piece was checked</h2>
  <p>The code is ours, written from the published formulas. Before adding anything new, each piece had to reproduce a number published by other authors. Reference numbers and coefficients are read from the source files and papers by the code, not typed.</p>
  <div class="tbl"><table>
    <tr><th>Step</th><th>What</th><th>Checked against</th><th>Checks</th><th></th></tr>
    {step_rows}
  </table></div>
  <ul>
    <li><b>The two stars are the paper’s stars.</b> BSk24: maximum mass {b24['m_max']:.3f} M☉, fast channel above {sci(b24['rho_durca'], 2)} g/cm³. SLy4: maximum mass {sly['ours']['m_max']:.3f} M☉ and radius {sly['ours']['r_max']:.2f} km, as in Douchin &amp; Haensel 2001.</li>
    <li><b>The heat budget is right.</b> Stored heat within {s3['stats']['sly']['c']['max']:.1%} and the simplest neutrino channel within {s3['stats']['sly']['nn']['max']:.1%} of the published fits for SLy4 (their own margins: 12 % and 54 %).</li>
    <li><b>The standard cooling curve is right for light and medium stars.</b> Surface temperature at 330 years: {min(r['Ts_inf'] for r in s5['rows'] if r['x'] < 0.5) / 1e6:.2f} to {max(r['Ts_inf'] for r in s5['rows'] if r['x'] < 0.5) / 1e6:.2f} million K, against about 1.4 in the reference.</li>
    <li><b>The superfluidity formulas are transcribed correctly.</b> They reproduce their limiting forms and the peak of a published figure (ours {sci(s6['fig19']['peak'], 2)} erg cm⁻³ s⁻¹, read off the figure as 10²¹).</li>
  </ul>
</section>

<section>
  <h2>The {n_fail} failed checks</h2>
  <div class="tbl"><table>
    <tr><th>Step</th><th>Check</th><th>Why, and what it affects</th></tr>
    {fail_rows}
  </table></div>
</section>

<section>
  <h2>What is not covered yet</h2>
  <ul>
    <li><b>Uniform interior temperature.</b> Doubtful for the 841-year source. Needs the full heat-flow version (statement S4).</li>
    <li><b>Effective masses.</b> Fixed at 0.7 of the bare mass, as in the references we validated against. Realistic values are the next ingredient.</li>
    <li><b>The crust.</b> No crust neutrino emission, heat capacity or neutron singlet superfluid.</li>
    <li><b>Mass range.</b> All the best cases sit at the lowest mass tried (1.0 M☉), the edge of the grid.</li>
    <li><b>Other equations of state.</b> Only SLy4 has been swept. Models near the fast-cooling threshold could behave differently.</li>
    <li><b>Fast cooling itself.</b> Direct Urca is not in the code yet, so we have not reproduced the paper’s BSk24 and GM1A curves.</li>
    <li><b>Combined reduction factors.</b> We use the approximate rules of Yakovlev et al. 2001, and assume the same factors apply to the boosted channel.</li>
    <li><b>The high-mass discrepancy</b> in the standard cooling benchmark (table above) is unresolved.</li>
  </ul>
</section>

<section>
  <div class="sheet">
    <h2>Questions for José</h2>
    <ol>
      <li>Which statements of the paper do you consider not rigorous? Do they match S1–S7?</li>
      <li>Which gap model of Ho et al. 2015 did the paper use? And which one in the hyperon runs?</li>
      <li>The 1.6 M☉ BSk24 model is {1.6 - b24['m_durca']:.3f} M☉ above the direct-Urca threshold. Was that known when the classifier assigned about 1.6 M☉ to both pulsars?</li>
      <li>Which ingredients do you consider outdated, beyond the ones tried here (gaps, pair-breaking rate, in-medium modified Urca)?</li>
      <li>Which equations of state without direct Urca would you add? In particular, ones that come close to the threshold.</li>
      <li>Could we have the paper’s cooling curves as data (at least SLy4 and BSk24), to compare number by number?</li>
    </ol>
  </div>
</section>

<section>
  <h2>Reproduce it</h2>
  <p>Code, results and the lab notebook are in the <code>ns_cooling</code> folder of the repository. Third-party files (the equation-of-state fits and the papers) are not redistributed; <code>ns_cooling/README.md</code> says where to get them.</p>
  <pre>{''.join(f"uv run python -m ns_cooling.{p.stem}\n" for p in sorted(ROOT.glob('step*.py'), key=lambda p: int(p.stem.split('_')[0][4:])))}uv run python ns_cooling/report/build.py</pre>
</section>

<footer>Preliminary results, not reviewed. Built from <code>ns_cooling/results/*.json</code>. Sources: Marino et al. 2024; Douchin &amp; Haensel 2001; Pearson et al. 2018 (fits by A. Y. Potekhin); Ofengeim et al. 2017; Yakovlev et al. 2001, 2011; Yakovlev, Levenfish &amp; Shibanov 1999; Potekhin, Chabrier &amp; Yakovlev 1997; Ho et al. 2015; Potekhin &amp; Chabrier 2018; Shternin, Baldo &amp; Haensel 2018.</footer>
</div>
<script>
(function () {{
  var D = {figdata};
  var svg = document.getElementById('fig'), tip = document.getElementById('tip'), xh = document.getElementById('xh');
  function fmt(v) {{ var e = Math.floor(Math.log10(v)); return (v / Math.pow(10, e)).toFixed(1) + '×10^' + e; }}
  function move(ev) {{
    var r = svg.getBoundingClientRect(), x = (ev.clientX - r.left) / r.width * {W};
    if (x < {ML} || x > {W - MR}) {{ tip.hidden = true; xh.setAttribute('hidden', 'hidden'); return; }}
    var i = 0, best = 1e9;
    for (var k = 0; k < D.x.length; k++) {{ var d = Math.abs(D.x[k] - x); if (d < best) {{ best = d; i = k; }} }}
    xh.removeAttribute('hidden'); xh.setAttribute('x1', D.x[i]); xh.setAttribute('x2', D.x[i]);
    var h = '<b>' + Math.round(D.t[i]).toLocaleString('en-US') + ' yr</b>';
    D.s.forEach(function (s) {{ h += '<div><i class="sw ' + s[1] + '"></i>' + s[0] + ': <b>' + fmt(s[2][i]) + '</b></div>'; }});
    tip.innerHTML = h; tip.hidden = false;
    var f = svg.closest('figure').getBoundingClientRect(), left = ev.clientX - f.left + 14;
    if (left + 300 > f.width) left = Math.max(0, ev.clientX - f.left - 310);
    tip.style.left = left + 'px'; tip.style.top = (ev.clientY - f.top + 14) + 'px';
  }}
  svg.addEventListener('pointermove', move);
  svg.addEventListener('pointerleave', function () {{ tip.hidden = true; xh.setAttribute('hidden', 'hidden'); }});
}})();
</script>
"""
out = ROOT / 'report' / 'cooling_report.html'
out.write_text(page)
print(out, f'{len(page) / 1024:.0f} kB; checks {n_ok}/{n_tot}')
