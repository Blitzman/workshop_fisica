"""Genera informe_axiones_completo.html (en report/ y results/): el informe principal en un solo fichero, sin
dependencias externas para el contenido. Se ejecuta después de report/build.py.

- Las fórmulas se convierten a SVG al generar (report/mathjax_svg.cjs, MathJax 3.2.2 en Node), así que la página no
  carga MathJax ni necesita conexión. Requisito: node y `npm ci --prefix report`.
- El código y los datos van incrustados (mismo formato que el informe radial completo: el extractor sirve para los
  dos). Con ellos, `uv run python report/build.py` regenera results/informe_axiones.html idéntico.
- Los enlaces a ficheros de GitHub apuntan a su copia incrustada.
Las tipografías de Google Fonts siguen enlazadas; sin conexión se usan las del sistema.
"""
import html
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GH = 'https://github.com/Blitzman/workshop_fisica/blob/main/'

STEPS = ['step1_background', 'step2_perturbations', 'step3_toy_model', 'step4_numerics', 'step5_nonradial',
         'step6_nonradial_numerics', 'step7_fullgr', 'step8_numerics', 'step9_axion_cowling', 'step10_axion_gr',
         'step14_radial_branches', 'step15_radial_star', 'step16_radial_stability', 'step17_radial_validation',
         'step18_rydberg', 'step19_avoided_crossing', 'step20_damped_count', 'step21_deep_artifacts',
         'step22_exterior_contour', 'step23_eos_consistency', 'step24_scattering_width',
         'step25_exterior_convergence']
CODE = ([f'axion_ns/{m}.py' for m in ('__init__', 'geometry', 'model', 'paper', 'compare', 'pipeline', 'nonradial',
                                       'fullgr', 'fmodes', 'bsk26', 'axstar', 'axgr', 'radial', 'figdata', 'run_all')]
        + [f'axion_ns/{s}.py' for s in STEPS]
        + ['report/build.py', 'report/template.html', 'report/charts.js', 'report/anim.js',
           'report/guia_departamento.html', 'report/build_axiones_standalone.py', 'report/mathjax_svg.cjs', 'report/package.json',
           'report/package-lock.json',
           'pyproject.toml', 'uv.lock'])
DATA = (['results/checks.json'] + [f'results/fig_{n}.json' for n in ('toy', 'levels', 'modes', 'meff', 'background',
                                                                     'axcow', 'axgr')]
        + [f'report/baseline_pre_h29/{n}' for n in ('checks.json', 'fig_axcow.json', 'fig_axgr.json')]
        + sorted(str(p.relative_to(ROOT)) for p in (ROOT / 'results').glob('radial_*.json'))
        + ['results/exterior_convergence.json'])
EXTRACT = """import html, re, sys
from pathlib import Path
src, dst = Path(sys.argv[1]).read_text(), Path(sys.argv[2])
for path, body in re.findall(r'<pre class="code" data-path="([^"]+)"><code>(.*?)</code></pre>', src, flags=re.S):
    f = dst / html.unescape(path)
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(html.unescape(body))
    print(f)"""


def anchor(path):
    return 'f-' + re.sub(r'[^A-Za-z0-9_.-]', '-', path)


def main():
    page = (ROOT / 'results' / 'informe_axiones.html').read_text()
    # 1. sin MathJax en el navegador: fórmulas a SVG al generar
    mj = re.compile(r'<script>\s*window\.MathJax = .*?</script>\s*<script src="[^"]*mathjax[^"]*"[^>]*></script>\n?', re.S)
    assert len(mj.findall(page)) == 1, 'no encuentro el bloque de MathJax'
    page = mj.sub('', page)
    # la página publicada es un fragmento (la plataforma añade el esqueleto): cabecera hasta <main>, cuerpo desde ahí
    head, body = page.split('<main class="wrap">', 1)
    r = subprocess.run(['node', str(ROOT / 'report' / 'mathjax_svg.cjs')], input='<main class="wrap">' + body,
                       capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f'mathjax_svg.cjs falló: {r.stderr.strip()}')
    mjx_style, body = r.stdout.split('\n<!--MJX-->\n', 1)
    n_math = r.stderr.strip()
    page = ('<!doctype html>\n<html lang="es">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
            f'{head.strip()}\n{mjx_style}\n</head>\n<body>\n{body}\n</body>\n</html>\n')
    assert '\\(' not in re.sub(r'title="[^"]*"', '', page), 'quedan fórmulas sin convertir'
    # 2. enlaces a GitHub -> copia incrustada
    files = CODE + DATA
    missing = [f for f in files if not (ROOT / f).exists()]
    assert not missing, missing
    linked = set(re.findall(re.escape(GH) + r'([^"]+)"', page))
    assert linked <= set(files), sorted(linked - set(files))
    page = re.sub(re.escape(GH) + r'([^"]+)" target="_blank" rel="noopener"', lambda m: f'#{anchor(m.group(1))}"', page)
    page = page.replace('href="guia_departamento.html"', f'href="#{anchor("report/guia_departamento.html")}"')
    # 3. ficheros incrustados, al final de «Reproducirlo»
    blocks = []
    for f in files:
        text = (ROOT / f).read_text()
        blocks.append(f'<details class="src" id="{anchor(f)}"><summary><code>{html.escape(f)}</code> '
                      f'<span class="small">({len(text.splitlines())} líneas)</span></summary>'
                      f'<pre class="code" data-path="{html.escape(f)}"><code>{html.escape(text)}</code></pre></details>')
    embedded = (
        '\n    <h3>Todo el código y los datos, dentro de esta página</h3>\n'
        '    <p>Esta copia del informe lleva incrustados todos los ficheros necesarios para repetir el cálculo y regenerar la página, '
        'sin acceso a GitHub. Para extraerlos a una carpeta (con Python 3, sin dependencias), guarda esto como <code>extraer.py</code> '
        'y ejecuta <code>python3 extraer.py informe_axiones_completo.html carpeta</code>. Dentro de <code>carpeta</code> funcionan '
        'los pasos de arriba desde <code>uv sync --frozen</code>. El paso 4 regenera <code>results/informe_axiones.html</code> idéntico al publicado.</p>\n'
        f'    <pre class="code"><code>{html.escape(EXTRACT)}</code></pre>\n'
        f'    <div class="srcs">{"".join(blocks)}</div>\n')
    marker = '<section class="stack repro" id="reproducirlo">'
    i = page.index(marker)
    j = page.index('</section>', i)
    page = page[:j] + embedded + '  ' + page[j:]
    style = ('.srcs { display: grid; gap: 4px; margin-top: 0.5rem; }\n'
             '.src summary { cursor: pointer; padding: 4px 0; }\n'
             '.src pre.code { max-height: 70vh; overflow: auto; font-size: 12px; }\n</style>')
    page = page.replace('</style>', style, 1)
    # abrir el fichero al seguir un enlace a él
    js = ('<script>(function(){function o(){var e=location.hash&&document.getElementById(location.hash.slice(1));'
          'if(e&&e.tagName==="DETAILS"){e.open=true;e.scrollIntoView();}}addEventListener("hashchange",o);o();})();</script>\n</body>')
    page = page.replace('</body>', js, 1)
    for dst in (ROOT / 'report' / 'informe_axiones_completo.html', ROOT / 'results' / 'informe_axiones_completo.html'):
        dst.write_text(page)
        print(dst, f'{len(page) // 1024} KB, {len(files)} ficheros incrustados, {n_math}')


if __name__ == '__main__':
    main()
