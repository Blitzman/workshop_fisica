"""Genera report/informe_radial_completo.html: el informe radial en un solo fichero, con el código y los datos
incrustados (sin ficheros auxiliares). Se ejecuta después de report/build_radial.py."""
import html
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
page = (ROOT / 'report' / 'informe_radial.html').read_text()
files = json.loads((ROOT / 'report' / 'radial_files.json').read_text())

head, body = page.split('<div class="wrap">', 1)
body = '<div class="wrap">' + body
blocks = []
for key in sorted(files):
    name = key.split('/', 1)[1].removesuffix('.txt').replace('__', '/').replace('//init/.py', '/__init__.py')
    text = (ROOT / files[key]).read_text()
    blocks.append(f'<details class="src"><summary><code>{html.escape(name)}</code> '
                  f'<span class="small">({len(text.splitlines())} líneas)</span></summary>'
                  f'<pre class="code" data-path="{html.escape(name)}"><code>{html.escape(text)}</code></pre></details>')
EXTRACT = """import html, re, sys
from pathlib import Path
src, dst = Path(sys.argv[1]).read_text(), Path(sys.argv[2])
for path, body in re.findall(r'<pre class="code" data-path="([^"]+)"><code>(.*?)</code></pre>', src, flags=re.S):
    f = dst / html.unescape(path)
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(html.unescape(body))
    print(f)"""
extract_html = ('<p><b>Extraer todos los ficheros</b> a una carpeta (con Python 3, sin dependencias): guarda esto como '
                '<code>extraer.py</code> y ejecuta <code>python3 extraer.py informe_radial_completo.html carpeta</code>. '
                'Después, dentro de <code>carpeta</code>, funcionan los comandos de arriba.</p>'
                f'<pre class="code"><code>{html.escape(EXTRACT)}</code></pre>')
embedded = extract_html + '<div class="srcs">' + ''.join(blocks) + '</div>'
body = re.sub(r'<ul class="files">.*?</ul>', lambda m: embedded, body, count=1, flags=re.S)
body = body.replace('El código y los datos están publicados junto a esta página (los <code>.py</code> van como <code>.py.txt</code>).',
                    'El código y los datos van incrustados abajo: despliega cada fichero para verlo.')
head = head.replace('</style>', '.srcs { display: grid; gap: 6px; }\n'
                    '.src summary { cursor: pointer; padding: 6px 0; }\n'
                    '.src pre.code { max-height: 70vh; overflow: auto; font-size: 12px; }\n</style>')
out = ('<!doctype html>\n<html lang="es">\n<head>\n<meta charset="utf-8">\n'
       '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
       f'{head}\n</head>\n<body>\n{body}\n</body>\n</html>\n')
for dst in (ROOT / 'report' / 'informe_radial_completo.html', ROOT / 'results' / 'informe_radial_completo.html'):
    dst.write_text(out)   # las dos copias, como build_radial.py con informe_radial.html
    print(dst, f'{len(out) // 1024} KB, {len(blocks)} ficheros incrustados')
