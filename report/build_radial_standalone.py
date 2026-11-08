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
    if key.endswith('uv.lock.txt'):
        continue                                   # 130 KB de dependencias fijadas: no es código del cálculo
    name = key.split('/', 1)[1].replace('__', '/').removesuffix('.txt')
    text = (ROOT / files[key]).read_text()
    blocks.append(f'<details class="src"><summary><code>{html.escape(name)}</code> '
                  f'<span class="small">({len(text.splitlines())} líneas)</span></summary>'
                  f'<pre class="code"><code>{html.escape(text)}</code></pre></details>')
embedded = '<div class="srcs">' + ''.join(blocks) + '</div>'
body = re.sub(r'<ul class="files">.*?</ul>', lambda m: embedded, body, count=1, flags=re.S)
body = body.replace('El código y los datos están publicados junto a esta página (los <code>.py</code> van como <code>.py.txt</code>).',
                    'El código y los datos van incrustados abajo: despliega cada fichero para verlo.')
head = head.replace('</style>', '.srcs { display: grid; gap: 6px; }\n'
                    '.src summary { cursor: pointer; padding: 6px 0; }\n'
                    '.src pre.code { max-height: 70vh; overflow: auto; font-size: 12px; }\n</style>')
out = ('<!doctype html>\n<html lang="es">\n<head>\n<meta charset="utf-8">\n'
       '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
       f'{head}\n</head>\n<body>\n{body}\n</body>\n</html>\n')
dst = ROOT / 'report' / 'informe_radial_completo.html'
dst.write_text(out)
print(dst, f'{len(out) // 1024} KB, {len(blocks)} ficheros incrustados')
