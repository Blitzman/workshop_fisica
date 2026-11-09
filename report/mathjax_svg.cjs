// Convierte las fórmulas TeX de una página HTML en SVG incrustado (sin MathJax en el navegador).
// Uso: node report/mathjax_svg.cjs < fragmento.html > salida.html
// Entra un fragmento de <body>; sale la hoja de estilo de MathJax (<style>), una línea <!--MJX-->, y el fragmento.
// Necesita mathjax-full 3.2.2 (fijado en report/package-lock.json): npm ci --prefix report
const path = require('path');
const MJ = process.env.MATHJAX_DIR || path.join(__dirname, 'node_modules', 'mathjax-full');
const { mathjax } = require(path.join(MJ, 'js/mathjax.js'));
const { TeX } = require(path.join(MJ, 'js/input/tex.js'));
const { AllPackages } = require(path.join(MJ, 'js/input/tex/AllPackages.js'));
const { SVG } = require(path.join(MJ, 'js/output/svg.js'));
const { liteAdaptor } = require(path.join(MJ, 'js/adaptors/liteAdaptor.js'));
const { RegisterHTMLHandler } = require(path.join(MJ, 'js/handlers/html.js'));

const html = require('fs').readFileSync(0, 'utf8');
const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);
const doc = mathjax.document('<!DOCTYPE html><html><head></head><body>' + html + '</body></html>', {
  InputJax: new TeX({ packages: AllPackages, inlineMath: [['\\(', '\\)']], displayMath: [['\\[', '\\]']] }),
  OutputJax: new SVG({ fontCache: 'global' }),
});
doc.render();
const n = Array.from(doc.math).length;
const errors = Array.from(doc.math).filter((m) => adaptor.outerHTML(m.typesetRoot).includes('data-mjx-error')).length;
process.stdout.write(adaptor.outerHTML(doc.outputJax.styleSheet(doc)) + '\n<!--MJX-->\n' + adaptor.innerHTML(adaptor.body(doc.document)));
process.stderr.write(`${n} fórmulas, ${errors} con error\n`);
if (errors) process.exit(2);
