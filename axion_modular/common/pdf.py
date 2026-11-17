"""Typeset an equation artifact to PDF (human-readable companion of the .eq.json).

The PDF is OUTPUT ONLY: numerical code must load the .eq.json (compiled with lambdify), never
retype equations from the PDF. The header carries the content hash of the JSON so the two
cannot drift apart unnoticed.

Long expressions are broken into lines of a few terms. `max_terms` limits how many terms of one
equation are typeset (the rest is announced explicitly: "N more terms in the JSON").
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import sympy as sp

_TEX_SPECIAL = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#', '_': r'\_', '{': r'\{',
                '}': r'\}', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}'}


def tex_escape(s: str) -> str:
    return ''.join(_TEX_SPECIAL.get(c, c) for c in str(s))


def _terms(expr):
    return list(sp.Add.make_args(expr))


def latex_lines(expr, per_line=3, max_terms=None):
    """Return (lines, omitted, expanded): LaTeX strings, one per line of the typeset equation.

    A single long term (e.g. 1/r^2 * (long sum)) would be shrunk to illegibility; for DISPLAY ONLY it is
    expanded into a sum so it can be broken into lines (`expanded` = True; the artifact is unchanged).
    """
    terms = _terms(expr)
    expanded = False
    if len(terms) == 1 and len(sp.latex(expr)) > 300:
        ex = sp.expand(expr)
        if len(_terms(ex)) > 1:
            terms, expanded = _terms(ex), True
    omitted = 0
    if max_terms is not None and len(terms) > max_terms:
        omitted = len(terms) - max_terms
        terms = terms[:max_terms]
    # greedy filling of lines by a character budget (per_line is kept for API compatibility: budget = 110 * per_line)
    budget = 110 * per_line
    lines, cur, cur_len = [], [], 0
    for t_ in terms or [sp.Integer(0)]:
        tl = sp.latex(t_)
        if cur and cur_len + len(tl) > budget:
            lines.append(cur)
            cur, cur_len = [], 0
        cur.append(tl)
        cur_len += len(tl) + 3
    lines.append(cur)
    out_lines = []
    for i, ch in enumerate(lines):
        s_ = ' + '.join(ch).replace('+ -', '- ')
        if i > 0 and not s_.startswith('-'):
            s_ = '+ ' + s_
        out_lines.append(s_)
    lines = out_lines
    return lines, omitted, expanded


_HEAD = r'''\documentclass[10pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,adjustbox,xcolor,fancyhdr,hyperref}
\usepackage[margin=16mm]{geometry}
\setlength{\parindent}{0pt}\setlength{\parskip}{3pt}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small %(short)s}\fancyhead[R]{\small sha256 %(sha)s}\fancyfoot[C]{\small\thepage}
\hypersetup{colorlinks=true,urlcolor=blue}
\begin{document}
{\Large\bfseries %(title)s}\par
\smallskip
%(meta)s
\par\medskip\hrule\medskip
'''


def _build_tex(title, header, entries, per_line, max_terms, sha):
    meta = r'\begin{tabular}{@{}lp{0.62\textwidth}@{}}' + ''.join(
        f'\\textbf{{{tex_escape(k)}}} & {tex_escape(v)}\\\\' for k, v in header.items()) + r'\end{tabular}'
    out = [_HEAD % {'title': tex_escape(title), 'meta': meta, 'sha': tex_escape(sha[:16]),
                    'short': tex_escape(header.get('artifact', title))}]
    out.append(r"\textit{Generated from the JSON artifact; do not retype equations from this file. "
               r"Each expression is the left side of an equation ``$= 0$'' unless stated otherwise.}\par" + "\n")
    for e in entries:
        lines, omitted, expanded = latex_lines(e['expr'], per_line, max_terms)
        sub = f"\\textbf{{{tex_escape(e['key'])}}}"
        if e.get('label'):
            sub += f"\\quad {tex_escape(e['label'])}"
        out.append(r'\bigskip\par\noindent' + sub + r'\par')
        if expanded:
            out.append(r"{\small\textit{display:} expanded for typesetting (same expression).}\par")
        for field in ('ref', 'note'):
            if e.get(field):
                out.append(f"{{\\small\\textit{{{field}:}} {tex_escape(e[field])}}}\\par")
        for j, ln in enumerate(lines):
            tail = r'\;=\;0' if (j == len(lines) - 1 and not omitted) else ''
            out.append(r'\noindent\adjustbox{max width=\textwidth}{$\displaystyle ' + ln + tail + r'$}\par')
        if omitted:
            out.append(r'{\small\color{red}\ldots\ %d more terms are in the JSON artifact (not typeset here).}\par' % omitted)
    out.append(r'\end{document}')
    return '\n'.join(out)


def write_equations_pdf(pdf_path, title, header: dict, entries: list, content_sha256: str, per_line=3, max_terms=None,
                        timeout=900):
    """entries: list of dict(key, expr, label, ref, note). Returns dict(pdf, tex, max_terms_used, ok)."""
    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    attempts = [max_terms]
    if max_terms is None:
        attempts += [300, 80, 20]
    else:
        attempts += [max(10, max_terms // 3), 10]
    last_err = ''
    for mt in attempts:
        tex = _build_tex(title, header, entries, per_line, mt, content_sha256)
        with tempfile.TemporaryDirectory() as td:
            tf = Path(td) / 'doc.tex'
            tf.write_text(tex, encoding='utf-8')
            try:
                proc = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', 'doc.tex'], cwd=td,
                                      capture_output=True, text=True, timeout=timeout)
            except subprocess.TimeoutExpired:
                last_err = f'timeout with max_terms={mt}'
                continue
            if proc.returncode == 0 and (Path(td) / 'doc.pdf').exists():
                shutil.copy(Path(td) / 'doc.pdf', pdf_path)
                pdf_path.with_suffix('.tex').write_text(tex, encoding='utf-8')
                return {'pdf': str(pdf_path), 'tex': str(pdf_path.with_suffix('.tex')), 'max_terms_used': mt, 'ok': True}
            last_err = (proc.stdout or '')[-600:]
    return {'pdf': None, 'tex': None, 'max_terms_used': None, 'ok': False, 'error': last_err}


def pdf_from_artifact(pdf_path, meta: dict, exprs: dict, per_line=3, max_terms=None):
    """Convenience: build the PDF of an artifact given its meta dict (from save_equations) and the expressions."""
    entries = [{'key': k, 'expr': exprs[k], 'label': v.get('label', ''), 'ref': v.get('ref', ''), 'note': v.get('note', '')}
               for k, v in meta['equations'].items()]
    header = {'artifact': f"{meta['name']}.eq.json", 'kind': meta['kind'], 'producer': meta['producer']['module'],
              'git': str(meta['producer']['git'].get('commit')), 'created': meta['producer']['created'],
              'content sha256': meta['content_sha256'][:32] + '...'}
    for i in meta.get('inputs', []):
        header[f"input {i['name']}"] = i['sha256'][:16] + '...'
    if meta.get('description'):
        header['description'] = meta['description']
    return write_equations_pdf(pdf_path, meta.get('title') or meta['name'], header, entries, meta['content_sha256'],
                               per_line, max_terms)
