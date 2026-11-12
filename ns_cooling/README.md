# ns_cooling

Neutron-star cooling: a 1D check of Marino et al. 2024 (Nat. Astron. 8, 1020; arXiv:2404.05371). See `LAB_NOTEBOOK.md`.

Third-party files are not committed. Before running, download:

- `data/src/bskfit18.f` from http://www.ioffe.ru/astro/NSG/BSk/bskfit18.f (Potekhin's fits for BSk22-26; version 13.02.23, sha256 `4bd3b716f04e40c69165fa83b5cd8ecb03c30aa1f2b45342fb070cef09a3a99c`). Only its coefficients are read.
- `data/src/sly4.dat` from http://www.ioffe.ru/astro/NSG/NSEOS/sly4.dat (Douchin & Haensel's SLy4 table).
- `paper/dh2001/dh2001.txt`: `pdftotext -layout` of arXiv:astro-ph/0111092 (Tables 4 and 6 are parsed from it).

Run a step from the repository root: `uv run python -m ns_cooling.step1_structure`.
