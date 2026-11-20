# Agent Protocol for the archived axion package

Before writing code or running calculations in `archive/axion_ns/`, read `CONVENTIONS.md` and the current state in `LAB_NOTEBOOK.md`.

## 1. Verification Rules (Zero Hallucinations)
- **Everything from the action:** No intermediate equation is accepted without deriving it in SymPy from the action (`axion_ns/model.py`) and reducing its difference against the reference (`axion_ns/paper.py` or `axion_ns/nonradial.py`) to exact `0` modulo the background equations (`axion_ns/compare.py`).
- **Strict separation:** The file that derives equations (`model.py`, `fullgr.py`) must never import or know about the equations transcribed from the paper (`paper.py`).
- **Mandatory negative controls:** Every new step must include at least 2 deliberate mutations (sign flip, `4π → 8π`, `2l → 2(l+1)`) and verify that the comparator returns `False`.
- **Direct numerical bridge:** Never transcribe equations by hand into the numerical integrator; compile them from verified SymPy expressions (`lambdify` with `cse=True`, `ccode`, or `fcode`, as in `axion_ns/axgr.py`).

## 2. Token and Compute Economy
- **Never print raw SymPy expressions:** Never `print(expr)` expanded tensors or equations to the console. If an equation does not vanish, print only the list of variables with non-zero coefficients (`[v for v in PERT if expr.coeff(v) != 0]`) or evaluate the numerical residual at a fixed point.
- **Iterate in `--fast` mode (1 witness point):** During numerical development and debugging, run a single star and a single witness frequency (`ν_a = 6.0 kHz`, `l = 2`, < 3 s). Only launch full sweeps (`run_all.py`) once the witness point passes all tests.
- **Do not re-read full `.tex` or HTML files:** Extract equations from LaTeX once into `paper.py`. When editing reports, use line ranges (`StartLine`/`EndLine`).
- **No blind `simplify()` / `trigsimp()`:** On expressions with spherical harmonics, behavior depends on `PYTHONHASHSEED` and can hang for hours. Use exact polynomial reduction in `(cos θ, sin θ)` (`nonradial._divide`).

## 3. Lab Notebook (`LAB_NOTEBOOK.md`)
- Before testing a non-trivial idea (a limit, an approximation, or an integrator change), record the **Hypothesis** in `LAB_NOTEBOOK.md`.
- Once you have the result, record whether it was **Confirmed** or **Refuted** (and why) so dead ends are not repeated after a `/compact`.
- Any discrepancy with the paper or unexpected physical behavior is logged in the **Insights físicos y preguntas abiertas para los autores** section (`LAB_NOTEBOOK.md` §2) with a link to the file and line.

## 4. Rules Added After the Fabricated Report Incident (`LAB_NOTEBOOK.md` §3)
- **External validation before the axion:** no new result counts as validated unless it first reproduces, with the axion turned off, a published result by other authors. A `*_ref.py` written in the same session is not an independent reference.
- **Zero free constants in integrators:** introducing tuning coefficients (`CHI_*`, scale factors, "approximate" sensitivity formulas) is forbidden. Every quantity must come from a SymPy-derived equation or a published value with an exact citation, and the value in code must match the citation.
- **Derive, do not define:** an "identity" only counts if it comes from a subtraction yielding `0`. Writing `x = f(y)` and calling it an "identification" does not prove it.
- **Report numbers generated from data:** every number in the HTML or notebook is inserted from `results/*.json` in the build script. Do not type them by hand.
- **Physical premise reviewed by the group:** before writing code for a new line of work, record the premise and its order of magnitude in `LAB_NOTEBOOK.md` and validate it with José. Frame results as preliminary ("we obtain", "suggests"), not as discoveries.
