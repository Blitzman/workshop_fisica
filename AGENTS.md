# Protocolo agéntico para astrofísica relativista

Antes de escribir código o lanzar cálculos, lee `CONVENTIONS.md` y el estado actual en `LAB_NOTEBOOK.md`.

## 1. Reglas de verificación (cero alucinaciones)
- **Todo desde la acción:** Ninguna ecuación intermedia se acepta sin deducirla en SymPy desde la acción (`axion_ns/model.py`) y reducir su diferencia con la referencia (`axion_ns/paper.py` o `axion_ns/nonradial.py`) a `0` exacto módulo las ecuaciones de fondo (`axion_ns/compare.py`).
- **Separación estricta:** El fichero que deriva (`model.py`, `fullgr.py`) jamás importa ni conoce las ecuaciones transcritas del paper (`paper.py`).
- **Controles negativos obligatorios:** Todo paso nuevo debe incluir al menos 2 mutaciones deliberadas (cambio de signo, `4π → 8π`, `2l → 2(l+1)`) y verificar que el comparador devuelve `False`.
- **Puente numérico directo:** Prohibido transcribir ecuaciones a mano al integrador numérico; compílalas desde las expresiones de SymPy verificadas (`lambdify` con `cse=True`, `ccode` o `fcode`, como en `axion_ns/axgr.py`).

## 2. Economía de tokens y tiempo de cómputo
- **Prohibido imprimir expresiones de SymPy crudas:** Nunca hagas `print(expr)` de tensores o ecuaciones desplegadas en consola. Si una ecuación no anula, imprime solo la lista de variables cuyo coeficiente es no nulo (`[v for v in PERT if expr.coeff(v) != 0]`) o evalúa el residuo numérico en un punto fijo.
- **Iterar en modo `--fast` (1 punto testigo):** Durante el desarrollo y depuración numérica, ejecuta una sola estrella y una sola frecuencia testigo (`ν_a = 6.0 kHz`, `l = 2`, < 3 s). Solo lanza barridos completos (`run_all.py`) cuando el punto testigo ya pasa todos los tests.
- **No releer `.tex` ni HTML completos:** Extrae las ecuaciones del LaTeX una sola vez a `paper.py`. Para editar informes, usa rangos de líneas (`StartLine`/`EndLine`).
- **Prohibido `simplify()` / `trigsimp()` ciego:** En expresiones con armónicos esféricos depende de `PYTHONHASHSEED` y puede colgarse horas. Usa reducción polinómica exacta en `(cos θ, sin θ)` (`nonradial._divide`).

## 3. Cuaderno de laboratorio (`LAB_NOTEBOOK.md`)
- Antes de probar una idea no trivial (un límite, una aproximación o un cambio de integrador), anota la **Hipótesis** en `LAB_NOTEBOOK.md`.
- Cuando obtengas el resultado, anota si fue **Confirmada** o **Refutada** (y por qué) para no repetir caminos muertos tras un `/compact`.
- Cualquier discrepancia con el paper o comportamiento físico inesperado se registra en la sección **Insights y preguntas para el grupo** con enlace al fichero y línea.
