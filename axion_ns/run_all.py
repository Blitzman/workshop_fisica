"""Ejecuta todo el pipeline: acción -> fondo -> perturbaciones -> modelo simplificado -> números."""
import json
import platform
import time
from importlib.metadata import version
from pathlib import Path

from . import (step1_background, step2_perturbations, step3_toy_model, step4_numerics, step5_nonradial,
               step6_nonradial_numerics, step7_fullgr, step8_numerics,
               step9_axion_cowling, step10_axion_gr, step14_radial_branches,
               step15_radial_star, step16_radial_stability,
               step17_radial_validation, step18_rydberg,
               step19_avoided_crossing, step20_damped_count,
               step21_deep_artifacts, step22_exterior_contour,
               step23_eos_consistency, step24_scattering_width)

STEPS = [('Fondo (TOV + KG)', step1_background), ('Perturbaciones radiales', step2_perturbations),
         ('Modelo simplificado (Sec. IV)', step3_toy_model), ('Números y figura', step4_numerics),
         ('Paper II: no radiales en Cowling', step5_nonradial),
         ('Paper II: números y niveles ligados', step6_nonradial_numerics),
         ('Más allá de Cowling (nuevo)', step7_fullgr),
         ('Más allá de Cowling: números y figuras', step8_numerics),
         ('Con axión en Cowling: Tablas I, II y IV del Paper II', step9_axion_cowling),
         ('Con axión y métrica dinámica (nuevo)', step10_axion_gr),
         ('Modos radiales: umbral y ramas (modelo simplificado)', step14_radial_branches),
         ('Modos radiales de la estrella con condensado (GR)', step15_radial_star),
         ('Estabilidad radial: principio del argumento', step16_radial_stability),
         ('Modos radiales sin axión: validación externa (Font et al. 2002)', step17_radial_validation),
         ('Átomo gravitatorio: ley de Rydberg', step18_rydberg),
         ('Cruce evitado F - átomo gravitatorio', step19_avoided_crossing),
         ('Recuento en el semiplano amortiguado', step20_damped_count),
         ('Ceros profundos frente al radio de empalme', step21_deep_artifacts),
         ('Condición exterior frente al rayo complejo', step22_exterior_contour),
         ('Consistencia termodinámica de la EOS', step23_eos_consistency),
         ('Anchura desde la dispersión a frecuencia real', step24_scattering_width)]


def _cpu():
    try:
        return next(l.split(':', 1)[1].strip() for l in open('/proc/cpuinfo') if l.startswith('model name'))
    except (OSError, StopIteration):
        return platform.processor() or platform.machine()


def main():
    t0, out = time.time(), []
    for title, step in STEPS:
        t1 = time.time()
        checks = [{'name': n, 'ok': bool(ok), 'note': note} for n, ok, note in step.run()]
        out.append({'step': title, 'module': step.__name__.split('.')[-1], 'seconds': round(time.time() - t1, 1),
                    'checks': checks})
    allc = [c for s in out for c in s['checks']]
    bad = [c['name'] for c in allc if not c['ok']]
    secs = time.time() - t0
    print(f'\n{len(allc) - len(bad)}/{len(allc)} comprobaciones OK en {secs:.0f} s')
    for name in bad:
        print('  FAIL:', name)
    path = Path(__file__).resolve().parent.parent / 'results' / 'checks.json'
    env = {'python': platform.python_version(), 'platform': f'{platform.system()} {platform.machine()}',
           'cpu': _cpu(), 'date': time.strftime('%Y-%m-%d'),
           **{pkg: version(pkg) for pkg in ('sympy', 'mpmath', 'numpy', 'scipy', 'matplotlib')}}
    path.write_text(json.dumps({'seconds': round(secs), 'env': env, 'steps': out}, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
