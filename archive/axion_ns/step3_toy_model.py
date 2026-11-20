"""Paso 3: modelo simplificado (Sec. IV del paper) derivado desde las Ec. (3.3) ya verificadas.

Aproximaciones del paper: sin perturbaciones métricas, fondo sin gradientes, simetría plana, y
    dm*/da ~ -beta sigma_N / (4 f_a),     V'' + n m*'' ~ m_a^2 (1 - n/n_c).
Además se comprueba en qué punto del potencial son válidas esas dos aproximaciones.
"""
import sympy as sp

from . import paper
from .compare import is_zero
from .model import r, w, V, mstar, eps0, mN, Phi0, lam0, a0, n0, dPhi, dlam, da1, dn1, xi1
from .pipeline import report

beta, sigma, fa, ma, eps_, mpi, fpi = sp.symbols('beta sigma_N f_a m_a epsilon m_pi f_pi', positive=True)
A, N, m2 = sp.symbols('A N m2')
nc = eps_ * mpi**2 * fpi**2 / sigma
ma2 = eps_ * beta * mpi**2 * fpi**2 / (4 * fa**2)  # m_a^2 del paper


def planar_limit(expr):
    """Sin métrica perturbada, fondo homogéneo (a0 = A, n0 = N), simetría plana (2/r -> 0)."""
    e = expr.subs({dPhi: 0, dlam: 0}).doit()
    e = e.subs({Phi0: 0, lam0: 0, a0: A, n0: N}).doit()
    e = e.xreplace({1 / r: 0})
    return sp.expand(e)


def run():
    res = []
    d = lambda f, k=1: sp.diff(f, r, k)
    # ---------------------------------------------------------- A) ecuaciones del toy desde (3.3)
    approx = {sp.diff(mstar(A), A): -beta * sigma / (4 * fa),
              sp.diff(mstar(A), A, 2): m2,
              sp.diff(V(A), A, 2): ma2 * (1 - N / nc) - N * m2}
    E8a, E8d, E8e = (planar_limit(paper.PERTURBATIONS[k]) for k in ('(3.3a) KG', '(3.3d) Euler', '(3.3e) continuidad'))
    xi_p = sp.solve(E8e, d(xi1))[0]
    toy1 = sp.expand(d(E8d).subs(d(xi1), xi_p)).subs(approx)
    toy2 = (E8a / fa).subs(approx)

    # variables del paper: delta p = (dp/dn) delta n,  delta a = f_a delta theta
    th = sp.Function('deltatheta')(r)
    dpdn = paper.dpdn.subs(n0, N).doit()
    epp = (paper.eps + paper.p).subs({n0: N, a0: A}).doit()
    cs2 = dpdn * N / epp
    dP = sp.Function('deltap')(r)
    to_paper = {dn1: dP / dpdn, da1: fa * th}
    toy1, toy2 = (sp.expand(e.subs(to_paper).doit()) for e in (toy1, toy2))

    paper_toy = lambda cs2, epp: (  # Ec. (4.2a,b) transcritas
        d(dP, 2) + w**2 / cs2 * dP - beta * eps_ * mpi**2 * fpi**2 / 4 * (N / nc) * d(th, 2),
        d(th, 2) + (w**2 - ma2 * (1 - N / nc)) * th + ma2 * (N / nc) * (dP / epp) / cs2)
    paper_toy1, paper_toy2 = paper_toy(cs2, epp)
    res.append(('Toy (Ec. 4.2a) se sigue de (3.3d)+(3.3e) con las aproximaciones del paper', is_zero(toy1 - paper_toy1), ''))
    res.append(('Toy (Ec. 4.2b) se sigue de (3.3a) con las aproximaciones del paper', is_zero(toy2 - paper_toy2), ''))

    # ---------------------------------------------------------- B) forma adimensional y matriz M
    al, la2, wb, k2, cs = sp.symbols('alpha lambda2 omegabar k2 c_s2', positive=True)
    EP, CS, pt, tt = sp.symbols('EP CS ptilde thetatilde')  # (eps+p), c_s^2 y amplitudes adimensionales
    P0 = beta * eps_ * mpi**2 * fpi**2 / 4                 # delta p = P0 * delta ptilde
    plane = {d(dP, 2): -k2 * ma**2 * P0 * pt, dP: P0 * pt, d(th, 2): -k2 * ma**2 * tt, th: tt}  # xbar = m_a x
    t1, t2 = paper_toy(CS, EP)
    t1, t2 = t1.subs(plane) / (ma**2 * P0), t2.subs(plane) / ma**2
    dimless = {w: wb * ma, N: al * nc, CS: cs, EP: ma**2 * fa**2 / (cs * la2)}  # lambda^2 = m_a^2 f_a^2/((eps+p) c_s^2)
    t1, t2 = (sp.simplify(e.subs(dimless).subs(fa, sp.sqrt(eps_ * beta) * fpi * mpi / (2 * ma))) for e in (t1, t2))
    Mmat = sp.Matrix([[wb**2 / cs + al**2 * la2, al * (wb**2 - 1 + al)], [al * la2, wb**2 - 1 + al]])
    pl = (Mmat - k2 * sp.eye(2)) * sp.Matrix([pt, tt])     # paper: v'' = -M v  ->  (M - k^2) v = 0
    # la fila 1 del paper es (fila 1 derivada) + alpha * (fila 2): mismo sistema
    ok_B = is_zero(t1 - (pl[0] - al * pl[1])) and is_zero(t2 - pl[1])
    res.append(('Forma adimensional -> matriz M de la Ec. (4.6)', ok_B, 'fila 1 = combinación de filas'))

    # ---------------------------------------------------------- C) relación de dispersión
    char = sp.factor((Mmat - k2 * sp.eye(2)).det())
    paper_det = (wb**2 - 1 + al - k2) * (wb**2 - k2 * cs) - al**2 * la2 * k2 * cs
    res.append(('det(M - k^2) == Ec. (4.7) (salvo factor c_s^2)', is_zero(char * cs - paper_det), ''))
    W = sp.Symbol('W')  # omegabar^2
    quad = sp.expand(paper_det.subs(wb**2, W))
    paper_quad = W**2 - W * (k2 * (1 + cs) + (1 - al)) + cs * k2 * (k2 + (1 - al) - al**2 * la2)
    res.append(('Cuadrática en omegabar^2 == Ec. (4.8)', is_zero(quad - paper_quad), ''))
    sols = sp.solve(paper_quad, W)
    disc = (k2 + 1 - al - k2 * cs)**2 + 4 * cs * k2 * al**2 * la2
    paper_sols = [(k2 + 1 - al + k2 * cs + s * sp.sqrt(disc)) / 2 for s in (1, -1)]
    res.append(('Raíces omegabar^2(k) == expresión del paper',
                all(any(is_zero(sp.expand(x - y)) for y in sols) for x in paper_sols), ''))
    ka2, ks2 = wb**2 - 1 + al, wb**2 / cs
    Delta = (ka2 - ks2 + al**2 * la2)**2 + 4 * al**2 * la2 * ks2
    kpm = [(ka2 + ks2 + al**2 * la2 + s * sp.sqrt(Delta)) / 2 for s in (1, -1)]
    res.append(('k_+-^2 y Delta (Ec. 4.12) son raíces de det(M - k^2) = 0',
                all(is_zero(sp.expand(char.subs(k2, kk))) for kk in kpm), ''))

    # ---------------------------------------------------------- D) ecuación trascendente (Ec. 4.13)
    kp, km, kq, kout, x0 = sp.symbols('k_p k_m k_a k_out x_0')
    # autovectores de M: v = (k^2 - k_a^2, alpha lambda^2)  (de la segunda fila)
    v = lambda k: sp.Matrix([k**2 - kq**2, al * la2])
    c, s_ = sp.cos, sp.sin
    # interior: sum_pm A_pm v_pm cos(k_pm x)  (regular en el centro);  exterior: dtheta = C e^{i kout x}, dp = 0
    # condiciones en x0: dp = 0, dtheta y dtheta' continuas  -> eliminando C:
    row_p = [v(k)[0] * c(k * x0) for k in (kp, km)]
    row_t = [v(k)[1] * (k * s_(k * x0) + sp.I * kout * c(k * x0)) for k in (kp, km)]
    det = sp.Matrix([row_p, row_t]).det() / (al * la2 * c(kp * x0) * c(km * x0))
    f_paper = (kq**2 - km**2) * kp * sp.tan(kp * x0) - (kq**2 - kp**2) * km * sp.tan(km * x0) \
        + sp.I * kout * (kp**2 - km**2)
    det = det.subs({s_(kp * x0): sp.tan(kp * x0) * c(kp * x0), s_(km * x0): sp.tan(km * x0) * c(km * x0)})
    ok_D = sp.simplify(sp.expand(det - f_paper)) == 0
    res.append(('Ajuste interior/exterior con dp(x0)=0 reproduce f(omega) de la Ec. (4.13)', ok_D,
                'sqrt(Delta) = k_+^2 - k_-^2'))
    # límite lambda -> 0: {k_+, k_-} -> {k_a, k_s}; f -> (k_s^2 - k_a^2)(k_a tan(k_a x0) + i k_out)
    ks = sp.Symbol('k_s')
    f0 = sp.simplify(f_paper.subs({kp: ks, km: kq}))
    res.append(('Límite lambda->0: f ∝ k_a tan(k_a x0) + i k_out  (Ec. 4.11)',
                sp.simplify(f0 / (ks**2 - kq**2) - (kq * sp.tan(kq * x0) + sp.I * kout)) == 0, ''))

    # ---------------------------------------------------------- E) ¿dónde valen las aproximaciones?
    th_ = sp.Symbol('theta', real=True)
    g = sp.sqrt(1 - beta * sp.sin(th_ / 2)**2)
    m_star = mN - sigma * (1 - g)               # Ec. (2.4)
    Vfun = eps_ * mpi**2 * fpi**2 * (1 - g)     # Ec. (2.1a)
    dm = sp.simplify(sp.diff(m_star, th_) / fa)                       # dm*/da
    U2 = sp.simplify((sp.diff(Vfun, th_, 2) + al * nc * sp.diff(m_star, th_, 2)) / fa**2)
    U2_0 = sp.simplify(U2.subs(th_, 0) / ma2)
    U2_pi = sp.simplify(U2.subs(th_, sp.pi) / ma2)
    res.append(('V\'\' + n m*\'\' en theta=0 vale m_a^2 (1 - alpha)  [la aproximación del paper]',
                sp.simplify(U2_0 - (1 - al)) == 0, ''))
    res.append(('dm*/da = -beta sigma sin(theta)/(4 f_a g): se ANULA en theta = 0 y theta = pi',
                dm.subs(th_, 0) == 0 and sp.simplify(dm.subs(th_, sp.pi)) == 0,
                f'dm*/da = {sp.simplify(dm)}'))
    res.append(('En theta=pi (fondo del toy si alpha>1): m_eff^2/m_a^2 = (alpha-1)/sqrt(1-beta) > 0',
                sp.simplify(U2_pi - (al - 1) / sp.sqrt(1 - beta)) == 0,
                'el toy usa 1-alpha < 0 (taquiónico), que corresponde a theta=0'))
    return report('Paso 3: modelo simplificado (Sec. IV)', res)


if __name__ == '__main__':
    run()
