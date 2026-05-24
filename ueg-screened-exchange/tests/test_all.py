"""Fast self-tests.  Run: python -m pytest tests -q   (or: python tests/test_all.py)"""
import os
import sys
import math
import numpy as np
from scipy.integrate import quad, dblquad

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from ueg_exchange import lowdensity as LD
from ueg_exchange import constants as C, lindhard as L, lorentzian as B, dispersion as T
from ueg_exchange import highdensity as H, ringsum as R


def test_A0_bookkeeping():
    val = 2 * 2 * (4 * np.pi) ** 2 * (2 * np.pi) ** -9 * (3 * np.pi ** 2) * 0.25
    assert abs(val - C.A0) < 1e-18


def test_onsager():
    assert abs(C.E2X - 0.0483583178362888) < 1e-14


def test_lindhard_limits():
    assert abs(L.lindhard_shape(1e-6, 0.0) - 1.0) < 1e-6           # lindhard_shape(0,0) = 1
    assert abs(L.brm_lindhard_g(1e-6, 0.0) - 2.0) < 1e-6       # BRM normalization
    # large-q expansion at u=0:  lindhard_shape = 4/(3q^2) + 16/(15 q^4) + O(q^-6)
    for q in (10.0, 20.0):
        ref = 4 / (3 * q ** 2) + 16 / (15 * q ** 4)
        assert abs(L.lindhard_shape(q, 0.0) - ref) < 2e-3 * ref


def test_coupling_screening_factor():
    for c in [1e-9, 1e-3, 1.0, 50.0]:
        ref = 2 * quad(lambda lam: lam / (1 + lam * c), 0, 1)[0]
        assert abs(L.coupling_screening_factor(np.array([c]))[0] - ref) < 1e-10


def test_kappa_closed_form():
    for Lam in [0.05, 0.7, 5.0]:
        assert abs(B.kappa(Lam) - B.kappa_quadrature(Lam)) < 1e-12


def test_A_shell_zero_sum():
    """int_0^inf A(z) dz = 0: the bare same-sign frequency kernel cancels."""
    assert abs(quad(B.A_shell, 0, np.inf, limit=400)[0]) < 1e-9


def test_block_large_q_kernel():
    """p_loc(s_lambda) = 2 pi s_lambda/(1+2 s_lambda)^2 against its defining frequency integral."""
    for s_lambda in [0.3, 0.5, 1.0, 3.0]:
        ref = quad(lambda u: s_lambda ** 2 * (0.25 - u ** 2) /
                   ((u ** 2 + s_lambda ** 2) * (u ** 2 + 0.25) ** 2), 0, np.inf, limit=300)[0]
        assert abs(B.p_loc(s_lambda) - ref) < 1e-10


def test_dispersion_reproduces_A1():
    assert abs(T.A1_dispersion() - C.A1) < 1e-12


def test_spin_identity():
    """e_kite(rs,1) = e_kite(2^-4/3 rs, 0) for the three-term expansion."""
    for rs in [2.0, 5.0, 1e4]:
        a = C.asymptotic_kite(rs, 1.0)
        b = C.asymptotic_kite(rs * 2 ** (-4 / 3), 0.0)
        assert abs(a - b) < 1e-14 * max(1.0, abs(a))


def test_spin_moments():
    assert abs(C.spin_moment(4, 1.0) - 2 ** (1 / 3)) < 1e-14
    assert abs(C.spin_moment(6, 1.0) - 2.0) < 1e-14
    assert abs(C.spin_moment(3, 0.4) - 1.0) < 1e-14  # fixed total density for every zeta
    x = [(1 + 0.4) ** (1 / 3), (1 - 0.4) ** (1 / 3)]
    for n in (4, 5, 6):
        assert abs(C.spin_moment(n, 0.4) - 0.5 * sum(xx ** n for xx in x)) < 1e-14


def test_calA_and_Cln2():
    num, closed = H.calA()
    assert abs(num - closed) < 1e-9
    assert abs(H.C_ln2_from_calA() - C.C_LN2) < 1e-15
    assert abs(H.int_h2() - np.pi / 3 * (1 - np.log(2))) < 1e-12


def test_ward_cancellation():
    """Pi_vtx^(1)^log + Pi_SE^(1)^log = 0 for every z (Prop. ward)."""
    for z in (0.0, 0.3, 1.0, 5.0):
        a, b = H.response_log_coefficients(z, 1.0)
        assert abs(a + b) < 1e-15


def test_spin_identity_ln2():
    """C_ln2(1) = 2^{-4/3} C_ln2(0); Phi_s(1) = 1 - log 2; I(theta) closed form vs quadrature."""
    assert abs(H.rho_zeta(1.0) - 2 ** (-4 / 3)) < 1e-12
    assert abs(H.rho_zeta(0.0) - 1.0) < 1e-12
    assert abs(H.Phi_spin(1.0) - (1 - np.log(2))) < 1e-14
    h = lambda z: 1 - z * np.arctan(1 / z)
    for th in (0.5, 2.0):
        num = quad(lambda z: h(th * z) * H.alpha_shell(z), 0, np.inf, limit=400)[0]
        assert abs(num - H.I_theta(th)) < 1e-9


def test_ring_against_BRM_tableVIII():
    assert abs(1000 * R.eps_r(1.0) + 157.6) < 0.15
    assert abs(1000 * R.eps_r(5.0) + 85.0) < 0.15


def test_ring_c0():
    """eps_r - c_L log rs -> c0 = -0.1423 Ry (BRM Table II)."""
    rs = 1e-5
    c0 = R.eps_r(rs) - C.C_L_RING * np.log(rs)
    assert abs(c0 + 0.1423) < 2e-4


def test_lindhard_static_and_small_nu():
    """Closed form vs f_L at nu=0 (q<2 included) and vs adaptive quadrature at small nu (B1-01)."""
    for q in (0.3, 0.5, 1.0, 1.048, 0.983723255150397, 1.999, 3.0):
        assert abs(L.lindhard_shape(q, 0.0) - L.f_L(q / 2)) < 1e-10, q
        for nu in (1e-6, 1e-3, 0.1, 1.0):
            a = 0.5 * q * q
            f = lambda k: k * np.log(((q * k + a) ** 2 + nu ** 2) / ((q * k - a) ** 2 + nu ** 2))
            pts = [0.0, 1.0]
            if q < 2:
                k0 = q / 2
                pts += [k0] + [k0 + s * d * nu / q for s in (-1, 1) for d in (1, 10, 100, 1000) if 0 < k0 + s * d * nu / q < 1]
            pts = sorted(set(pts))
            ref = sum(quad(f, pts[i], pts[i + 1], limit=400)[0] for i in range(len(pts) - 1)) / (2 * q)
            assert abs(L.lindhard_shape(q, nu) - ref) < 1e-7 * max(1.0, abs(ref)), (q, nu, L.lindhard_shape(q, nu), ref)
    assert np.isfinite(L.lindhard_shape(0.983723255150397, 0.0))
    # scalar input on every branch (M1-7), and the exact q -> 0 limit
    assert isinstance(L.lindhard_shape(0.01, 10.0), float)          # far branch, scalar
    assert isinstance(L.lindhard_shape(0.5, 0.0), float)            # main branch, scalar
    assert abs(L.lindhard_shape(0.0, 0.0) - 1.0) < 1e-14
    assert abs(L.lindhard_shape(1e-10, 0.0) - 1.0) < 1e-12
    h = lambda z: 1 - z * np.arctan(1 / z)
    for q, nu in ((1e-5, 1e-5), (1e-6, 1e-2)):
        assert abs(L.lindhard_shape(q, nu) - h(nu / q)) < 1e-9 * max(1.0, h(nu / q))


def test_A3_from_profiles():
    """A3 from the profiles: closed 1-D reductions (App. A) and raw profile quadratures."""
    cl = LD.A3_pieces_claimed()
    # exact Beta-function reductions
    assert abs(LD.A3_cross_closed() - cl["cross"]) < 1e-12
    assert abs(LD.A3_recoil_closed() - cl["recoil"]) < 1e-12
    assert abs(LD.A3_lindhard_closed() - cl["lindhard"]) < 1e-8
    assert abs(LD.A3_from_profiles() - C.A3) < 1e-8
    # raw multi-dimensional profiles agree with them
    assert abs(LD.A3_cross_numeric() - cl["cross"]) < 1e-8
    assert abs(LD.A3_recoil_numeric() - cl["recoil"]) < 2e-4 * abs(cl["recoil"])
    assert abs(LD.A3_lindhard_numeric() - cl["lindhard"]) < 2e-4 * abs(cl["lindhard"])


def test_A3_frequency_identities():
    """The two contour identities behind the recoil and Lindhard reductions (App. A)."""
    for m in (0.3, 1.0, 7.0):
        rec = quad(lambda u: (1 + 4 * u * u) * ((1j * u - 0.5) ** -4).real / (1 + 4 * u * u + m),
                   0, np.inf, limit=400)[0] / np.pi
        assert abs(rec - LD.I_recoil(m)) < 1e-8 * max(1.0, abs(rec)), m
        lin = -16.0 / 15 * quad(lambda u: (12 * u * u - 1) * ((1j * u - 0.5) ** -2).real
                                / ((1 + 4 * u * u) * (1 + 4 * u * u + m) ** 2), 0, np.inf, limit=400)[0] / np.pi
        assert abs(lin - LD.U_lindhard(m)) < 1e-8 * max(1.0, abs(lin)), m


def test_large_z_series():
    """h(z) and alpha(z) must not cancel to zero at large z (C0-1, C0-2)."""
    for z in (1e3, 1e4, 1e8, 1e12):
        assert abs(L.h_shell(z) * 3 * z * z - 1.0) < 1e-6, z
    assert abs(L.lindhard_shape(1e-8, 1.0) * 3 * (1e8) ** 2 - 1.0) < 1e-6
    # alpha(z) = -1/(3z^2) + 3/(5z^4) - 5/(7z^6) + ...
    for z in (1e3, 1e4, 1e8):
        a = B.A_shell(z) / (8 * np.pi ** 2)
        ser = -1 / (3 * z ** 2) + 3 / (5 * z ** 4) - 5 / (7 * z ** 6)
        assert abs(a - ser) < 1e-6 * abs(ser), z            # no cancellation, matches the series
        assert a < 0, z
    assert abs(H.alpha_shell(1e8) + 1 / (3 * 1e16)) < 1e-22


def test_block_endpoints():
    assert B.kappa(0.0) == 0.0 and B.kappa(np.inf) == 0.0 and B.p_loc(np.inf) == 0.0
    assert abs(B.kappa(1e4) - B.kappa_quadrature(1e4)) < 1e-14
    assert abs(B.A_shell(1e8) + 8 * np.pi ** 2 / (3e16)) < 1e-18
    assert L.coupling_screening_factor(np.array([0.0, np.inf]))[1] == 0.0



def test_ward_from_the_insertion_integral():
    """Prop. ward, checked from the defining integral rather than assumed.

    self_energy_response_log_slope evaluates the exchange self-energy insertion directly and
    fits its log(1/q) coefficient; it must equal +q_tf_sq^2 alpha(z)/4, the
    negative of the vertex coefficient, for every z.
    """
    for z in (0.0, 0.3, 1.0, 3.0):
        h = 1 - z * np.arctan(1 / z) if z > 0 else 1.0
        hp = (-np.arctan(1 / z) + z / (1 + z * z)) if z > 0 else -np.pi / 2
        got = H.self_energy_response_log_slope(z)
        ref = (h + z * hp) / 4.0
        assert abs(got - ref) < 2e-6 * max(abs(ref), 1e-3), (z, got, ref)


def test_mellin_closed_forms():
    """Xi, K(s) on -1<Re s<0, and b0 = gamma/2 - 3/8."""
    from scipy.integrate import quad
    for z in (0.3, 1.0, 5.0):
        ref = quad(lambda t: z * np.exp(-z * t) / (1 + t * t), 0, np.inf, limit=400)[0]
        assert abs(B.Xi_closed(z) - ref) < 1e-10
    for s in (-0.9, -0.5, -0.2):
        A = 60.0
        core = quad(lambda z: z ** (s - 1) * B.Xi_closed(z), 0, A, limit=800)[0]
        tail = quad(lambda z: z ** (s - 1) * (B.Xi_closed(z) - 1.0), A, np.inf,
                    limit=800)[0] - A ** s / s
        assert abs(B.K_mellin(s) - (core + tail)) < 1e-3 * abs(B.K_mellin(s))
    assert abs(B.b0_closed() - (0.5 * np.euler_gamma - 0.375)) < 1e-15

def test_helpers_at_extreme_arguments():
    """The three helpers that carry a cancellation must survive large arguments.

    f_L and Phi_spin are both written as a difference that loses every digit of
    the answer at large argument; yukawa_pole has to bracket a root whose
    position grows with r_s. Each is compared against a form that does not
    cancel, computed here.
    """
    from ueg_exchange.lindhard import f_L
    from ueg_exchange.highdensity import Phi_spin, yukawa_pole

    # f_L(x) = sum_{k>=1} x^{-2k}/((2k-1)(2k+1)) at large x, and must stay positive
    for x in (30.0, 1e3, 1e6, 1e8):
        ref = sum(x ** (-2 * k) / ((2 * k - 1) * (2 * k + 1)) for k in range(1, 12))
        got = float(f_L(np.array([x]))[0])
        assert got > 0.0, f"f_L({x}) = {got} is not positive"
        assert abs(got - ref) <= 1e-12 * ref, (x, got, ref)

    # Phi_spin(a) -> log a - 2/3 + 3/(2a) - 9/(10 a^2) + ...
    for a in (1e3, 1e6, 1e10):
        ref = (math.log(a) - 2.0 / 3.0 + 1.5 / a - 0.9 / a ** 2 + (2.0 / 3.0) / a ** 3)
        got = float(Phi_spin(np.array([a]))[0])
        assert abs(got - ref) <= 1e-12 * abs(ref), (a, got, ref)
    # and still reproduces the printed spin ratios
    from ueg_exchange.highdensity import rho_zeta
    for zeta, want in ((0.0, 1.0), (0.4, 0.956279), (0.8, 0.781213), (1.0, 2 ** (-4 / 3))):
        assert abs(rho_zeta(zeta) - want) < 2e-6, (zeta, rho_zeta(zeta), want)

    # the screened pole must be found at every density and approach pi q_tf_sq / 8
    aL = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
    for rs in (1.0, 200.0, 1e5, 1e8):
        kappa, Z = yukawa_pole(rs)
        q_tf_sq = 4.0 * aL * rs / np.pi
        assert kappa > 0.0 and 1.0 < Z <= 2.0 + 1e-12, (rs, kappa, Z)
        if rs >= 1e5:
            assert abs(kappa / (np.pi * q_tf_sq / 8.0) - 1.0) < 1e-5, (rs, kappa)
            assert abs(Z - 2.0) < 1e-5, (rs, Z)


def test_auxiliary_functions_at_extreme_arguments():
    """The three helpers outside the production path, where cancellation bites.

    f_L, Phi_spin and yukawa_pole are each written so that the printed closed
    form is replaced by a series where it would subtract two nearly equal
    quantities. These are the arguments at which the unstabilized forms fail.
    """
    from ueg_exchange.lindhard import f_L
    from ueg_exchange.highdensity import Phi_spin, yukawa_pole

    # f_L(x) -> 1/(3x^2) + 1/(15x^4) + ...; the closed form subtracts 1/2 from
    # a quantity tending to 1/2 and returns garbage, and then a negative value.
    for x in (1e4, 1e6, 1e8, 1e12):
        got = float(f_L(x))
        ref = 1.0 / (3 * x * x) + 1.0 / (15 * x ** 4)
        assert got > 0.0, (x, got)
        assert abs(got / ref - 1.0) < 1e-9, (x, got, ref)
    # the two branches agree at the switch: compare the series used above 30
    # with the closed form evaluated at the same argument just below it
    x = 30.0
    closed = 0.5 + (1 - x * x) / (4 * x) * np.log(abs((1 + x) / (1 - x)))
    assert abs(float(f_L(x)) / closed - 1.0) < 1e-10, (float(f_L(x)), closed)

    # Phi_spin(a) -> log a - 2/3 + 3/(2a) - ...; the direct form returned zero
    # for large a because 2a^2 - a cancels against the expansion of -2a^3 log(1+a).
    import math
    for a in (1e4, 1e6, 1e9):
        got = float(Phi_spin(a))
        ref = math.log(a) - 2.0 / 3.0 + 1.5 / a - 0.9 / a ** 2
        assert abs(got / ref - 1.0) < 1e-9, (a, got, ref)
    assert abs(float(Phi_spin(1.0)) - (1.0 - math.log(2.0))) < 1e-12
    a = 150.0
    closed = 2 * a ** 2 - a + (1 - 2 * a ** 3) * math.log1p(a) + 2 * a ** 3 * math.log(a)
    assert abs(float(Phi_spin(a)) / closed - 1.0) < 1e-9, (float(Phi_spin(a)), closed)

    # yukawa_pole must bracket its root at every density, and reach the strong
    # screening limits kappa -> pi q_tf_sq / 8 and Z -> 2.
    aL = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
    for rs in (1.0, 200.0, 1000.0, 1e5, 1e8):
        kappa, Z = yukawa_pole(rs)
        q_tf_sq = 4.0 * aL * rs / np.pi
        assert kappa > 0.0 and 1.0 < Z <= 2.0 + 1e-12, (rs, kappa, Z)
        if rs >= 1e5:
            assert abs(kappa / (np.pi * q_tf_sq / 8.0) - 1.0) < 1e-5, (rs, kappa)


def test_phi_spin_against_high_precision():
    """Phi_spin vs 90-digit evaluation of the printed form, both sides of each
    switch point included; also the exact large-a coefficients 3(n+1)/(n(n+3))."""
    from fractions import Fraction
    import mpmath as mp
    from ueg_exchange.highdensity import Phi_spin
    assert [Fraction((-1) ** (n + 1) * 3 * (n + 1), n * (n + 3)) for n in (5, 6)] == \
        [Fraction(9, 20), Fraction(-7, 18)]
    mp.mp.dps = 90
    exact = lambda a: 2 * a ** 2 - a + (1 - 2 * a ** 3) * mp.log1p(a) + 2 * a ** 3 * mp.log(a)
    pts = list(np.logspace(-12, 15, 200)) + [0.5, np.nextafter(0.5, 1), np.nextafter(2.0, 0),
                                             2.0, 30.0, 100.0, 149.999, 150.0]
    for a in pts:
        e = exact(mp.mpf(float(a)))
        rel = abs(float((mp.mpf(float(Phi_spin(a))) - e) / e))
        assert rel < 1e-14, (a, rel)


def test_exchange_shape_against_high_precision():
    """One stable shape serves f_L and both self-energy helpers.

    The reference uses the printed expression at 90 decimal digits, including
    both sides of each branch and the cusp. The old direct self-energy helper
    was negative at k=1e6 instead of approaching the positive 1/(3k^2) tail.
    """
    import mpmath as mp
    points = list(np.logspace(-12, 12, 100)) + [0.0, 1.0]
    for switch in (1e-3, 1.0, 30.0):
        points += [np.nextafter(switch, 0.0), switch,
                   np.nextafter(switch, np.inf), switch * (1 - 1e-8),
                   switch * (1 + 1e-8)]
    points = np.array(points)
    with mp.workdps(90):
        def exact(k):
            x = abs(mp.mpf(float(k)))
            if x == 0:
                return mp.mpf(1)
            if x == 1:
                return mp.mpf('0.5')
            return mp.mpf('0.5') + (1 - x*x)/(4*x)*mp.log(abs((1+x)/(1-x)))
        refs = [exact(k) for k in points]
        for fn in (L.f_L, L.exchange_self_energy_shape, H.sigma_x_shape, H.f_L_static):
            vals = fn(points.reshape(3, -1))
            assert vals.shape == (3, len(points) // 3)
            vals = vals.ravel()
            for k, got, want in zip(points, vals, refs):
                assert got > 0, (fn.__name__, k, got)
                assert abs((mp.mpf(float(got)) - want) / want) < 2e-12, (fn.__name__, k, got, want)
                assert float(fn(-k)) == got, (fn.__name__, k)
            assert float(fn(np.inf)) == 0.0


def test_dg_dnu_large_q_matches_integral():
    """dg/dnu (external/benites_direct.py) against a 50-digit integral, on both sides of the
    q = 3 switch and at large q with small nu, where the closed form cancels."""
    import mpmath as mp
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "external"))
    import benites_direct as bd
    mp.mp.dps = 50
    def ref(q, nu):
        q, nu = mp.mpf(q), mp.mpf(nu); a = q * q / 2
        return -8 * a * nu * mp.quad(lambda k: k * k / (((q * k + a) ** 2 + nu * nu) * ((q * k - a) ** 2 + nu * nu)), [0, 1])
    for q, nu in [(2.999999, 0.7), (3.0, 0.7), (3.000001, 0.7), (30.0, 1e-3), (1e3, 1e-6), (1e4, 1e-3), (1e4, 5.0)]:
        got = float(bd.dg_dnu_brm(np.array([q]), np.array([nu]))[0])
        want = float(ref(q, nu))
        assert got * want > 0 and abs(got - want) <= 1e-12 * abs(want), (q, nu, got, want)


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"  PASS  {name}")
            except Exception as exc:  # noqa
                fails += 1
                print(f"  FAIL  {name}: {exc}")
    print("all passed" if not fails else f"{fails} failure(s)")
    sys.exit(1 if fails else 0)
