"""Physical identities for the common v0 / inverse-v0 notation convention."""
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ueg_exchange import lindhard as L, constants as C, highdensity as H


def test_screened_line_and_coupling_average():
    for rs in (0.0, 0.1, 5.0, 1000.0):
        for zeta in (0.0, 0.5, 1.0):
            for q, nu in ((0.1, 0.0), (0.7, 0.4), (4.0, 3.0)):
                v = 1.0 / q ** 2
                assert L.screened_interaction(q, nu, rs, 0.0, zeta) == 0.0
                eps = L.rpa_dielectric(q, nu, rs, 0.3, zeta)
                Pi0 = L.noninteracting_response(q, nu, rs, zeta)
                chi = L.rpa_density_response(q, nu, rs, 0.3, zeta)
                W = L.screened_interaction(q, nu, rs, 0.3, zeta)
                assert np.isclose(eps, 1 - 0.3 * v * Pi0, rtol=2e-15)
                assert np.isclose(chi * eps, Pi0, rtol=2e-15)
                assert np.isclose(W, 0.3 * v + (0.3 * v) ** 2 * chi,
                                  rtol=2e-10, atol=1e-12)
                average = quad(lambda coupling: L.screened_interaction(
                    q, nu, rs, coupling, zeta), 0, 1, epsabs=1e-12, epsrel=1e-11)[0]
                assert np.isclose(L.coupling_averaged_interaction(q, nu, rs, zeta),
                                  average, rtol=2e-10, atol=1e-12)
                if rs == 0:
                    assert L.coupling_averaged_interaction(q, nu, rs, zeta) == v / 2


def test_interaction_limits_and_polarized_free_response():
    for zeta in (0.0, 0.7, 1.0):
        assert np.isclose(L.lindhard_shape_spin(0.0, 0.0, zeta), C.spin_moment(1, zeta))
        assert L.lindhard_shape_spin(0.0, 1.0, zeta) == 0.0
        for coupling in (0.1, 0.5, 1.0):
            assert np.isclose(L.screened_interaction(0.7, 1e7, 100.0, coupling, zeta),
                              coupling / 0.7 ** 2, rtol=1e-11)
        assert np.isclose(L.screened_interaction(0.7, 0.4, 2.0, 1e-12, zeta) / 1e-12,
                          1 / 0.7 ** 2, rtol=1e-11)


def test_common_response_normalization_preserves_energy():
    for rs in (0.1, 1.0, 10.0):
        screening_sq = float(C.q_tf_sq(rs))
        for z in (0.0, 0.3, 1.0, 3.0):
            vtx, se = H.response_log_coefficients(z, rs)
            legacy_vertex = -screening_sq * H.alpha_shell(z) / (4 * np.pi ** 2)
            assert np.isclose(vtx, np.pi ** 2 * screening_sq * legacy_vertex)
            assert vtx + se == 0.0
            old_closure = -3 / (np.pi * screening_sq) * legacy_vertex
            new_closure = -3 / (np.pi ** 3 * screening_sq ** 2) * vtx
            assert np.isclose(new_closure, old_closure, rtol=3e-15)
        new_vertex_prefactor = -screening_sq ** 2 / (32 * np.pi ** 2)
        assert np.isclose(-3 / (np.pi ** 3 * screening_sq ** 2) * new_vertex_prefactor,
                          C.A0, rtol=3e-15)


def test_brm_derivative_at_kohn_edge():
    """g_BRM=2L has the right derivative -pi/4 at exactly q=2,nu=0."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'external'))
    from benites_direct import dg_dnu_brm
    edge = float(dg_dnu_brm(2.0, 0.0))
    assert edge == -np.pi / 4.0
    values = np.asarray([dg_dnu_brm(2.0, nu) for nu in (1e-4, 1e-6, 1e-8)])
    assert np.all(np.diff(np.abs(values-edge)) < 0)
    assert abs(values[-1]-edge) < 1e-7


if __name__ == '__main__':
    for name, function in sorted(globals().copy().items()):
        if name.startswith('test_') and callable(function):
            function()
            print('PASS', name)
