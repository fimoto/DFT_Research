"""Public Lorentzian selectors, occupation sectors and independent frequency integrals."""
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'external/kite-numerics/mc'))
from ueg_exchange import lorentzian
import xc_mc


class LorentzianContract(unittest.TestCase):
    def test_sector_sum_and_bare_limit(self):
        for q, width in [(0.4, 0.2), (2.3, 1.0), (8.0, 4.0)]:
            def samples(kernel, lam=width):
                return xc_mc.xc_samples(q, 4096, np.random.default_rng(2919),
                                        kernel=kernel, lam=lam)
            total=samples('lorentzian')
            same=samples('lorentzian_same'); opposite=samples('lorentzian_opp')
            self.assertTrue(np.any(total != 0))
            np.testing.assert_array_equal(total, same + opposite)
            np.testing.assert_allclose(samples('lorentzian', 1e12), samples('static'),
                                       rtol=1e-9, atol=1e-9)

    def test_frequency_integral(self):
        for width, a, b in [(0.7, 0.8, 1.4), (0.7, 0.8, -1.4), (1.0, 1.0, 1.0)]:
            direct=quad(lambda nu: (a*b-nu**2)/((a*a+nu**2)*(b*b+nu**2))
                        * width**2/(width**2+nu**2)/np.pi,
                        0, np.inf, epsabs=1e-12, epsrel=1e-12)[0]
            closed=width**2*xc_mc._J_pole(width, np.array([a]), np.array([b]))[0]
            self.assertAlmostEqual(direct, closed, places=12)

    def test_invalid_selector_and_width_raise(self):
        with self.assertRaisesRegex(ValueError, 'Unknown kernel'):
            xc_mc.xc_samples(1.0, 16, np.random.default_rng(1), kernel='unsupported')
        for width in (None, 0.0, -1.0, np.inf, np.nan):
            with self.assertRaisesRegex(ValueError, 'positive width'):
                xc_mc.xc_samples(1.0, 16, np.random.default_rng(1),
                                 kernel='lorentzian', lam=width)

    def test_radial_energy_prefactor(self):
        q=np.array([2.0, 10.0, 100.0]); width=3.0
        expected=(3/np.pi**3)*lorentzian.p_loc(width/q)/q**2
        np.testing.assert_allclose(lorentzian.radial_energy_large_q(q,width),
                                   expected, rtol=5e-16, atol=0)


if __name__ == '__main__':
    unittest.main()
