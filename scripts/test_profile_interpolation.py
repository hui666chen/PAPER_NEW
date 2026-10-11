"""Regression checks against NumPy's complex endpoint-clamped interpolation."""
import sys
import unittest
import numpy as np

sys.path.insert(0, '<LOCAL_USER>/chh/GPU/High')
from audit_isolated_delay_profiles import interp_uniform


class BoundaryTests(unittest.TestCase):
    def test_endpoints_and_outside(self):
        y = np.array([1+2j, 4-3j, 9+7j])
        q = np.array([-100, -0.1, 0, 2, 2.1, 100])
        np.testing.assert_array_equal(interp_uniform(y, 0, 1, q), y[[0,0,0,2,2,2]])

    def test_random_batched_complex_and_knots(self):
        rng = np.random.default_rng(42)
        x = np.linspace(-2, 2, 401)
        y = rng.normal(size=401) + 1j*rng.normal(size=401)
        q = np.concatenate((x, rng.uniform(-5,5,10000)))[None, :]
        np.testing.assert_allclose(interp_uniform(y,-2,.01,q), np.interp(q,x,y), atol=1e-11, rtol=1e-11)


if __name__ == '__main__':
    unittest.main()
