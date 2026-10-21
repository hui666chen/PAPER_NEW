import unittest
import numpy as np
from .core import Configuration, estimate


class Tests(unittest.TestCase):
    x = np.linspace(-2, 2, 4001)

    def curve(self, tau=0):
        return np.maximum(1-abs(self.x-tau), 0).astype(complex)

    def test_clean_and_shift(self):
        for tau in (-.04, 0, .04):
            r = estimate(self.x, self.curve(tau))
            self.assertTrue(r['valid'])
            self.assertAlmostEqual(r['estimate_chip'], tau, places=8)

    def test_complex_gain_invariance(self):
        z = self.curve()+.3j*self.curve(.08)
        a = estimate(self.x, z)['estimate_chip']
        b = estimate(self.x, z*3*np.exp(2.1j))['estimate_chip']
        self.assertAlmostEqual(a, b, places=10)

    def test_classical_dd_suppresses_resolved_echo(self):
        z = self.curve()+.3*self.curve(.3)
        eml = estimate(self.x, z, config=Configuration((.1,), (1.,)))
        dd = estimate(self.x, z)
        self.assertGreater(abs(eml['estimate_chip']), .001)
        self.assertLess(abs(dd['estimate_chip']), 1e-8)

    def test_short_delay_not_magically_removed(self):
        z = self.curve()+.3*self.curve(.01)
        self.assertGreater(abs(estimate(self.x, z)['estimate_chip']), .001)

    def test_invalid_support_and_flat(self):
        self.assertFalse(estimate(self.x, np.zeros_like(self.x))['valid'])
        self.assertFalse(estimate(self.x, np.ones_like(self.x))['valid'])
        x = np.linspace(-.1, .1, 21)
        self.assertFalse(estimate(x, np.ones(21))['valid'])

    def test_reject_invalid_configuration(self):
        with self.assertRaises(ValueError):
            Configuration((.1,), (1., 2.))


if __name__ == '__main__':
    unittest.main()
