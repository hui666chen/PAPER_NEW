import unittest
import numpy as np
from fit import Grid, gain_residual, huber


class FitTest(unittest.TestCase):
    def test_grid_matches_direct_lstsq(self):
        rng = np.random.default_rng(7)
        d = rng.normal(size=(17, 5))+1j*rng.normal(size=(17, 5))
        y = rng.normal(size=(4, 17))+1j*rng.normal(size=(4, 17))
        weights = rng.uniform(.3, 2, 4)
        pixel = rng.uniform(.2, 1, y.shape)
        grid = Grid(d, [0, 1, 2], [2, 3, 4], [0, .2+.3j, -.5j])
        scores = [np.sum(weights[:, None]*pixel*
                         abs(gain_residual(y, grid.template(k), pixel))**2)
                  for k in range(9)]
        self.assertEqual(grid.search(y, weights, pixel), np.argmin(scores))

    def test_huber_monotone(self):
        rng = np.random.default_rng(8)
        d = rng.normal(size=(21, 4))+1j*rng.normal(size=(21, 4))
        grid = Grid(d, [0, 1], [2, 3], [0, .3j])
        y = np.tile(grid.template(2), (8, 1))+.02*rng.normal(size=(8, 21))
        y[0, 5] += 20
        fit = grid.fit(y, np.ones(8), .1)
        self.assertTrue(np.all(np.diff(fit["objective"]) <= 1e-8))
        self.assertEqual(fit["index"], 2)
        self.assertAlmostEqual(float(huber(np.array(2.), 1)), 3)

    def test_amplitude_veto_after_persistence(self):
        from score import decisions
        records = []
        for w, rho in enumerate([.2, .2, .2, .2, .5, .4]):
            f = dict(rho=rho, converged=True, boundary=False, gain=.1, tau_chip=.1)
            records.append(dict(window=w, huber=[dict(f), dict(f)]))
        rows = decisions(records, "huber_cap", 0)
        self.assertFalse(rows[4]["active"])
        self.assertEqual(rows[4]["correction_m"], 0.)
        self.assertTrue(rows[5]["active"])


if __name__ == "__main__":
    unittest.main()
