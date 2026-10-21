import unittest
import csv
import json
import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np

from scripts.gss_multipath.gss import (
    Config, SampledCorrelation, estimate, golden_maximum, triangular_acf,
)


class GSSTests(unittest.TestCase):
    def test_golden_matches_known_maximum(self):
        for target in (-0.8, -0.03, 0.41, 1.3):
            result = golden_maximum(lambda x: -(x - target)**2, Config())
            self.assertTrue(result.converged)
            self.assertLessEqual(abs(result.delay - target), 0.0005)

    def test_clean_shift_arm_and_sign_invariance(self):
        for delay in (-0.1, 0.0, 0.1):
            for gain in (1, -3, 2j, -2j, 1+1j):
                result = estimate(lambda x: gain * triangular_acf(x - delay))
                self.assertLess(abs(result["estimate_chip"] - delay), 0.002)
                self.assertEqual(len(result["rounds"]), 4)
                self.assertFalse(result["pseudorange_correction_applied"])

    def test_short_single_reflection_phase_sweep(self):
        for phase in np.deg2rad([0, 45, 90, 120, 180, 270]):
            result = estimate(lambda x: triangular_acf(x)
                              + 0.2*np.exp(1j*phase)*triangular_acf(x-0.05))
            self.assertLess(abs(result["estimate_chip"]), 0.003)

    def test_three_weak_reflections(self):
        def curve(x):
            return triangular_acf(x) + sum(
                a*np.exp(1j*p)*triangular_acf(x-d)
                for a, p, d in [(0.15, 0, 0.05), (0.1, 1.2, 0.15),
                                 (0.1, 2.8, 0.4)])
        self.assertLess(abs(estimate(curve)["estimate_chip"]), 0.003)

    def test_destructive_and_over_subtraction_branches(self):
        for delay, phase, flag in [
            (0.2, 3*np.pi/4, "shifted_plus_one"),
            (0.6, np.pi, "one_chip_alias_rejected"),
        ]:
            result = estimate(lambda x: triangular_acf(x)
                              + 0.5*np.exp(1j*phase)*triangular_acf(x-delay))
            self.assertTrue(any(row[flag] for row in result["rounds"]))
            self.assertLess(abs(result["estimate_chip"]), 0.003)

    def test_strong_destructive_failure_is_not_labeled_confident(self):
        # Preserve a counterexample: prompt sign inversion can select the
        # reflection. This test documents behavior, not successful mitigation.
        result = estimate(lambda x: triangular_acf(x)-1.2*triangular_acf(x-0.1))
        self.assertGreater(abs(result["estimate_chip"]), 0.09)
        self.assertFalse(result["physical_confidence_calibrated"])

    def test_existing_csv_adapter(self):
        from scripts.gss_multipath.run_scans import main
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory)/"scan.csv", Path(directory)/"result.json"
            with source.open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(["window", "epoch", "delay_chip", "mixed_i", "mixed_q"])
                for x in np.linspace(-1, 2.5, 701):
                    writer.writerow([2, 3, x, triangular_acf(x), 0])
            with patch("sys.argv", ["run_scans", "--input", str(source),
                                   "--output", str(output), "--expected-taps", "701",
                                   "--axis", "physical-delay"]):
                main()
            result = json.loads(output.read_text())["results"][0]
            self.assertEqual((result["window"], result["epoch"]), (2, 3))
            self.assertLess(abs(result["estimate_chip"]), 0.002)

    def test_sampling_checks_and_interpolation(self):
        x = np.linspace(-1, 2.5, 701)
        correlation = SampledCorrelation(x, np.maximum(0, 1-abs(x)))
        correlation.validate_support(Config())
        self.assertLess(abs(estimate(correlation)["estimate_chip"]), 0.002)
        with self.assertRaises(ValueError):
            correlation(2.51)
        with self.assertRaises(ValueError):
            SampledCorrelation(np.linspace(-2, 2, 401), np.ones(401)).validate_support(Config())

    def test_bad_and_flat_zero_inputs(self):
        with self.assertRaises(ValueError):
            estimate(lambda x: 0j)
        with self.assertRaises(ValueError):
            estimate(lambda x: complex(float("nan")))
        with self.assertRaises(ValueError):
            estimate(lambda x: triangular_acf(x), reference=lambda x: 2.0)
        with self.assertRaises(ValueError):
            estimate(lambda x: triangular_acf(x), config=Config(max_iterations=1))


if __name__ == "__main__":
    unittest.main()
