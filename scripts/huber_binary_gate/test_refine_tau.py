import unittest
import numpy as np
from refine_tau import refine

class RefineTests(unittest.TestCase):
    def test_interior(self):
        r=refine(lambda t:(t-.083456)**2,.08,.01)
        self.assertAlmostEqual(r['tau_chip'],.083456,places=6)
        self.assertLessEqual(r['train_rss_after'],r['train_rss_before'])
    def test_clipped_global_domain(self):
        r=refine(lambda t:(t-.5)**2,.4,.01)
        self.assertEqual(r['tau_chip'],.4)
        self.assertTrue(r['local_edge'])
    def test_local_not_global(self):
        r=refine(lambda t:(t-.2)**2,.08,.01)
        self.assertAlmostEqual(r['tau_chip'],.09)
    def test_existing_bank_off_grid_recovery(self):
        from run import BANK,ResponseSpec,load_response_model
        from fit import evaluate
        model=load_response_model(BANK,expected_spec=ResponseSpec(prn=9))
        scale=np.linalg.norm(model.exact([0])[0])
        def template(t):
            d=model.exact([-t,.15-t])/scale
            return d[0]+.3j*d[1]
        truth=.083456
        y=np.asarray([template(truth),(.8+.2j)*template(truth)])
        r=refine(lambda t:evaluate(y,np.ones(2),template(t)),.08,.01)
        self.assertAlmostEqual(r['tau_chip'],truth,places=6)
        self.assertLess(r['train_rss_after'],1e-10)

if __name__=='__main__':unittest.main()
