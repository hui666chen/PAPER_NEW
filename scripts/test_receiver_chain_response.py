import unittest

import numpy as np
from scipy.signal import fftconvolve

from receiver_chain_response_control import ChainResponse, H, KMIN, KMAX, OMEGA
from run_prn14_bandlimited_interval import delayed_block
from audit_isolated_delay_profiles import interp_uniform, profile_for_window


class ChainTests(unittest.TestCase):
    def test_integer_bank_and_fractional_response_against_direct_samples(self):
        rng = np.random.default_rng(140)
        raw = rng.normal(size=(1800,2))
        sig = raw[:,0]+1j*raw[:,1]
        start, n = 600, 350
        w = rng.normal(size=(n,4))+1j*rng.normal(size=(n,4))
        extended = sig[start-KMAX:start+n-KMIN]
        bank = fftconvolve(extended[:,None],w[::-1],mode='valid',axes=0)[::-1]
        direct = np.array([sig[start-k:start+n-k]@w for k in range(KMIN,KMAX+1)])
        np.testing.assert_allclose(bank,direct,atol=1e-11)
        model = ChainResponse(bank)
        for d in [-1.2, -.06138, 0, .28644, .301785, .3069, .32736, 2.4]:
            samples = delayed_block(raw,start,n,d/H)
            expected = np.exp(1j*OMEGA*d/H)*(samples@w)
            np.testing.assert_allclose(model.exact([d])[0],expected,atol=1e-10)

    def test_cached_queries(self):
        rng = np.random.default_rng(4)
        model = ChainResponse(rng.normal(size=(KMAX-KMIN+1,4))+1j*rng.normal(size=(KMAX-KMIN+1,4)))
        model.prepare()
        d = np.array([-1.2,0,.0005,.3,2.4])
        np.testing.assert_allclose(model.response(d),model.exact(d),atol=1e-11)
        d = np.array([.301785,.301285])
        np.testing.assert_allclose(model.response(d),model.exact(d),atol=1e-11)

    def test_complex_ls_matches_lstsq(self):
        rng=np.random.default_rng(5)
        model=ChainResponse(rng.normal(size=(KMAX-KMIN+1,11))+1j*rng.normal(size=(KMAX-KMIN+1,11)))
        x=np.linspace(-2,2,11);y=rng.normal(size=11)+1j*rng.normal(size=11)
        seps=np.array([.3,.4]);taus=np.array([-.01,0,.01])
        p,t,_=model.profile(x,None,y,seps,taus)
        for i,s in enumerate(seps):
            costs=[]
            for tau in taus:
                a=np.column_stack((model.exact([-tau])[0],model.exact([s-tau])[0],np.ones_like(x),x))
                r=y-a@np.linalg.lstsq(a,y,rcond=None)[0]
                costs.append(np.vdot(r,r).real)
            self.assertAlmostEqual(p[i],min(costs),places=10)
            self.assertEqual(t[i],taus[np.argmin(costs)])

    def test_only_response_columns_change(self):
        rng=np.random.default_rng(6)
        x=np.linspace(-2,2,401)
        primary=rng.normal(size=401)+1j*rng.normal(size=401)
        y=rng.normal(size=401)+1j*rng.normal(size=401)
        model=ChainResponse(None)
        model.response=lambda d:interp_uniform(primary,x[0],np.median(np.diff(x)),x[None,:]+np.asarray(d)[:,None])
        sep=np.array([.04,.3,1.2]);tau=np.array([-1.2,-.005,0,.005,1.2])
        old=profile_for_window(x,primary,y,sep,tau)
        new=model.profile(x,primary,y,sep,tau)
        for a,b in zip(old,new):np.testing.assert_allclose(a,b,rtol=1e-12,atol=1e-12)


if __name__=='__main__':unittest.main()
