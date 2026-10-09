"""Growth recovery, forecast meaning, and endpoint behavior."""
import numpy as np
import pytest
import epidemic_modeling as em


@pytest.mark.parametrize('causal,wlen',[(True,7),(False,8)])
def test_log_fit_known_exponential(causal,wlen):
    x=25*np.exp(.035*np.arange(40)); rt,a,g,fit=em.rt_exp_fit_log_lin_reg(x,wlen,1,causal)
    indices=slice(wlen-1,None) if causal else slice(wlen//2,-wlen//2)
    np.testing.assert_allclose(g[indices],.035,atol=1e-12)
    np.testing.assert_allclose(a[indices],x[indices])
    np.testing.assert_allclose(fit[indices],x[indices]*np.exp(.035))


def test_generation_ratio_known_exponential():
    x=25*np.exp(.035*np.arange(40)); rt,g,rs,gs=em.rt_exp_fit_gen_ratios(x,7,5,1)
    np.testing.assert_array_equal(g[:5],0)
    np.testing.assert_allclose(g[5:],.035)
    np.testing.assert_allclose(gs[11:],.035)


def test_nonlinear_fit_and_zero_windows():
    x=25*np.exp(.035*np.arange(40)); rt,a,g,fit=em.rt_exp_fit_nonlin_ls(x,7,1)
    np.testing.assert_allclose(g[6:],.035,atol=1e-6)
    np.testing.assert_allclose(a[6:],x[6:],rtol=1e-6)
    assert np.all(a[:6]==0)
    x[10]=0; result=em.rt_exp_fit_nonlin_ls(x,7,1)
    assert result[2][10]==0 and result[1][10]==0
    centered=em.rt_exp_fit_nonlin_ls(x,8,1,False)
    assert centered[2][9]==0 and centered[1][9]==x[9]
    with pytest.raises(ValueError): em.rt_exp_fit_log_lin_reg(x,7,1)
