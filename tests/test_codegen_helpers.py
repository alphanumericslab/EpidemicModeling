"""Standalone Coder helpers agree with the six-state analytic model."""
from types import SimpleNamespace
import numpy as np
from epidemic_modeling import codegen, default_si_params
from epidemic_modeling.kalman import _si_handles


def test_standalone_helper_protocol():
    p=default_si_params([3,3]); p.update(a=np.array([.06,.04]),b=.02)
    handles=_si_handles(controlled=True,nonnegative_tie=True)
    np.testing.assert_array_equal(codegen.obs_hard_margins(np.array([-.1]),p),[-.1])
    s=np.array([.8,.1,.3,.05,.02,.03]); u=np.array([1.,2.]); w=np.zeros(6)
    for name,args in [('state_hard_margins',(s,)),
        ('nlin_state_update',(u,s,w)),('nlin_obs_update',(u,s,0)),('state_jacobians',(u,s,w)),
        ('obs_jacobian',(u,s,0)),('state_hessian_terms',(u,s,np.eye(6),w,np.eye(6))),
        ('obs_hessian_terms',(u,s,np.eye(6),0,1))]:
        result=getattr(codegen,name)(*args,p)
        reference=getattr(handles,name)(*args,SimpleNamespace(**p),0)
        if isinstance(result,tuple):
            for actual,expected in zip(result,reference): np.testing.assert_allclose(actual,expected)
        else: np.testing.assert_allclose(result,reference)
