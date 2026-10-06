"""Hand-written analytical API example; no biological/simulation claim."""
import numpy as np
from pycoloc import coloc_abf, plot_coloc
d1=dict(snp=['a','b','c'],position=[100,200,300],chromosome='22',
    beta=[.1,.5,-.2],varbeta=[.01,.02,.03],type='cc')
d2=dict(d1,beta=[0.,.4,.1],type='quant',sdY=1.2)
result=coloc_abf(d1,d2,p1=1e-4,p2=1e-4,p12=1e-5,
    prior_weights1={'a':1.,'b':2.,'c':0.})
assert np.isclose(result.results['SNP.PP.H4'].sum(),1.)
assert result.results.loc[result.results.snp=='c','SNP.PP.H4'].item()==0.
fig,axes=plot_coloc(result,d1,d2)
assert len(axes)==2
print(result.summary.to_string(index=False))
