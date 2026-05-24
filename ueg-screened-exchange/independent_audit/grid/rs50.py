import numpy as np
from hd_split import energy_split
from math import sqrt
from a7_mc import aTF
rs=50.0
st=sqrt(aTF*rs); qs=max(1.0,(4*aTF*rs/3)**0.25)
br=(0.0,min(st,0.5),max(2*qs,2.0),None); ns=(20,36,28)
vals=[]
for i in range(8):
    v=1000*energy_split(rs,br,ns,nnu=40,N=200_000,seed=3000+i)
    vals.append(v); print(f'  replicate {i}: {v:.3f}',flush=True)
m=np.mean(vals); se=np.std(vals,ddof=1)/np.sqrt(len(vals))
print(f'\n  rs=50   this work {m:.3f} +- {se:.3f} mRy   PRB direct -19.854(0.122)   BRM -19.45(0.65)')
print(f'  vs PRB {(m+19.854)/np.hypot(se,0.122):+.1f} sigma   vs BRM {(m+19.45)/np.hypot(se,0.65):+.1f} sigma')
