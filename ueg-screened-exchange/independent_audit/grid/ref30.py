import numpy as np
from hd_split import energy_split
from math import sqrt
from a7_mc import aTF
rs=30.0; st=sqrt(aTF*rs); qs=max(1.0,(4*aTF*rs/3)**0.25)
for tag,br,ns in [('as used',(0.0,min(st,0.5),max(2*qs,2.0),None),(20,36,28)),
                  ('finer/wider',(0.0,0.2,1.0,3.0,8.0,None),(20,28,36,36,40))]:
    v=[1000*energy_split(rs,br,ns,nnu=56,N=200_000,seed=7000+i) for i in range(6)]
    print(f'  {tag:<12} {np.mean(v):8.3f} +- {np.std(v,ddof=1)/np.sqrt(6):5.3f}  (PRB -15.789)',flush=True)
