"""Independent MC against BRM Table IX (zeta=0) and PRB Table II direct column."""
import numpy as np
from hd_split import energy_split
from math import sqrt,pi
from a7_mc import aTF
BRM={0.1:(43.53,0.18),1:(27.62,0.84),2:(20.31,0.33),3:(14.85,0.68),5:(8.04,0.34),
     7:(2.23,0.41),10:(-2.11,0.41),20:(-11.13,0.90),30:(-14.85,0.28),50:(-19.45,0.65)}
PRB={0.1:(43.424,0.017),1:(27.828,0.045),2:(19.515,0.058),3:(14.011,0.065),5:(6.715,0.076),
     7:(1.886,0.083),10:(-3.072,0.091),20:(-11.669,0.106),30:(-15.789,0.114),50:(-19.854,0.122)}
print('  rs     this work (mRy)     PRB direct        BRM Table IX     (mine-PRB)/sig  (mine-BRM)/sig')
for rs in [0.1,1,2,3,5,7,10,20,30,50]:
    st=sqrt(aTF*rs); qs=max(1.0,(4*aTF*rs/3)**0.25)
    br=(0.0,min(st,0.5),max(2*qs,2.0),None); ns=(20,36,28)
    vals=[1000*energy_split(rs,br,ns,nnu=40,N=250_000,seed=3000+i) for i in range(8)]
    m=np.mean(vals); se=np.std(vals,ddof=1)/np.sqrt(8)
    p,pe=PRB[rs]; b,be=BRM[rs]
    print(f'  {rs:<6} {m:8.3f} +- {se:5.3f}   {p:8.3f}({pe:.3f})  {b:8.2f}({be:.2f})   '
          f'{(m-p)/np.hypot(se,pe):+6.1f}        {(m-b)/np.hypot(se,be):+6.1f}',flush=True)
