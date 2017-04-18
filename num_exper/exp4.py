import os
import numpy as np

np   = 4
fmin = 1.0
fmax = 9.0
df   = 0.5

for i in np.arange(fmin+df, fmax, df): 
    
              
    str1 = 'python3 elast_wedge.py --ndims=2 --dx=100.0 --dy=100.0 --dz=100.0 --df=0.1 --degree=1 --damping=0.5 --maxit=300 \
                                   --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                                   --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8'
               
    
    w_tick3 = np.logspace(np.log10(i), np.log10(fmax), num=np, endpoint=True)
    int1 = ' --freq=[1,'+str(i)+']'
    int2 = ' --freq=['+str(i)+','+str(w_tick3[1])+']'
    int3 = ' --freq=['+str(w_tick3[1])+','+str(w_tick3[2])+']'
    int4 = ' --freq=['+str(str(w_tick3[2]))+',9]'
    
    os.system(str1+int1)
    os.system(str1+int2)
    os.system(str1+int3)
    os.system(str1+int4)                  