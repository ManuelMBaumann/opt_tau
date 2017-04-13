import os
import numpy as np
import matplotlib.pyplot as plt
from plot_misc import opt_tau_anal
from math import pi

out_file = 'experm/it_tmp.txt'
cc = list('bgrcmy')

dg_pp = 5
freq = np.linspace(1.0, 5.0, 5, endpoint=True)
dampings = [0.5, 0.3, 0.1, 0.05]
n_exp = 30
#dampings = [0.8, 0.5]
#n_exp = 2
    
tau_im     = np.empty((len(dampings),n_exp+1))
tau_im_ind = np.empty(len(dampings), dtype='int')
iters      = np.empty((len(dampings),n_exp+1))
kk         = 0
    
for damping in dampings:
     
    om           = 2.0*pi*freq*(1.0-1j*damping)    
    tau_anal     = opt_tau_anal(damping, min(om.real), max(om.real))
    tau_re       = tau_anal.real
    tau_im_tmp   = -np.linspace(0.01,1.0,n_exp)*max(om.real)    
    tau_im[kk,:] = np.sort(np.concatenate([tau_im_tmp,[tau_anal.imag]]))

    # find index of analytic tau
    for ii in range(len(tau_im[kk,:])):
        if tau_im[kk,ii] > tau_anal.imag:
            tau_im_ind[kk] = ii-1
            break
    
    for ii in range(n_exp+1):
   
        damp_str   = ' --damping='+str(damping)
        tau_re_str = ' --tau_re='+str(tau_re) 
        tau_im_str = ' --tau_im='+str(tau_im[kk,ii]) 
        dg_pp_str  = ' --dg_pp='+str(dg_pp) 
       
        str1 = 'python3 elast_wedge.py --ndims=2 --dx=10.0 --dz=10.0 --freq=[1,5] --Nom=5 --block=True --solver_flag=0' 
        os.system(str1+damp_str+tau_re_str+tau_im_str+dg_pp_str+' --nprocs=8')

        with open(out_file, "r") as f:
            iters[kk,ii] = int( f.read() )
    
    kk = kk+1 
    
    
for kk in range(len(dampings)):
    plt.plot( tau_im[kk,:]/max(om.real), iters[kk,:], cc[kk], label=r'$\epsilon$ = '+str(dampings[kk]), linewidth=2.0 )           
    plt.plot(tau_im[kk,tau_im_ind[kk]]/max(om.real), iters[kk,tau_im_ind[kk]], cc[kk]+'x' , markersize=15, mew=2 )
    #plt.title('Optimally preconditioned multi-shift GMRES')
    plt.xlabel('Seed shift (scaled imag part)', fontsize=20)
    plt.ylabel('Number of iterations', fontsize=20)
    plt.ylim((0,np.max(iters)+10))
    plt.legend(fontsize=20)
    
plt.savefig('experm/exp2_dp'+str(int(dg_pp))+'.png')
plt.savefig('experm/exp2_dp'+str(int(dg_pp))+'.pdf')
plt.show()