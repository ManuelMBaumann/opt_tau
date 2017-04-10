#! /bin/bash
. ~/.bashrc

(for fmid in 1.5 `seq 2 0.5 8.5`; do 
     python3 elast_wedge.py --ndims=2 --dx=100.0 --dy=100.0 --dz=100.0 --freq=[1,$fmid] --df=0.1 --degree=1 --damping=0.5 --maxit=300 \
                            --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8

     python3 elast_wedge.py --ndims=2 --dx=100.0 --dy=100.0 --dz=100.0 --freq=[$fmid,9] --df=0.1 --degree=1 --damping=0.5 --maxit=300 \
                            --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8
done) 


python3 elast_wedge.py --ndims=2 --dx=100.0 --dy=100.0 --dz=100.0 --freq=[1,9] --df=0.1 --degree=1 --damping=0.5 --maxit=300 \
                       --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                       --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8