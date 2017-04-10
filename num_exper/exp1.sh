#! /bin/bash
. ~/.bashrc


#Exp 1
# python3 -u elast_wedge.py --ndims=2 --dx=5.0 --dy=5.0 --dz=5.0 --freq=[1,5] --Nom=10 --degree=1 --damping=0.07 --maxit=300 \
#                             --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
#                             --plots=True --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1.txt

python3 -u elast_wedge.py --ndims=2 --dx=100.0 --dy=100.0 --dz=100.0 --freq=[1,3] --Nom=7 --degree=1 --damping=0.05 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=False \
                          --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1.txt
                            
                            
# python3 -u elast_wedge.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[20,26] --Nom=7 --degree=1 --damping=0.05 --maxit=30 \
#                           --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=False \
#                           --plots=True --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1.txt
                            
                            
                            