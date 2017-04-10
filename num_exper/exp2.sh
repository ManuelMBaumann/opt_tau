#! /bin/bash
. ~/.bashrc


#Exp 2
python3 -u elast_wedge.py --ndims=2 --dx=5.0 --dy=5.0 --dz=5.0 --freq=[1,10] --Nom=10 --degree=1 --damping=0.07 --maxit=300 \
                            --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                            --plots=True --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1.txt
python3 -u elast_wedge.py --ndims=2 --dx=5.0 --dy=5.0 --dz=5.0 --freq=[1,10] --Nom=10 --degree=1 --damping=0.07 --maxit=300 \
                            --tol=1e-8 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1.txt
python3 -u elast_wedge.py --ndims=2 --dx=5.0 --dy=5.0 --dz=5.0 --freq=[1,10] --Nom=10 --degree=1 --damping=0.07 --maxit=300 \
                            --tol=1e-8 --dg_pp=2 --tau_re=-100 -tau_im=-0.7 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1.txt
python3 -u elast_wedge.py --ndims=2 --dx=5.0 --dy=5.0 --dz=5.0 --freq=[1,10] --Nom=10 --degree=1 --damping=0.07 --maxit=300 \
                            --tol=1e-8 --dg_pp=3 --tau_re=-100 -tau_im=-0.7 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1.txt                            
python3 -u elast_wedge.py --ndims=2 --dx=5.0 --dy=5.0 --dz=5.0 --freq=[1,10] --Nom=10 --degree=1 --damping=0.07 --maxit=300 \
                            --tol=1e-8 --dg_pp=4 --tau_re=-100 -tau_im=-0.7 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1.txt
python3 -u elast_wedge.py --ndims=2 --dx=5.0 --dy=5.0 --dz=5.0 --freq=[1,10] --Nom=10 --degree=1 --damping=0.07 --maxit=300 \
                            --tol=1e-8 --dg_pp=5 --tau_re=-100 -tau_im=-0.7 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1.txt                            
python3 -u elast_wedge.py --ndims=2 --dx=5.0 --dy=5.0 --dz=5.0 --freq=[1,10] --Nom=10 --degree=1 --damping=0.07 --maxit=300 \
                            --tol=1e-8 --dg_pp=10 --tau_re=-100 -tau_im=-0.7 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1.txt                            
                            
                            
                            
                            
                            