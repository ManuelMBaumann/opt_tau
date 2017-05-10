#! /bin/bash
. ~/.bashrc


                          
python3 -u elast_wedge.py --ndims=2 --dx=5 --dz=5 --freq=[5,10] --Nom=10 --degree=1 --damping=0.0 --maxit=300 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 --tau_im=-0.7 --block=True \
                          --plots=False --plot_resnrm=True --plot_ritz=False --solver_flag=0 --rot=False --nprocs=8  | tee experm/exp5_topt.txt

python3 -u elast_wedge.py --ndims=2 --dx=5 --dz=5 --freq=[5,10] --Nom=10 --degree=1 --damping=0.0 --maxit=300 \
                          --tol=1e-8 --dg_pp=0 --tau_re=0.3 --tau_im=-0.7 --block=True \
                          --plots=False --plot_resnrm=True --plot_ritz=False --solver_flag=0 --rot=False --nprocs=8  | tee experm/exp5_t1.txt

python3 -u elast_wedge.py --ndims=2 --dx=5 --dz=5 --freq=[5,10] --Nom=10 --degree=1 --damping=0.0 --maxit=300 \
                          --tol=1e-8 --dg_pp=0 --tau_re=1.0 --tau_im=-0.5 --block=True \
                          --plots=False --plot_resnrm=True --plot_ritz=False --solver_flag=0 --rot=False --nprocs=8  | tee experm/exp5_t2.txt                          
                          
python3 -u elast_wedge.py --ndims=2 --dx=5 --dz=5 --freq=[5,10] --Nom=10 --degree=1 --damping=0.0 --maxit=300 \
                          --tol=1e-8 --dg_pp=0 --tau_re=0.91 --tau_im=-0.0035 --block=True \
                          --plots=False --plot_resnrm=True --plot_ritz=False --solver_flag=0 --rot=False --nprocs=8  | tee experm/exp5_t3.txt                            
                          
                          
                          
                          
                          
