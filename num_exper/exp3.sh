#! /bin/bash
. ~/.bashrc


#Exp 3
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=5 --degree=1 --damping=0.5 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=False \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=False --nprocs=8  | tee experm/exp3_f510_Nom5_e05_norot_C0.txt
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=15 --degree=1 --damping=0.5 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=False \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=False --nprocs=8  | tee experm/exp3_f510_Nom15_e05_norot_C0.txt  

python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=5 --degree=1 --damping=0.5 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=False \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=True --nprocs=8  | tee experm/exp3_f510_Nom5_e05_rot_C0.txt
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=15 --degree=1 --damping=0.5 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=False \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=True --nprocs=8  | tee experm/exp3_f510_Nom15_e05_rot_C0.txt                           
                          
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=5 --degree=1 --damping=0.05 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=False \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=False --nprocs=8  | tee experm/exp3_f510_Nom5_e005_norot_C0.txt
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=15 --degree=1 --damping=0.05 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=False \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=False --nprocs=8  | tee experm/exp3_f510_Nom15_e005_norot_C0.txt  

python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=5 --degree=1 --damping=0.05 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=False \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=True --nprocs=8  | tee experm/exp3_f510_Nom5_e005_rot_C0.txt
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=15 --degree=1 --damping=0.05 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=False \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=True --nprocs=8  | tee experm/exp3_f510_Nom15_e005_rot_C0.txt  





python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=5 --degree=1 --damping=0.5 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=False --nprocs=8  | tee experm/exp3_f510_Nom5_e05_norot.txt
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=15 --degree=1 --damping=0.5 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=False --nprocs=8  | tee experm/exp3_f510_Nom15_e05_norot.txt  

python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=5 --degree=1 --damping=0.5 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=True --nprocs=8  | tee experm/exp3_f510_Nom5_e05_rot.txt
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=15 --degree=1 --damping=0.5 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=True --nprocs=8  | tee experm/exp3_f510_Nom15_e05_rot.txt                           
                          
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=5 --degree=1 --damping=0.05 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=False --nprocs=8  | tee experm/exp3_f510_Nom5_e005_norot.txt
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=15 --degree=1 --damping=0.05 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=False --nprocs=8  | tee experm/exp3_f510_Nom15_e005_norot.txt  

python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=5 --degree=1 --damping=0.05 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=True --nprocs=8  | tee experm/exp3_f510_Nom5_e005_rot.txt
python3 -u elast_wedge.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[5,10] --Nom=15 --degree=1 --damping=0.05 --maxit=1 \
                          --tol=1e-8 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 --block=True \
                          --plots=False --plot_resnrm=true --plot_ritz=False --solver_flag=1 --rot=True --nprocs=8  | tee experm/exp3_f510_Nom15_e005_rot.txt                            
                            
                            
                            
                            
                            
                            
                            
                            
                            