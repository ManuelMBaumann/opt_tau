#! /bin/bash
. ~/.bashrc


# Exp. 1
# python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,5.0] --Nom=5 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp1_Nom5_f15.txt  
# python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,5.0] --Nom=10 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp1_Nom10_f15.txt  
# python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,5.0] --Nom=20 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp1_Nom20_f15.txt  
# 
# python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,10.0] --Nom=5 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp1_Nom5_f110.txt  
# python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,10.0] --Nom=10 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp1_Nom10_f110.txt  
# python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,10.0] --Nom=20 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp1_Nom20_f110.txt 
# 
# 
# 
# python3 -u elast_wedge.py --tau_re=-100 --dx=2.5 --dz=2.5 --freq=[1.0,5.0] --Nom=5 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp1_Nom5_f15_h2.txt  
# python3 -u elast_wedge.py --tau_re=-100 --dx=2.5 --dz=2.5 --freq=[1.0,5.0] --Nom=10 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp1_Nom10_f15_h2.txt  
# python3 -u elast_wedge.py --tau_re=-100 --dx=2.5 --dz=2.5 --freq=[1.0,5.0] --Nom=20 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp1_Nom20_f15_h2.txt  
# 
# 
# python3 -u elast_wedge.py --tau_re=-100 --dx=2.5 --dz=2.5 --freq=[1.0,10.0] --Nom=10 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp1_Nom10_f110_h2.txt  
                            

# python3 elast_wedge.py --tau_re=-100 --dx=10 --dz=10 --freq=[1.0,10.0] --Nom=5 --damping=0.7 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8                    
python3 elast_wedge.py --tau_re=0.3 --tau_im=-0.7 --dx=10 --dz=10 --freq=[1.0,10.0] --Nom=5 --damping=0.7 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8                                  
   
                            