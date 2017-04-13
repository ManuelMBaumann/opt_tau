#! /bin/bash
. ~/.bashrc


#Exp 2
python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,10.0] --Nom=10 --damping=0.05 --dg_pp=0 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp2_dgpp0.txt                    
python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,10.0] --Nom=10 --damping=0.05 --dg_pp=1 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp2_dgpp1.txt     
python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,10.0] --Nom=10 --damping=0.05 --dg_pp=2 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp2_dgpp2.txt     
python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,10.0] --Nom=10 --damping=0.05 --dg_pp=3 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp2_dgpp3.txt     
python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,10.0] --Nom=10 --damping=0.05 --dg_pp=4 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp2_dgpp4.txt     
python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,10.0] --Nom=10 --damping=0.05 --dg_pp=5 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp2_dgpp5.txt     
python3 -u elast_wedge.py --tau_re=-100 --dx=5 --dz=5 --freq=[1.0,10.0] --Nom=10 --damping=0.05 --dg_pp=10 --block=True --plots=False --solver_flag=0 --nprocs=8 | tee experm/exp2_dgpp10.txt     
                            
                            
                            
                            