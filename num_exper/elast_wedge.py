from __future__ import print_function, division
import matplotlib
matplotlib.use('agg')
from os import path, remove
from nutils import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.io import mmwrite
from mskrylov import poly_driver
from mekrylov import me_driver
import time
from math import pi
import scipy.sparse.linalg as spla
from plot_misc import plot_opt_splitting2

out_file = 'experm/split_f19_np4.txt'
ppw = 20

def test_orig_problem(K, C, M, b, freq, damping, x, solver_flag):
    Nom    = len(freq)
    trace_C = np.sum(C.diagonal())
    
    if abs(trace_C) > 0 & solver_flag==0:
        om = (1.0-1j*damping)*2.0*np.pi*freq         # 'our' damping model
    else:
        om2 = (1.0-1j*damping)*(2.0*np.pi*freq)**2   # Mulders' damping model
        om = np.sqrt(om2)                             
    
    normb  = np.linalg.norm(b)
    relerr = np.zeros((Nom,))
    for j in range(Nom):
        xs = x[:,j]
        r  = b - (K*xs + 1j*om[j]*(C*xs) - om[j]**2*(M*xs))
        relerr[j] = np.linalg.norm(r)/normb
    print('Relative residual of original problem:' +str(relerr))
    

@log.title
def makeplots( domain, geom, Lx, Lz, value, name, title, ndigits=0, index=None, clim=None, lineOn=False, imgtype=None,):
  points, colors = domain.elem_eval( [ geom, value ], ischeme='bezier3', separate=True )

  with plot.PyPlot( name, ndigits=ndigits, figsize=(5,6), index=index, imgtype=imgtype ) as plt:
    plt.mesh( points, colors, triangulate='bezier', edgecolors='none' )
    plt.title(title)
    plt.xlabel('x [m]')
    plt.ylabel('z [m]')

    plt.xticks( [0, Lx/2.0, Lx], ['0', '300', '600'] )
    plt.yticks( [-0, -0.4*Lz, -0.8*Lz, -Lz], ['0', '400', '800', '1000'] )

    if clim is not None:
      plt.clim(*clim)
    plt.colorbar()

    if lineOn:
      # Only for wedge problem in 2D
      plt.plot( [0, 600],[-400, -500],'k' )
      plt.plot( [0, 600],[-800, -600],'k' )

def makevtk(domain, geom, rho, lam, mu, cp, cs, sol, freq, vec_basis, name):
  Nom = len(freq)
  vtk_geom, vtk_rho, vtk_lam, vtk_mu, vtk_cp, vtk_cs = domain.simplex.elem_eval( [ geom, rho, lam, mu, cp, cs ], ischeme='vtk', separate=True )
  with plot.VTKFile( name ) as vtk:
      vtk.unstructuredgrid( vtk_geom )
      vtk.pointdataarray( 'rho', vtk_rho )
      vtk.pointdataarray( 'lambda', vtk_lam )
      vtk.pointdataarray( 'mu', vtk_mu )
      vtk.pointdataarray( 'cp', vtk_cp )
      vtk.pointdataarray( 'cs', vtk_cs )
      for k in range(0,Nom):
          disp = vec_basis.dot( sol[k,:] ).real
          vtk_disp = domain.simplex.elem_eval( disp, ischeme='vtk', separate=True )
          vtk.pointdataarray( 'disp_f'+str(freq[k]), vtk_disp )

def makespyplot( matrix, name, imgtype=None ):
  if not scipy.sparse.isspmatrix( matrix ):
      matrix = matrix.toscipy()
  with plot.PyPlot( name, ndigits=0, imgtype=imgtype ) as plt:
    plt.spy( matrix, markersize=0.8, color='blue')
    plt.title( name+', nnz = '+str(matrix.nnz) )

def point_eval(func, domain, geom, point):
  domain = domain[tuple(slice(0, p) if p > 0 else slice(None) for p in point)]
  for p in point:
      domain = domain.boundary['right' if p > 0 else 'left']
  return numpy.asarray(domain.integrate( func, geometry=geom, ischeme='gauss2' ).toscipy().todense())

def elast_mat(rho, cp, cs, lam, mu, ndims, nx, ny, nz, vec_basis, domain, geom, block):
  # define PDE
  stress = lambda u: lam*u.div(geom)[:,_,_]*function.eye(ndims) + 2.0*mu*u.symgrad(geom)
  elasticity = function.outer( stress(vec_basis), vec_basis.grad(geom) ).sum([2,3])

  w_mass = lambda u: rho*u
  mass = function.outer( w_mass(vec_basis), vec_basis ).sum(-1)

  # define Sommerfeld BC
  n = geom.normal()
  t = np.eye(ndims)
  t = t-(t*n[_,:]).sum(1)
  B_bc = cp*n[:,_]*n[_,:]+cs*(t[:,:,_]*t[:,_,:]).sum(0)

  bc_fun = lambda u: rho*(B_bc*u[:,_,:]).sum(-1)
  sommerfeld = function.outer( bc_fun(vec_basis), vec_basis ).sum(-1)

  if ndims == 2:
      sommerfeld_boundary = 'left,right,bottom'
      source_position = nx//2, nz
  else:
      sommerfeld_boundary = 'left,right,bottom,front,back'
      source_position = nx//2, ny//2, nz

  # Build matrices
  K, M  = domain.integrate( [elasticity, mass], geometry=geom, ischeme='gauss2' )
  C     = domain.boundary[sommerfeld_boundary].integrate( sommerfeld, geometry=geom, ischeme='gauss2' )
  
  # Build RHS
  if not block:
      C = 0.0*C
      source_position = nx//2, nz//2
      rhs = point_eval(vec_basis, domain, geom, source_position)[:,-1] #+ point_eval(vec_basis, domain, geom, source_position)[:,0] 
  else:
      rhs = point_eval(vec_basis, domain, geom, source_position)[:,-1] 
             
  return K, C, M, rhs

    
def main( ndims=2,           # problem dimension (2,3) 
          dx=100.0,          # grid size in x-direction 
          dy=100.0,          # grid size in y-direction          
          dz=100.0,          # grid size in z-direction  
          freq=[1.0,9.0],    # frequencies in Hz 
          Nom=7,             # number of freq's
          #df=.5,            # number of equally-spaced freq's
          degree=1,          # degree of FEM splines
          damping=0.7,       # viscous damping param    
          maxit=300,         # max no of iterations
          tol=1e-8,          # residual norm tolerance
          dg_pp=0,           # degree of poly preconditioner
          rot = False,       # rotation in MEqn approach
          tau_re=0.7,        # real(seed), if tau.real<0: take 'optimal' tau
          tau_im=-0.3,       # imag(seed)
          block=True,        # C=0 if False
          plots=False,       # plots on/off
          plot_resnrm=False, # display residual norm live
          plot_ritz=False,   # plot ritz values
          iLU=False,
          fill_factor=10,
          solver_flag=0):    # -1(python's built-in), 0(poly_pre), 1(matr_eqn)

  tau = tau_re+1j*tau_im

  # domain size
  Lx = 600.0
  Ly = 600.0
  Lz = 1000.0

  # problem parameters
  freq = np.linspace(freq[0],freq[-1],Nom)
  #freq = np.arange(freq[0], freq[-1]+df, df)
  Nom  = len(freq)

  # define physical params
  rho0 = 1800.0
  rho1 = 2100.0
  rho2 = 1950.0

  cp0  = 2000.0
  cp1  = 3000.0
  cp2  = 2300.0

  cs0  = 800.0
  cs1  = 1600.0
  cs2  = 1100.0

  # define Cartesian grid
  nx = int(np.round(Lx/dx))+1
  nz = int(np.round(Lz/dz))+1
  verts_x = np.linspace( 0, Lx, nx )
  verts_z = np.linspace( -Lz, 0, nz )

  if ndims == 2:
      ny = 1
      dy = 0.
      verts = [verts_x, verts_z]
  elif ndims == 3:
      ny = int(np.round(Ly/dy))+1
      verts_y = np.linspace( 0, Ly, ny )
      verts = [verts_x, verts_y, verts_z]

  domain, geom = mesh.rectilinear(verts)
  vec_basis    = domain.splinefunc( degree=degree ).vector( ndims )

  # define wedge problem
  rho = function.select(
      [function.greater(geom[-1]+0.4*Lz+geom[0]/6, 0), function.greater(geom[-1]+0.8*Lz-geom[0]/3, 0)],
      [rho0, rho1], rho2)
  cp  = function.select(
      [function.greater(geom[-1]+0.4*Lz+geom[0]/6, 0), function.greater(geom[-1]+0.8*Lz-geom[0]/3, 0)],
      [cp0, cp1], cp2)
  cs  = function.select(
      [function.greater(geom[-1]+0.4*Lz+geom[0]/6, 0), function.greater(geom[-1]+0.8*Lz-geom[0]/3, 0)],
      [cs0, cs1], cs2)
  mu   = cs**2 * rho
  lam  = rho * (cp**2 - 2.0*cs**2) 
    
  # problem summary
  print( '----     WEDGE PROBLEM     ----' )  
  print( 'problem size   : ' + str(nx-1+degree)+' x '+str(ny-1+degree)+' x '+str(nz-1+degree) )
  print( '# dofs         : ' + str(len(vec_basis)) )
  print( 'max. frequency : ' + str( min(cs0,cs1,cs2,cp0,cp1,cp2)/(ppw*max(dx,dy,dz)) ) )
  print( '-------------------------------\n' )

  # Create discretization matrices using nutils
  K, C, M, rhs = elast_mat(rho, cp, cs, lam, mu, ndims, nx, ny, nz, vec_basis, domain, geom, block)

  # plot_opt_splitting2(damping, 2*pi*1.0, 2*pi*9.0)


  t0 = time.time()
  if solver_flag==0:
      print('Use poly_msgmres of degree '+str(dg_pp))
      sol, it = poly_driver(K.toscipy().tocsc(), C.toscipy().tocsc(), M.toscipy().tocsc(), rhs, freq, tau, damping, tol, maxit, dg_pp, plot_ritz=plot_ritz)
      
      #if path.exists(out_file):
          #remove(out_file)
      with open(out_file, "a") as myfile:
          print(freq)
          myfile.write(''+str(round(min(freq),1))+' '+str(round(max(freq),1))+' '+str(it)+'\n')

  elif solver_flag==1:           
      print('Use megmres')
      sol, it = me_driver(K.toscipy().tocsc(), C.toscipy().tocsc(), M.toscipy().tocsc(), rhs, freq, tau, damping, tol, maxit,  iLU=iLU, fill_factor=fill_factor, rot=rot, plot_ritz=plot_ritz)
            
  else:
      print('Use pythons built-in solver...')
      sol = np.zeros((Nom, len(vec_basis)), dtype=complex)
      it = -1
      for k in range(0,Nom):
          om = 2.0*np.pi*freq[k]*(1.0-1j*damping)
          matrix = K + 1j*om*C - om**2*M
          A = matrix.toscipy().tocsc()
          if ndims==2:
              t0_lu = time.time()
              lu = spla.splu(A)
              print('LU decomposition:'+str(time.time()-t0_lu))
              t0_solve = time.time()
              sol[k,:] = lu.solve(rhs)
              print('solve:'+str(time.time()-t0_solve))
          else:
              print('Use ILU+GMRES')
              class gmres_counter(object):
                  def __init__(self, disp=True):
                      self._disp = disp
                      self.resvec=[]
                      self.niter = 0
                  def __call__(self, rk=None):
                      self.niter += 1
                      self.resvec.append(rk)
                      if self._disp:
                          print('iter %3i\trk = %s' % (self.niter, str(rk)))            
              t0_lu = time.time()
              invA = spla.spilu( A, fill_factor=10.0)
              invA_x = lambda x: invA.solve(x)
              ilu = spla.LinearOperator(A.shape, invA_x)
              print('ilu setup:'+str(time.time()-t0_lu))
              t0_solve = time.time()
              counter = gmres_counter(disp=True)
              sol[k,:], info = spla.gmres(A, rhs, tol=1e-16, restart=200, maxiter=20, M=ilu, callback=counter)
              it = info
              print('GMRES time:'+str(time.time()-t0_solve))
              print('GMRES info:'+str(counter.niter)+' -- '+str(counter.resvec[-1]))
              
  te = time.time()
  print('No iterations: '+str(it)+'     CPU time: '+str(te-t0))        
  test_orig_problem(K.toscipy(), C.toscipy(), M.toscipy(), rhs, freq, damping, sol.T, solver_flag) 
  
  if plots:
      if(ndims ==2):
          makeplots( domain, geom, Lx, Lz, rho, 'rho', 'Density' )
          #makeplots( domain, geom, Lx, Lz, rho, 'rho', 'Density', imgtype='eps' )
          #makeplots( domain, geom, Lx, Lz, cp, 'cp', 'c_p [m/s]' )
          #makeplots( domain, geom, Lx, Lz, cs, 'cs', 'c_s [m/s]' )
          
          for k in range(0,Nom):
              disp     = vec_basis.dot( sol[k,:] )  # FEM summation
              #disp_x   = disp[0].real               # Plot Re(u_x)  
              disp_z   = disp[-1].real              # Plot Re(u_z)
              #makeplots( domain, geom, Lx, Lz, disp_x, 'disp_x'+str(k), 'u_x at {} Hz'.format(freq[k]), lineOn=True )
              makeplots( domain, geom, Lx, Lz, disp_z, 'disp_z'+str(k), 'u_z at {} Hz'.format(freq[k]), lineOn=True )
              #makeplots( domain, geom, Lx, Lz, disp_x, 'disp_x'+str(k), 'u_x at {} Hz'.format(freq[k]), lineOn=True, imgtype='eps' )
              #makeplots( domain, geom, Lx, Lz, disp_z, 'disp_z'+str(k), 'u_z at {} Hz'.format(freq[k]), lineOn=True, imgtype='eps' )
          makevtk(domain, geom, rho, lam, mu, cp, cs, sol, freq, vec_basis, 'wedge2d')

      elif(ndims==3):
          makevtk(domain, geom, rho, lam, mu, cp, cs, sol.T, freq, vec_basis, 'wedge3d')

util.run( main )
