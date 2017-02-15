from numpy import pi, sin, cos, sqrt
from math import atan2, atan
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

axis_color = 'lightgoldenrodyellow'
Nom_min = 3
Nom_max = 10

def opt_tau_anal(e,w,W):    
    r  = sqrt(w*W*(1.0+e**2))
    th = atan(-sqrt( (e**2*(W+w)**2+(W-w)**2) /(4.0*w*W) ))
    return r*cos(th) + 1j*(r*sin(th))

def J_opt(e, w, W):
    tau = opt_tau_anal(e,w,W)
    r     = 0.5*np.sqrt(1.0 + (tau.real/tau.imag)**2)
    c1_im = tau.real/(2.0*tau.imag) - ((tau.imag+e*tau.real)*w)/((w-tau.real)**2+(e*w+tau.imag)**2)
    cN_im = tau.real/(2.0*tau.imag) - ((tau.imag+e*tau.real)*W)/((W-tau.real)**2+(e*W+tau.imag)**2)
    R     = np.sqrt(tau.real**2+tau.imag**2)*np.sqrt((e**2+1.0))/(2.0*abs(tau.real*e+tau.imag))
    C_im  = e*(tau.real**2+tau.imag**2)/(2.0*tau.imag*(tau.real*e+tau.imag))
    return np.sqrt(r**2/(R**2-C_im**2+2.0*C_im*c1_im))

def J(om, tau):
    eta    = om/(om-tau)
    Jval = np.zeros((len(om),))
    r = abs((tau - 0)/(tau - np.conj(tau)))
    for k in range(len(om)):
        ck = ((0.0 - np.conj(tau))/(tau - np.conj(tau)) - eta[k])
        Jval[k] = r/abs(ck)
    return np.max(Jval)
    
  
def calc_circles(freq, tau, e):
    om  = 2.0*pi*freq*(1.0-1j*e)
    eta = om/(om-tau)
    #e   = -om[0].imag/om[0].real
    C   = 0.0 + 1j*( (e*abs(tau)**2)/(2.0*tau.imag*(tau.imag+e*tau.real)) )
    R   = sqrt( abs(tau)**2*(e**2+1.0)/(4.0*(tau.imag+e*tau.real)**2) )
    ck  = np.zeros((len(om),), dtype=complex)
    for k in range(0,len(om)):
        ck[k] = -np.conj(tau)/(tau-np.conj(tau)) - eta[k]
    r  = abs(tau/(tau-np.conj(tau)))
    return C, R, ck, r

def draw_circles(tau_re, tau_im, eps, freq):
    om  = 2.0*pi*freq*(1.0-1j*eps)
    NOP = 1000
    th  = np.linspace(0.0, 2.0*pi, NOP)
    
    C, R, c, r = calc_circles(freq, tau_re+1j*tau_im, eps)
    X   = R*np.cos(th)+C.real
    Y   = R*np.sin(th)+C.imag

    #fig.subplots_adjust(left=0.25, bottom=0.25)

    lC1, = ax.plot(X, Y, 'k')
    lC2, = ax.plot(C.real, C.imag, 'kx', markersize=10)
    
    lc1 = []
    lc2 = []
    for k in range(0,len(om)):
        x = r*np.cos(th)+c[k].real
        y = r*np.sin(th)+c[k].imag
        l1, = ax.plot(x, y, col[k]+'--')
        l2, = ax.plot(c[k].real, c[k].imag, color=col[k], marker='x', markersize=10)
        lc1.append(l1)
        lc2.append(l2)
        
    ax.axhline(linewidth=0.5, color='k')
    ax.axvline(linewidth=0.5, color='k')
    
    #Jopt = J_opt(eps, min(om.real), max(om.real))
    Jopt = J(om, tau_re+1j*tau_im)
    txt  = ax.text(0.7*ax.get_xlim()[1], 0.88*ax.get_ylim()[1], r'$\mathcal{J} = $'+str(round(Jopt,4)), fontsize=22)
    txt2  = ax.text(0.7*ax.get_xlim()[1], 0.8*ax.get_ylim()[1], r'$\mathcal{J}^\ast = $'+str(round(Jopt,4)), fontsize=15)
    
    return lC1, lC2, lc1, lc2, txt, txt2
    
def upd_circles(tau_re, tau_im, eps, freq, Jopt):
    om  = 2.0*pi*freq*(1.0-1j*eps)
    NOP = 1000
    th  = np.linspace(0.0, 2.0*pi, NOP)
    
    C, R, c, r = calc_circles(freq, tau_re+1j*tau_im, eps)
    X   = R*np.cos(th)+C.real
    Y   = R*np.sin(th)+C.imag

    lC1.set_xdata(X)
    lC1.set_ydata(Y)
    lC2.set_xdata(C.real)
    lC2.set_ydata(C.imag)
    for k in range(0,len(om)):
        x = r*np.cos(th)+c[k].real
        y = r*np.sin(th)+c[k].imag
        lc1[k].set_xdata(x)
        lc1[k].set_ydata(y)
        lc2[k].set_xdata(c[k].real)
        lc2[k].set_ydata(c[k].imag)
        
    #Jnew = J(eps, min(om.real), max(om.real), tau=tau_re+1j*tau_im)
    Jnew = J(om, tau_re+1j*tau_im)
    txt.set_text(r'$\mathcal{J} = $'+str(round(Jnew,4)))
    txt2.set_text(r'$\mathcal{J}^\ast = $'+str(round(J_opt(eps, min(om.real), max(om.real)),4)))
    if (Jopt-Jnew)>1e-6:
        txt.set_color('red')
    else:
        txt.set_color('black')    
    
    
# Default params
Nom  = 8
fmin = 1.0
fmax = 9.0
eps  = 0.7
freq = np.linspace(fmin,fmax,Nom)
om   = 2.0*np.pi*freq*(1.0-1j*eps)
tau  = opt_tau_anal(eps,om[0].real,om[-1].real) 
Jopt = J_opt(eps, min(om.real), max(om.real))

cc  = list('gbcmy')
col = list('r')
j = -1
for k in range(1,Nom-1):
    j=j+1
    if (j>4):
        j=0
    col.append(cc[j])
col.append('r')
    
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
ax.axis('equal')
fig.subplots_adjust(left=0.25, bottom=0.25)

lC1, lC2, lc1, lc2, txt, txt2 = draw_circles(tau.real, tau.imag, eps, freq)


# Add sliders for tweaking the parameters
eps_slider_ax    = fig.add_axes([0.04, 0.5, 0.15, 0.05], axisbg=axis_color)
eps_slider       = Slider(eps_slider_ax, r'$\epsilon$', 0.0, 1.0, valinit=eps)

tau_re_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axis_color)
tau_re_slider    = Slider(tau_re_slider_ax, r'Re($\tau$)', 0.0, 1.0, valinit=tau.real/(2*pi*fmax))
tau_re_slider.vline.color='b'
tau_im_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axis_color)
tau_im_slider    = Slider(tau_im_slider_ax, r'Im($\tau$)', -1.0, 0.0, valinit=tau.imag/(2*pi*fmax))
tau_im_slider.vline = None

fmin_slider_ax   = fig.add_axes([0.04, 0.8, 0.15, 0.05], axisbg=axis_color)
fmin_slider      = Slider(fmin_slider_ax, r'$f_{min}$', 0.5, 20.0, valinit=fmin)
fmax_slider_ax   = fig.add_axes([0.04, 0.7, 0.15, 0.05], axisbg=axis_color)
fmax_slider      = Slider(fmax_slider_ax, r'$f_{max}$', 1.0, 20.0, valinit=fmax)
Nom_slider_ax    = fig.add_axes([0.04, 0.6, 0.15, 0.05], axisbg=axis_color)
Nom_slider       = Slider(Nom_slider_ax, r'$N_f$', Nom_min, Nom_max, valinit=Nom, valfmt='%0.0f')

def sliders_on_changed1(val):
    fmin = fmin_slider.val  
    fmax = fmax_slider.val
    Nom  = Nom_slider.val
    freq = np.linspace(fmin,fmax,Nom)
    
    opt_tau = opt_tau_anal(eps_slider.val,2*pi*fmin,2*pi*fmax)
    tau_re_slider.set_val( opt_tau.real/(2*pi*fmax) )
    tau_im_slider.set_val( opt_tau.imag/(2*pi*fmax) )
    
    upd_circles(tau_re_slider.val*(2*pi*fmax), tau_im_slider.val*(2*pi*fmax), eps_slider.val, freq, Jopt)
    fig.canvas.draw_idle()

def sliders_on_changed2(val):
    fmin = fmin_slider.val  
    fmax = fmax_slider.val
    Nom  = Nom_slider.val
    freq = np.linspace(fmin,fmax,Nom)
    
    upd_circles(tau_re_slider.val*(2*pi*fmax), tau_im_slider.val*(2*pi*fmax), eps_slider.val, freq, Jopt)
    fig.canvas.draw_idle()

eps_slider.on_changed(sliders_on_changed1)
fmin_slider.on_changed(sliders_on_changed1)
fmax_slider.on_changed(sliders_on_changed1)
Nom_slider.on_changed(sliders_on_changed1)

tau_re_slider.on_changed(sliders_on_changed2)
tau_im_slider.on_changed(sliders_on_changed2)

# Add a button for resetting the parameters
reset_button_ax = fig.add_axes([0.06, 0.12, 0.1, 0.04])
reset_button = Button(reset_button_ax, r"Reset $\mathbf{\tau^\ast}$", color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    fmin = fmin_slider.val  
    fmax = fmax_slider.val
    Nom  = Nom_slider.val
    freq = np.linspace(fmin,fmax,Nom)
    opt_tau = opt_tau_anal(eps_slider.val,2*pi*fmin,2*pi*fmax)
    
    tau_re_slider.set_val( opt_tau.real/(2*pi*fmax) )
    tau_im_slider.set_val( opt_tau.imag/(2*pi*fmax) )

reset_button.on_clicked(reset_button_on_clicked)

plt.show()
