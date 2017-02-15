'''

Use the ``bokeh serve`` command to run the example by executing:

    bokeh serve opt_tau_bokeh.py

at your command prompt. Then navigate to the URL

    http://localhost:5006/opt_tau_online

in your browser.

'''

from bokeh.io import curdoc, vform
from bokeh.layouts import row, column, widgetbox
from bokeh.models import ColumnDataSource, Spacer, Label
from bokeh.models.widgets import Slider, TextInput, Button, RangeSlider
from bokeh.plotting import figure
from bokeh.models.widgets import Div

import numpy as np
from numpy import pi, sin, cos, sqrt
from math import atan2, atan

NOP = 50

def opt_tau_anal(e,w,W):    
    r  = sqrt(w*W*(1.0+e**2))
    th = atan(-sqrt( (e**2*(W+w)**2+(W-w)**2) /(4.0*w*W) ))
    return r*cos(th) + 1j*(r*sin(th))

def J_opt(e, w, W):
    tau = opt_tau_anal(e,w,W)
    r     = 0.5*np.sqrt(1.0 + (tau.real/tau.imag)**2)
    c1_im = tau.real/(2.0*tau.imag) - ((tau.imag+e*tau.real)*w)/((w-tau.real)**2+(e*w+tau.imag)**2)
    #cN_im = tau.real/(2.0*tau.imag) - ((tau.imag+e*tau.real)*W)/((W-tau.real)**2+(e*W+tau.imag)**2)
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
    C   = 0.0 + 1j*( (e*abs(tau)**2)/(2.0*tau.imag*(tau.imag+e*tau.real)) )
    R   = sqrt( abs(tau)**2*(e**2+1.0)/(4.0*(tau.imag+e*tau.real)**2) )
    ck  = np.zeros((len(om),), dtype=complex)
    for k in range(0,len(om)):
        ck[k] = -np.conj(tau)/(tau-np.conj(tau)) - eta[k]
    r  = abs(tau/(tau-np.conj(tau)))
    return C, R, ck, r


# Set up data
Nom  = 7
fmin = 1.0
fmax = 9.0
eps  = 0.7
freq = np.linspace(fmin,fmax,Nom)
om   = 2.0*np.pi*freq*(1.0-1j*eps)
tau  = opt_tau_anal(eps,om[0].real,om[-1].real) 
Jopt = J_opt(eps, min(om.real), max(om.real))
th   = np.linspace(0.0, 2.0*pi, NOP)

cc  = list(['green','blue','cyan','magenta','yellow'])
col = list(['red'])
j = -1
for k in range(1,Nom-1):
    j=j+1
    if (j>4):
        j=0
    col.append(cc[j])
col.append('red')

C, R, c, r = calc_circles(freq, tau, eps)
X = R*np.cos(th)+C.real
Y = R*np.sin(th)+C.imag
bigC     = ColumnDataSource(data=dict(x=X, y=Y))
bigC_cen = ColumnDataSource(data=dict(x=[C.real], y=[C.imag]))

circ = []
circ_cen = []
for k in range(0,len(om)):
    x = r*np.cos(th)+c[k].real
    y = r*np.sin(th)+c[k].imag
    circ.append(ColumnDataSource(data=dict(x=x, y=y)))
    circ_cen.append(ColumnDataSource(data=dict(x=[c[k].real], y=[c[k].imag])))

# Set up plot
plot = figure(plot_height=640, plot_width=640,
              tools="pan,box_zoom,reset,save,wheel_zoom", active_drag=None,
              x_range=[-2, 2], y_range=[-1, 3])

plot.line(x=[-10,10],y=[0,0], line_width=1, line_color='black')
plot.line(x=[0,0],y=[10,-10], line_width=1, line_color='black')

mytext  = Label(x=-1.8, y=2.7, text='J = '+str(round(Jopt,4)), text_font_size='22pt')
mytext2 = Label(x=-1.8, y=2.5, text='J* = '+str(round(Jopt,4)), text_font_size='16pt', text_color="grey")
plot.add_layout(mytext)
plot.add_layout(mytext2)

plot.line('x', 'y', source=bigC, line_width=3, line_alpha=0.6, line_color='black')
plot.circle('x', 'y', source=bigC_cen, size=10, color="black", alpha=0.6)

for k in range(0,len(om)):
    plot.line('x', 'y', source=circ[k], line_width=3, line_alpha=0.6, line_color=col[k])
    plot.circle('x', 'y', source=circ_cen[k], size=8, color=col[k], alpha=0.6)


# Set up widgets
tau_re_s = Slider(title="Real part", value=tau.real/(2*pi*fmax), start=0.0, end=1.0, step=0.001)
tau_im_s = Slider(title="Imag part", value=tau.imag/(2*pi*fmax), start=-1.0, end=0.0, step=0.001)

eps_s    = Slider(title="viscous damping", value=eps, start=0.0, end=1.0, step=0.01)
#Nom_text = TextInput(title="n_f = ", value=str(Nom))
fmin_s   = Slider(title="f_min [Hz]", value=fmin, start=0.0, end=15.0, step=0.1)
fmax_s   = Slider(title="f_max [Hz]", value=fmax, start=0.0, end=15.0, step=0.1)
frange_s = RangeSlider(start=0.0, end=15.0, range=(fmin,fmax), step=0.1, title="freq. range [Hz]")
#frange_s = DateRangeSlider(bounds=(0.0,15.0), value=(fmin,fmax), title="freq range [Hz]")


reset = Button(label="Reset to optimal tau")

# Set up callbacks
def update_title(attrname, old, new):
    plot.title.Nom_text = Nom_text.value
Nom_text.on_change('value', update_title)

def update_data1(attrname, old, new):

    # Get the current slider values
    #fmin = fmin_s.value
    #fmax = fmax_s.value
    fmin = frange_s.range[0]
    fmax = frange_s.range[1]
    
    Nom  = int(Nom_text.value)
    freq = np.linspace(fmin,fmax,Nom)
    
    C, R, c, r = calc_circles(freq, tau_re_s.value*(2*pi*fmax)+1j*tau_im_s.value*(2*pi*fmax), eps_s.value)
    X   = R*np.cos(th)+C.real
    Y   = R*np.sin(th)+C.imag

    # Generate the new curveeps_s
    bigC.data     = dict(x=X, y=Y)
    bigC_cen.data = dict(x=[C.real], y=[C.imag])
    for k in range(0,Nom):
        x = r*np.cos(th)+c[k].real
        y = r*np.sin(th)+c[k].imag
        circ[k].data     = dict(x=x, y=y)
        circ_cen[k].data = dict(x=[c[k].real], y=[c[k].imag])
        
    Jnew = J(2.0*pi*freq*(1.0-1j*eps_s.value), tau_re_s.value*(2*pi*fmax)+1j*tau_im_s.value*(2*pi*fmax))
    mytext.text = 'J = '+str(round(Jnew,4))
        

def update_data2(attrname, old, new):

    # Get the current slider values
    #fmin = fmin_s.value
    #fmax = fmax_s.value
    fmin = frange_s.range[0]
    fmax = frange_s.range[1]
    Nom  = int(Nom_text.value)
    freq = np.linspace(fmin,fmax,Nom)
    
    opt_tau = opt_tau_anal(eps_s.value,2*pi*fmin,2*pi*fmax)
    tau_re_s.value = opt_tau.real/(2*pi*fmax)
    tau_im_s.value = opt_tau.imag/(2*pi*fmax)
    
    C, R, c, r = calc_circles(freq, tau_re_s.value*(2*pi*fmax)+1j*tau_im_s.value*(2*pi*fmax), eps_s.value)
    X   = R*np.cos(th)+C.real
    Y   = R*np.sin(th)+C.imag

    # Generate the new curveeps_s
    bigC.data     = dict(x=X, y=Y)
    bigC_cen.data = dict(x=[C.real], y=[C.imag])
    for k in range(0,Nom):
        x = r*np.cos(th)+c[k].real
        y = r*np.sin(th)+c[k].imag
        circ[k].data = dict(x=x, y=y)
        circ_cen[k].data = dict(x=[c[k].real], y=[c[k].imag])
     
    Jnew = J(2.0*pi*freq*(1.0-1j*eps_s.value), tau_re_s.value*(2*pi*fmax)+1j*tau_im_s.value*(2*pi*fmax))
    Jopt = Jopt = J_opt(eps_s.value, 2*pi*fmin, 2*pi*fmax)
    mytext.text = 'J = '+str(round(Jnew,4))
    mytext2.text = 'J* = '+str(round(Jopt,4))
     
def click():
    #fmin = fmin_s.value
    #fmax = fmax_s.value
    fmin = frange_s.range[0]
    fmax = frange_s.range[1]
    opt_tau = opt_tau_anal(eps_s.value,2*pi*fmin,2*pi*fmax)
    tau_re_s.value = opt_tau.real/(2*pi*fmax)
    tau_im_s.value = opt_tau.imag/(2*pi*fmax)
    

for w in [tau_re_s, tau_im_s]:
    w.on_change('value', update_data1)
eps_s.on_change('value', update_data2)
frange_s.on_change('range', update_data2)    
reset.on_click(click)
    

# Set up layouts and add to document
inputs_1 = widgetbox(tau_re_s, tau_im_s)
#inputs_2 = widgetbox(Nom_text, fmin_s, fmax_s, eps_s)
inputs_2 = widgetbox(frange_s, eps_s)
inputs_2 = widgetbox(frange_s, eps_s)
inputs_3 = vform(reset)
spacer_1 = Spacer(width=20, height=40)
spacer_2 = Spacer(width=20, height=40)
spacer_3 = Spacer(width=20, height=130)


div1 = Div(text="""<font size="4"><b>Seed parameter tau</b></font> <br> (relative w.r.t. f_max)""",
width=300, height=30)
text1 = widgetbox(div1)

div2 = Div(text="""<font size="4"><b>Frequency range & damping</b></font>""",
width=300, height=20)
text2 = widgetbox(div2)

div3 = Div(text="""A link to our preprint that describes the present visualization can be found <a href="http://www.ewi.tudelft.nl/en/the-faculty/departments/applied-mathematics/reports/">[here]</a>.""",
width=300, height=20)
text3 = widgetbox(div3)


curdoc().add_root(row(column(text1, inputs_1, spacer_1, text2, inputs_2, spacer_2, text3, spacer_3, inputs_3), plot, width=1200))
curdoc().title = "Optimal tau applet"
