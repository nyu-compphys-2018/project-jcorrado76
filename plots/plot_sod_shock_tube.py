import numpy as np
import matplotlib.pyplot as plt
from eulerExact import riemann
from utils import *
from EulerSolver import EulerSolver

def main():
    order = 'low'
    tfinal = 0.1
    CFL = 0.5
    N = 400
    if order == 'low':
        title = "Low Order"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=1,spatial_order=1 )
    else:
        title = "High Order"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=3,spatial_order=2 )
    # set initial conditions
    left_states=[1,0,1]
    right_states=[0.125,0.0,0.1]
    e.setSod(left_states=left_states , right_states=right_states)
    title+=" Sod"
    winit = e.W.copy()
    e.evolve( tfinal )
    axes = e.plot(title=title)
    init_labels = ['Initial Density','Initial Velocity','Initial Pressure']
    exact_labels=['Exact Density','Exact Velocity','Exact Pressure']
    a = e.a
    b = e.b
    x0 = (e.b + e.a)/ 2.
    rhol = left_states[0]; vl = left_states[1]; pl = left_states[2];
    rhor = right_states[0]; vr = right_states[1]; pr = right_states[2];
    gamma = e.gamma
    xexact , rexact , vexact , pexact = riemann( a , b , x0 , e.Nx , tfinal ,rhol ,vl , pl , rhor , vr , pr , gamma )
    exact_primitives = [rexact , vexact , pexact]
    for i, axis in enumerate(axes):
        axis.plot( xexact , exact_primitives[i], label=exact_labels[i],color='red',alpha=0.4)
        # axis.plot( e.x , winit[i,:], label=init_labels[i],linestyle='dashed',alpha=0.7)
        axis.legend()

    plt.show()

main()
