import numpy as np
import matplotlib.pyplot as plt
from eulerExact import riemann
from utils import *
from EulerSolver import EulerSolver, plot_convergence


def main():
    order = 'low'
    tfinal = 0.2
    CFL = 0.01
    N = 100
    if order == 'low':
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=1,spatial_order=1 )
    else:
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=2,spatial_order=2 )
    # set initial conditions
    # relativistic sod conditions
    # left_states = [10.0,0.0,13.3]
    # right_states = [1.0,0.0,0.66e-6]
    # left_states = [10.0,0.0,13.3]
    # right_states = [1.0,0.0,0.66e-8]
    # e.setSod(left_states=left_states,right_states=right_states)
    e.setSod()
    # e.setSmoothWave()
    title="Sod Problem"
    # rho0 = 1.0; p0 = 0.6; alpha = 0.2; x0=0.5; sigma=0.4
    # e.setIsentropicWave(rho0,p0,alpha,f,x0,sigma)
    winit = e.W.copy()
    e.evolve( tfinal )
    axes = e.plot(title=title)
    init_labels = ['Initial Density','Initial Velocity','Initial Pressure']
    for i, axis in enumerate(axes):
        axis.plot( e.x , winit[i,:], label=init_labels[i],linestyle='dashed',alpha=0.7)
        axis.legend()

    # if order == 'low':
    #     plot_convergence(order='low')
    # else:
    #     plot_convergence(order='high')
    plt.show()

main()
