import numpy as np
import matplotlib.pyplot as plt
from utils import *
from EulerSolver import EulerSolver
from sod_shock_tube_parameters import *

def main():
    order = 'low'
    # order = 'high'
    tfinal = 0.4
    CFL = 0.4
    N = 500
    if order == 'low':
        title = "Low Order in Space Low order in Time"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL,\
                time_order=3,spatial_order=1)
    else:
        title = "High Order"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=3,spatial_order=2 )

    # params = nonrelativistic_parameters
    params = problem_2_relativistic_parameters

    e.setSod( params=params )
    title += " {} Sod".format( params['title'] )

    # SMOOTH WAVE 
    # e.setSmoothWave()
    # rho0 = 1.0; p0 = 0.6; alpha = 0.2; x0=0.5; sigma=0.4
    # e.setIsentropicWave(rho0,p0,alpha,f,x0,sigma)

    e.evolve( tfinal )
    axes = e.plot(title=title)

    plt.show()

main()
