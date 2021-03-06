import numpy as np
import matplotlib.pyplot as plt
from utils import *
from EulerSolver import EulerSolver
import sys
sys.path.insert(0,'../tests')
from parameters import *

def main():
    # order = 'low'
    order = 'high'
    tfinal = 0.4
    CFL = 0.3
    N = 400
    if order == 'low':
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL,\
                time_order=1,spatial_order=1)
    else:
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=3,spatial_order=2 )

    # problem_1_relativistic_parameters
    # problem_2_relativistic_parameters
    # problem_3_relativistic_parameters
    # nonrelativistic_parameters

    # SOD SHOCK TUBE 
    params =problem_1_relativistic_parameters 
    e.IC_manager.setSod( params=params )
    e.title += " {} Sod".format( params['title'] )

    # SMOOTH WAVE 
    # e.IC_manager.setSmoothWave()
    # e.title += " Smooth Wave"

    # ISENTROPIC WAVE 
    # rho0 = 1.0; p0 = 0.6; alpha = 0.2; x0=0.5; sigma=0.4
    # e.setIsentropicWave(rho0,p0,alpha,f,x0,sigma)
    # e.title += " Isentropic Wave"

    e.evolve( tfinal )
    axes = e.plot()

    plt.show()

main()
