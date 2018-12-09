import numpy as np
import matplotlib.pyplot as plt
from eulerExact import riemann
from utils import *
from EulerSolver import EulerSolver

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

    # SOD SHOCK TUBE PARAMETERS
    problem_1_relativistic_parameters = {"rho_l":10.0 ,"v_l":0.0 , "p_l":13.33,\
                               "rho_r":1.0 ,   "v_r":0.0, "p_r":1e-8 ,\
                               "gamma":5./3. , "title":"Relativistic Problem 1" }
    problem_2_relativistic_parameters = {"rho_l":1.0 ,"v_l":0.0 , "p_l":1000.0,\
                               "rho_r":1.0 ,   "v_r":0.0, "p_r":1e-2 ,\
                               "gamma":5./3. , "title":"Relativistic Problem 2" }
    problem_3_relativistic_parameters = {"rho_l":1.0 ,"v_l":0.9 , "p_l":1.0,\
                               "rho_r":1.0 ,   "v_r":0.0, "p_r":10.0 ,\
                               "gamma":4./3. , "title":"Relativistic Problem 3" }
    nonrelativistic_parameters = {"rho_l":1.0 ,"v_l":0.0 , "p_l":1.0,\
                               "rho_r":0.125 ,   "v_r":0.0, "p_r":0.1 ,\
                               "gamma":1.4 , "title":"NonRelativistic" }

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
