import numpy as np
import matplotlib.pyplot as plt
from eulerExact import riemann
from utils import *
from EulerSolver import EulerSolver
from plot_sod_shock_tube_convergence import plot_sod_shock_tube_convergence

def main():
    order = 'low'
    # order = 'high'
    tfinal = 0.4
    CFL = 0.4
    N = 500
    if order == 'low':
        title = "Low Order in Space High order in Time"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL,\
                time_order=3,spatial_order=1)
    else:
        title = "High Order"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=3,spatial_order=2 )

    relativistic_parameters = {"rho_l":10.0 ,"v_l":0.0 , "p_l":13.33,\
                               "rho_r":1.0 ,   "v_r":0.0, "p_r":1e-8 ,\
                               "gamma":5./3. , title:"Relativistic" }
    nonrelativistic_parameters = {"rho_l":1.0 ,"v_l":0.0 , "p_l":1.0,\
                               "rho_r":0.125 ,   "v_r":0.0, "p_r":0.1 ,\
                               "gamma":1.4 , title:"NonRelativistic" }

    params = relativistic_parameters

    e.setSod( parameters=params )
    title += " {} Sod".format( params['title'] )




    # SMOOTH WAVE 
    # e.setSmoothWave()
    # rho0 = 1.0; p0 = 0.6; alpha = 0.2; x0=0.5; sigma=0.4
    # e.setIsentropicWave(rho0,p0,alpha,f,x0,sigma)




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
