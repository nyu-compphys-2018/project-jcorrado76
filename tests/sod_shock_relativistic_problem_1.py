import matplotlib.pyplot as plt
from ..EulerSolver import EulerSolver

problem_1_relativistic_parameters = {"rho_l":10.0 ,"v_l":0.0 , "p_l":13.33,\
                                     "rho_r":1.0 ,   "v_r":0.0, "p_r":1e-8 ,\
                                     "gamma":5./3. , "title":"Relativistic Sod Problem 1" }

order = 'low'
# order = 'high'
tfinal = 0.4
CFL = 0.4
N = 500
if order == 'low':
    title = "Low Order in Space Low order in Time"
    e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL,\
            time_order=1,spatial_order=1)
else:
    title = "High Order"
    e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=3,spatial_order=2 )

e.setSod( params = problem_1_relativistic_parameters )
title += " {}".format( params['title'] )

e.evolve( tfinal )
axes = e.plot(title=title)

plt.show()
