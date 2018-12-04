import numpy as np
import matplotlib.pyplot as plt
from utils import *
from EulerSolver import EulerSolver
from eulerExact import riemann

def plot_sod_shock_tube_convergence( e , order='low'):
    """
    this function plots the convergence rate of then
    scheme
    """
    Ns = [ 32 , 64 , 128 , 256 , 512 ]
    t = 0.1
    b = e.b
    a = e.a
    x0 = ( b + a ) / 2
    rhol = 1.0; vl = 0.0; pl = 1.0
    rhor = 0.1; vr = 0.0; pr = 0.125
    gamma=e.gamma

    rerr = []
    verr = []
    perr = []

    L1 = np.zeros( len(Ns) )
    for i in range(len(Ns)):
        deltaX = (b-a)/Ns[i]
        if order == 'low':
            e = EulerSolver(Nx=Ns[i],spatial_order=1,time_order=1)
        else:
            e = EulerSolver(Nx=Ns[i], spatial_order=2,time_order=3)
        e.setSod(left_states=[rhol,vl,pl],right_states=[rhor,vr,pr])
        e.evolve(t)
        xexact , rexact , vexact , pexact = riemann( a , b , x0 , Ns[i] , t ,rhol ,vl , pl , rhor , vr , pr , gamma )
        rerr.append( compute_l1_error( e.W[0,e.physical] , rexact , deltaX ))
        verr.append( compute_l1_error( e.W[1,e.physical] , vexact , deltaX ))
        perr.append( compute_l1_error( e.W[2,e.physical] ,   pexact , deltaX ))

    print("L1 Errors on Density",rerr)
    print("L1 Errors on Velocity",verr)
    print("L1 Errors on Pressure",perr)

    logns = np.log( Ns )
    logrs = np.log( rerr )
    logvs = np.log( verr )
    logps = np.log( perr )

    rb , rm = fit_line( logns , logrs )
    vb , vm = fit_line( logns , logvs )
    pb , pm = fit_line( logns , logps )

    print("Density slope:" , rm)
    print("Velocity slope:" , vm)
    print("Pressure slope:" , pm)

    fig = plt.figure()

    plt.plot( logns, logrs , color='red' , marker='o',label="L1 Error on Density")
    plt.plot( logns, logvs , color='green' , marker='o',label="L1 Error on Velocity")
    plt.plot( logns, logps , color='blue' , marker='o',label="L1 Error on Pressure")

    plt.plot( logns , line(logns , rb , rm ) , linestyle='dashed', color='red' , label='L1 Error Fit on Density')
    plt.plot( logns , line(logns , vb , vm ) , linestyle='dashed',color='green' , label='L1 Error Fit on Velocity')
    plt.plot( logns , line(logns , pb , pm ) , linestyle='dashed',color='blue' , label='L1 Error Fit on Pressure')

    plt.title("Convergence Plot for Sod Shock Tube")
    plt.xlabel("$log(N)$")
    plt.ylabel("$log(L_1)$")
    plt.legend()


def main():
    order = 'low'
    CFL = 0.5
    N = 400

    if order == 'low':
        title = "Low Order"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=1,spatial_order=1 )
    else:
        title = "High Order"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=3,spatial_order=2 )
    plot_sod_shock_tube_convergence( e , order=order )
    plt.show()

main()
