import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../src')
from EulerSolver import EulerSolver

plot_extension = ".jpg"
dpi = 400

def run_sod_test( parameters=None , order='low' ):
    tfinal = 0.4
    CFL = 0.4
    N = 500
    sod_plots_dir = "../plots/sod_shock_tube_plots/"
    if parameters is None:
        sys.exit()
    if order == 'low':
        title = "Low Order in Space Low order in Time"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL,\
                time_order=1,spatial_order=1)
    else:
        title = "High Order"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=3,spatial_order=2 )

    e.setSod( params = parameters )
    title += " {}".format( parameters['title'] )

    e.evolve( tfinal )
    axes = e.plot(title=title)

    save_string = sod_plots_dir + title + plot_extension
    plt.savefig( save_string , bbox_inches='tight',format='jpg',dpi=dpi)
