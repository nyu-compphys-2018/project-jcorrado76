import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../src')
from EulerSolver import EulerSolver
from utils import f

plots_dir = "../plots/"
plot_extension = ".jpg"
dpi = 400
tfinal = 0.4
CFL = 0.4
N = 500

def run_sod_test( parameters=None , order='low' ):
    sod_plots_dir = plots_dir + "sod_shock_tube_plots/"
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

def run_isentropic_wave_test(order='low'):
    bc='periodic'
    if order == 'low':
        title = "Low Order in Space Low order in Time"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL,\
                time_order=1,spatial_order=1, bc=bc)
    else:
        title = "High Order"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL,
                time_order=3,spatial_order=2 , bc=bc)
    rho0 = 1.0; p0 = 0.6; alpha = 0.2; x0=0.5; sigma=0.4
    e.setIsentropicWave(rho0,p0,alpha,f,x0,sigma)
    title += " {}".format( "Isentropic Wave")

    e.evolve( tfinal )
    axes = e.plot(title=title)

    save_string = plots_dir + title + plot_extension
    plt.savefig( save_string , bbox_inches='tight',format='jpg',dpi=dpi)
