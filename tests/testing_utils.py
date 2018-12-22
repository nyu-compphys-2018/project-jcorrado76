import matplotlib.pyplot as plt
import sys
# need both of these because if you call from makefile or from local tests
# directory, need to make sure can see eulersolver
sys.path.insert(0,'../../src')
sys.path.insert(0,'../src')
from EulerSolver import EulerSolver
from utils import f

plots_dir = "../plots/"
plot_extension = ".jpg"
dpi = 400
tfinal = 0.4
CFL = 0.4
N = 500

def run_sod_test( parameters=None ):
    try:
        order = sys.argv[1]
    except:
        order = 'low'
    if order == 'low':
        sod_plots_dir = plots_dir + "low_order_sod_shock_tube_plots/"
    else:
        sod_plots_dir = plots_dir + "high_order_sod_shock_tube_plots/"
    if parameters is None:
        sys.exit()
    if order == 'low':
        title = "LowOrderSpaceLowOrderTime"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL,\
                time_order=1,spatial_order=1)
    else:
        title = "ThirdOrderTimeSecondOrderSpace"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL, time_order=3,spatial_order=2 )

    e.IC_manager.setSod( params = parameters )
    title += "{}".format( parameters['title'] )

    e.evolve( tfinal )
    axes = e.plot(title=title)

    save_string = sod_plots_dir + title + plot_extension
    plt.savefig( save_string , bbox_inches='tight',format='jpg',dpi=dpi)

def run_isentropic_wave_test(parameters=None):
    try:
        order = sys.argv[1]
    except:
        order = 'low'
    bc='periodic'
    # bc='outflow'
    if order == 'low':
        title = "LowOrderSpaceLowOrderTime"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL,\
                time_order=1,spatial_order=1, bc=bc)
    else:
        title = "ThirdOrderTimeSecondOrderSpace"
        e = EulerSolver( Nx=N , a=0.0 , b=1.0 , cfl=CFL,
                time_order=3,spatial_order=2 , bc=bc)
    x0=0.5; sigma=0.4
    e.IC_manager.setIsentropicWave(parameters , f,x0,sigma)
    title += "{}".format( "Isentropic Wave")

    e.evolve( 0.8 )
    axes = e.plot(title=title)

    isentropic_plots_dir = plots_dir + "isentropic_wave_plots/"
    save_string = isentropic_plots_dir + title + plot_extension
    plt.savefig( save_string , bbox_inches='tight',format='jpg',dpi=dpi)
