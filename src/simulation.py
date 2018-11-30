from grid_1d import *
import matplotlib.pyplot as plt
from params import sod_params

def minmod( a , b ):
    if abs(a) < abs(b) and a * b > 0.0:
        return a
    elif abs(b) < abs(a) and a * b > 0.0:
        return b
    else:
        return 0.0

def maxmod( a , b ):
    if abs(a) > abs(b) and a * b > 0.0:
        return a
    elif abs(b) < abs(a) and a * b > 0.0:
        return b
    else:
        return 0.0

class Simulation( object ):
    def __init__( self , params , grid ):
        self.grid = grid
        self.t = 0.0

    def set_ICs( self , type='tophat' ):
        if type == 'tophat':
            self.grid.U[ : , np.logical_and( self.grid.xs >= 0.333, \
                    self.grid.xs <= 0.666 ) ] = 1.0

        elif type == "sine":
            # initialize all state variables to 1.0
            self.grid.U[ : , : ] = 1.0

            # indices are the middle third of grid
            index = np.logical_and( self.grid.xs >= 0.333, \
                    self.grid.xs <= 0.666 )

            self.grid.U[ : , index ] += \
                    0.5*np.sin(2.0*np.pi*(self.grid.xs[index]-0.333)/0.333)

        elif type == "rarefaction":
            self.grid.U[:] = 1.0
            self.grid.U[ self.grid.xs > 0.5 ]  = 2.0

    # TODO: this needs to depend on eigenvalues; alphas, not just the velocities
    # but I will only put velocities here to get the structure correct
    def compute_timestep( self , CFL ):
        return( CFL * self.grid.dx / \
                max( abs( self.grid.U[ 1 , self.grid.ilo :\
                self.grid.ihi + 1 ] ) ) )

    def cons_to_prim( self , U ):
        q = self.g.get_scratch_array()
        gamma = self.params['gamma']

        WRHO = self.g.QRHO
        URHO = self.g.URHO
        WV = self.g.WV
        UMX = self.g.URHOV
        WP = self.g.WP
        UENER = self.g.UENER

        q[ WRHO , : ] = U[ URHO , : ]
        q[ WV , : ]   = U[ UMX , : ] / U[ URHO , : ]
        q[ WP , : ]   = (U[ UENER , : ] - 0.5 * q[ WRHO , : ] *\
                        q[ WV , : ]**2) * (gamma - 1.0 )

        return q

    def states( self , dt ):
        """ reconstruct left and right interface states
            implemented:
            minmod

            not implemented:
            godunov
            centered
            MC
            superbee
        """
        g = self.grid

        # compute piecewise linear slopes -- 2nd order MC limiter
        # pick a range of cells that includes 1 ghost cell on either side
        ie = g.ihi + 1

        U = g.U

        dc = g.get_scratch_array()
        dl = g.get_scratch_array()
        dr = g.get_scratch_array()

        dc[ : , g.ilo-1:ie+1] = 0.5* ( U[ : , g.ilo : ie+2 ] - U[ : , g.ilo-2 : ie ] )
        dl[ : , g.ilo-1:ie+1] = U[ : , g.ilo-1+1 : ie+2 ] - U[ : , g.ilo-1 : ie+1 ]
        dr[ : , g.ilo-1:ie+1] = U[ : , g.ilo-1 : ie+1 ] - U[ : , g.ilo-2 : ie ]


        # minmod
        d1 = 2.0 * np.where( np.fabs( dl ) < np.fabs( dr ) , dl , dr )
        d2 = np.where( np.fabs( dc ) < np.fabs( d1 ) , dc , d1 )
        ldeltau = np.where( dl * dr > 0.0 , d2 , 0.0 )

        # now interface states. there is one more interface than zones
        ul = g.get_scratch_array()
        ur = g.get_scratch_array()

        ur[ : , g.ilo-1 : ie + 2 ] = U[ : , g.ilo-1 : ie + 2 ] - \
                0.5 * ( 1.0 + U[ : , g.ilo-1 : ie+2 ] * dt / self.grid.dx ) * \
                ldeltau[ : , g.ilo-1 : ie+2 ]

        ul[ : , g.ilo : ie+2 ] = U[ : , g.ilo-1 : ie+1 ] +\
                0.5 * ( 1.0 - U[ : , g.ilo-1 : ie+1 ] * dt / self.grid.dx ) * \
                ldeltau[ : , g.ilo-1 : ie+1 ]

        return ul , ur

    # TODO: Need to write 3 separate flux functions in here
    # TODO: Need to compute HLL flux in here
    def riemann( self , ul , ur ):
        """
        solve the riemann problem given the left and right states
        returns flux at interfaces, given left and right states
        implemented:
        upwinding
        HLL

        """

        # TODO: need to compute the HLL flux here
        S = 0.5 * ( ul + ur )
        ushock = np.where( S>0.0 , ul , ur )
        ushock = np.where( S==0.0 , 0.0 , ushock )

        # rarefaction solution
        urare = np.where( ur <= 0.0 , ur , 0.0 )
        urare = np.where( ul >= 0.0 , ul , urare )

        us = np.where( ul > ur , ushock , urare )

        return( 0.5 * us * us )

    # TODO: try to use this inside riemann
    def get_HLL_Flux( self , U , ap , am ):
        """ compute update using HLLE flux  """
        FL = self.F[:,:-1]
        FR = self.F[:,1:]
        UL = U[:,:-1]
        UR = U[:,1:]
        FHLL = ( ap * FL + am * FR - ap * am * ( UR - UL )) / (ap+am)
        LU = np.zeros((3,self.Nx))
        LU[:,1:-1] = -(FHLL[:,1:]-FHLL[:,:-1]) / self.dx
        return LU

    # TODO: implement rk4_substep update as well
    def update( self , dt , flux ):
        """ perform the conservative update """

        g = self.grid

        unew = g.get_scratch_array()

        unew[ : , g.ilo : g.ihi+1 ] = \
                g.U[ : , g.ilo : g.ihi+1 ] + \
                dt / g.dx * \
                ( flux[ : , g.ilo : g.ihi+1 ] - flux[ : , g.ilo+1 : g.ihi+2 ] )

        return( unew )

    def evolve( self , C , tmax ):
        # set time to zero
        self.t = 0.0

        # alias for access to grid
        g = self.grid

        while( self.t < tmax ):
            # fill boundary conditions
            g.fill_BCs()

            # get the timestep
            dt = self.compute_timestep( C )

            if ( self.t + dt > tmax ):
                dt = tmax - self.t

            # compute interface states
            ul , ur = self.states( dt )

            # solve riemann problem
            flux = self.riemann( ul , ur )

            # do conservative update
            unew = self.update( dt , flux )
            self.grid.U[:] = unew[:]

            # TODO: update primitives after conservative update

            self.t += dt

if __name__ == "__main__":
    xmin = 0.0
    xmax = 1.0
    nx = 256
    ng = 2
    g = Grid1d_Euler( nx , ng , bc='outflow' )
    tmax = 0.1
    CFL = 0.8
    physical = slice(g.ilo,g.ihi+1)
    fig,ax = plt.subplots(3,1,sharex=True)
    s = Simulation( sod_params , g )

    initial_condition = "sine"
    s.set_ICs( initial_condition )

    uinit = s.grid.U.copy()

    s.evolve( CFL , tmax )

    g = s.grid
    for var in range(g.NVAR):
        ax[var].plot( g.xs[physical],\
                g.U[var,physical], color='k',\
                label='t={}'.format(tmax))
        ax[var].plot( g.xs[physical] ,uinit[var,physical],\
                ls=":",color="red",zorder=-1,\
                label='initial configuration')
        ax[var].legend(loc='best')
        ax[var].set_ylabel(g.cons_vars[var])

    plt.xlabel("$x$")
    # plt.ylabel("$u$")
    plt.savefig("fv-burger-{}.pdf".format( initial_condition ) )
    plt.suptitle("Solution for {} wave".format(initial_condition))
    plt.show()
