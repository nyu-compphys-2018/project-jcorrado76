from grid_1d import *
import matplotlib.pyplot as plt 


class Simulation( object ):
    def __init__( self , grid ):
        self.grid = grid 
        self.t = 0.0

    def set_ICs( self , type='tophat' ):
        if type == 'tophat':
            self.grid.U[ np.logical_and( self.grid.xs >= 0.333, \
                    self.grid.xs <= 0.666 ) ] = 1.0

        elif type == "sine":
            # initialize all state variables to 1.0
            self.grid.U[:] = 1.0

            # indices are the middle third of grid
            index = np.logical_and( self.grid.xs >= 0.333, \
                    self.grid.xs <= 0.666 ) 

            self.grid.U[ index ] += \
                    0.5*np.sin(2.0*np.pi*(self.grid.xs[index]-0.333)/0.333)

        elif type == "rarefaction":
            self.grid.U[:] = 1.0
            self.grid.U[ self.grid.xs > 0.5 ]  = 2.0 

    def compute_timestep( self , CFL ):
        return( CFL * self.grid.dx / \
                max( abs( self.grid.U[ self.grid.ilo :\
                self.grid.ihi + 1 ] ) ) )

    def states( self , dt ):
        """ compute left and right interface states """
        g = self.grid

        # compute piecewise linear slopes -- 2nd order MC limiter
        # pick a range of cells that includes 1 ghost cell on either side 
        ib = g.ilo - 1
        ie = g.ihi + 1

        U = g.U

        dc = g.get_scratch_array()
        dl = g.get_scratch_array()
        dr = g.get_scratch_array()

        dc[ib:ie+1] = 0.5* (U[ ib+1 : ie+2 ] - U[ ib-1 : ie ] )
        dl[ib:ie+1] = U[ ib+1 : ie+2 ] - U[ ib : ie+1 ]
        dr[ib:ie+1] = U[ ib : ie+1 ] - U[ ib-1 : ie ]


        # minmod 
        d1 = 2.0 * np.where( np.fabs( dl ) < np.fabs( dr ) , dl , dr )
        d2 = np.where( np.fabs( dc ) < np.fabs( d1 ) , dc , d1 )
        ldeltau = np.where( dl * dr > 0.0 , d2 , 0.0 )

        # now interface states. there is one more interface than zones
        ul = g.get_scratch_array()
        ur = g.get_scratch_array()


        ur[ ib : ie + 2 ] = U[ ib : ie + 2 ] - \
                0.5 * ( 1.0 + U[ ib : ie+2 ] * dt / self.grid.dx ) * \
                ldeltau[ ib : ie+2 ]

        ul[ ib+1 : ie+2 ] = U[ ib : ie+1 ] +\
                0.5 * ( 1.0 - U[ ib : ie+1 ] * dt / self.grid.dx ) * \
                ldeltau[ib : ie+1 ]

        return ul , ur 

    def riemann( self , ul , ur ):
        """ 
        returns flux at interfaces, given left and right states 
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

    def get_HLL_Flux( self , U , ap , am ):
        FL = self.F[:,:-1]
        FR = self.F[:,1:]
        UL = U[:,:-1]
        UR = U[:,1:]
        FHLL = ( ap * FL + am * FR - ap * am * ( UR - UL )) / (ap+am)
        LU = np.zeros((3,self.Nx))
        LU[:,1:-1] = -(FHLL[:,1:]-FHLL[:,:-1]) / self.dx
        return LU

    def update( self , dt , flux ):
        """ perform the conservative update """

        g = self.grid

        unew = g.get_scratch_array()

        unew[ g.ilo : g.ihi+1 ] = \
                g.U[ g.ilo : g.ihi+1 ] + \
                dt / g.dx * \
                ( flux[ g.ilo : g.ihi+1 ] - flux[ g.ilo+1 : g.ihi+2 ] )

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

            self.t += dt

if __name__ == "__main__":
    xmin = 0.0
    xmax = 1.0
    nx = 256
    ng = 2
    g = Grid1d( nx , ng , bc='periodic' )
    tmax = (xmax - xmin)/1.0

    C = 0.8

    plt.clf()

    s = Simulation( g )

    for i in range(10):
        tend = (i+1)*0.02*tmax

        s.set_ICs( "sine" )

        uinit = s.grid.U.copy()

        s.evolve( C , tend )

        c = 1.0 - ( 0.1 + i * 0.1 )

        g = s.grid
        plt.plot( g.xs[g.ilo:g.ihi+1],\
                g.U[g.ilo:g.ihi+1], color=str(c))

        g = s.grid

    g = s.grid
    plt.plot( g.xs[g.ilo:g.ihi+1] ,uinit[g.ilo:g.ihi+1],\
            ls=":",color="0.9",zorder=-1 )

    plt.xlabel("$x$")
    plt.ylabel("$u$")
    plt.savefig("fv-burger-sine.pdf")






