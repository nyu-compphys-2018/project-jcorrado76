from grid_1d import Grid1d_Euler
from reconstruction import State_Reconstructor,shifted
import numpy as np

def minmod( x , y , z ):
    return( 1./4. * np.fabs( np.sign(x) + np.sign(y)) * \
            (np.sign(x) + np.sign(z)) * \
            min(np.minimum(np.minimum(np.fabs(x),np.fabs(y)),np.fabs(z))))
class Simulation( object ):
    def __init__( self , grid , CFL=0.8 , gamma=1.4 ):
        self.grid = grid
        self.t = 0.0
        self.gamma=gamma
        self.CFL=CFL

    def compute_timestep( self ):
        return( self.CFL * self.grid.dx / self.max_lambda() )

    def max_lambda( self ):
        """ return maximum relativistic eigenvalue """
        rho = self.grid.U[0,:]
        v = self.grid.U[1,:]/rho
        p = (self.gamma-1.0)*(self.grid.U[2,:]-rho*v*v/2.)
        cs = self.compute_soundspeed( p , rho )
        return( max( (np.abs(v)+cs)/(1-v*cs) ) )

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

    def compute_soundspeed( self , p , rho ):
        return np.sqrt( self.gamma * p / rho )

    # lambdal and lambdar are computed after left and right primitive states are obtained
    def lambdaL( self , WL ):
        """ 1D relativistic maximal eigval """
        v = WL[1,:]
        cs = compute_soundspeed(WL[2,:],WL[0,:])
        return (v + cs)/(1+v*cs)

    def lambdaR( self , WR ):
        """ 1D relativistic minimal eigval """
        v = WR[1,:]
        cs = compute_soundspeed(WR[2,:],WL[0,:])
        return (v - cs)/(1+v*cs)

    def euler_flux( self , U ):
        """ compute basic euler flux from conservative variables """
        g = self.grid
        flux = g.get_scratch_array()
        rho = U[0,:]
        rhov = U[1,:]
        E = U[2,:]
        v = rhov / rho
        p = (self.gamma - 1.) * ( E - rho * v * v / 2. )

        flux[0,:] = rhov
        flux[1,:] = rhov * v + p
        flux[2,:] = ( E + p ) * v
        return flux

    def reconstruct_states( self ):
        """ reconstruct left and right states """
        g = self.grid
        UL = g.get_scratch_array()
        UR = g.get_scratch_array()
        i = slice(0,g.N-2)
        theta=1.5
        for var in range(g.NVAR):
            UL[var,i] = g.U[var,shifted(i,1)] + 0.5 *\
            minmod( theta * ( g.U[var,shifted(i,1)]-g.U[var,i]),\
            0.5 * ( g.U[var,shifted(i,2)]-g.U[var,i] ) ,\
            theta * (g.U[var,shifted(i,2)]-g.U[var,shifted(i,1)]))

            UR[var,i] = g.U[var,shifted(i,1)] - 0.5 *\
            minmod( theta * ( g.U[var,shifted(i,1)]-g.U[var,i]),\
            0.5 * ( g.U[var,shifted(i,2)]-g.U[var,i] ) ,\
            theta * (g.U[var,shifted(i,2)]-g.U[var,shifted(i,1)]))

        return UL,UR

    # TODO: Need to write 3 separate flux functions in here
    # TODO: Need to compute HLL flux in here
    def advective_flux( self , ul , ur ):
        """
        computes flux
        solve the riemann problem given the left and right states
        returns flux at interfaces, given left and right states
        implemented:
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
        unew[ : , g.physical ] = \
                g.U[ : , g.physical ] + \
                dt / g.dx * ( flux[ : , g.physical ] - flux[ : , shifted(g.physical,1) ] )
        return( unew )

    def rk3_update( self , dt , flux ):
        g = self.grid
        unew = g.get_scratch_array()
        unew[ : , g.physical ] = \
                g.U[ : , g.physical ] + \
                dt / g.dx * ( flux[ : , g.physical ] - flux[ : , shifted(g.physical,1) ] )
        return( unew )

    def states(self, dt):
        """ compute the left and right interface states """

        g = self.grid
        # compute the piecewise linear slopes -- 2nd order MC limiter
        # we pick a range of cells that includes 1 ghost cell on either
        # side
        ib = g.ilo-1
        ie = g.ihi+1

        u = g.U

        # this is the MC limiter from van Leer (1977), as given in
        # LeVeque (2002).  Note that this is slightly different than
        # the expression from Colella (1990)

        dc = g.get_scratch_array()
        dl = g.get_scratch_array()
        dr = g.get_scratch_array()

        dc[:,ib:ie+1] = 0.5*(u[:,ib+1:ie+2] - u[:,ib-1:ie  ])
        dl[:,ib:ie+1] = u[:,ib+1:ie+2] - u[:,ib  :ie+1]
        dr[:,ib:ie+1] = u[:,ib  :ie+1] - u[:,ib-1:ie  ]

        # these where's do a minmod()
        d1 = 2.0*np.where(np.fabs(dl) < np.fabs(dr), dl, dr)
        d2 = np.where(np.fabs(dc) < np.fabs(d1), dc, d1)
        ldeltau = np.where(dl*dr > 0.0, d2, 0.0)

        # now the interface states.  Note that there are 1 more interfaces
        # than zones
        ul = g.get_scratch_array()
        ur = g.get_scratch_array()

        # are these indices right?
        #
        #  --+-----------------+------------------+
        #     ^       i       ^ ^        i+1
        #     ur(i)     ul(i+1) ur(i+1)
        #
        ur[:,ib:ie+2] = u[:,ib:ie+2] - \
                      0.5*(1.0 + u[:,ib:ie+2]*dt/self.grid.dx)*ldeltau[:,ib:ie+2]

        ul[:,ib+1:ie+2] = u[:,ib:ie+1] + \
                        0.5*(1.0 - u[:,ib:ie+1]*dt/self.grid.dx)*ldeltau[:,ib:ie+1]

        return ul, ur

    def evolve( self , tmax ):
        self.t = 0.0
        g = self.grid # alias
        while( self.t < tmax ):
            dt = self.compute_timestep()
            if ( self.t + dt > tmax ): # if we're about to overstep
                dt = tmax - self.t     # don't
            self.t += dt

            # rk3 substeps
            # U0 = self.U
            # udot = self.LU()
            # self.U = U0 + dt * udot
            # self.update_primitive_variables()
            # U1 = self.U
            # udot = self.LU()
            # self.U = 3./4. * U0 + 1./4. * U1 + 1./4. * dt * udot
            # self.update_primitive_variables()
            # U2 = self.U
            # udot = self.LU()
            # self.U = 1./3. * U0 + 2./3. * U2 + 2./3. * dt * udot

            # THESE ARE DONE INSIDE LU------------
            # compute interface states
            ul , ur = self.states(dt)

            # solve riemann problem
            flux = self.advective_flux( ul , ur )

            # do conservative update
            unew = self.update( dt , flux )
            self.grid.U[:] = unew[:]

            # TODO: update primitives after conservative update

            g.fill_BCs()

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from initial_conditions import IC_Manager
    xmin = 0.0
    xmax = 1.0
    nx = 256
    ng = 2
    g = Grid1d_Euler( nx , ng , bc='outflow' )

    initial_condition = "rarefaction"
    ICs = IC_Manager( g )
    ICs.set_ICs(initial_condition)
    tmax = 0.1
    CFL = 0.8
    physical = slice(g.ilo,g.ihi+1)
    fig,ax = plt.subplots(3,1,sharex=True)
    s = Simulation( g )

    # s.set_ICs( initial_condition )

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
    plt.savefig("plots/fv-burger-{}.pdf".format( initial_condition ) )
    plt.suptitle("Solution for {} wave".format(initial_condition))
    plt.show()
