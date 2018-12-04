import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

from utils import *
import pdb

def specific_internal_energy( rho , pressure , gamma=1.4 ):
    """ assumes ideal gas law """
    return( pressure / ( rho * ( gamma-1. ) ) )

def specific_enthalpy( rho , pressure , e ):
    return( 1 + e * (pressure / rho) )

def check_if_negative_pressures( pressures ):
    if isinstance(pressures,np.ndarray):
        if (pressures < 0.0).any():
            print("Warning, negative pressure encountered when computing sound speed")
    else:
        if pressures < 0.0:
            print("Negative pressure encountered when computing sound speed")
            print(pressures)

def initialize_animation():
    line.set_data([], [])
    return line,

def animate(t):
    e = EulerSolver(Nx=400 , a=0.0 , b=1.0 , cfl=0.3, time_order=2,spatial_order=2 )
    e.setSod()
    e.evolve(t)
    line.set_data( e.x,e.W[2,:] )
    return line,

class EulerSolver:
    def __init__(self, Nx=10 ,  a=0.0 , b=1.0 ,cfl=0.5, spatial_order=1, time_order=1, bc='outflow'):
        self.Nx = Nx
        self.a = a
        self.b = b
        self.t = 0.0
        self.cfl=cfl
        self.bc=bc
        self.spatial_order = spatial_order
        self.time_order=time_order

        if self.spatial_order == 1:
            self.Ng = 1
        else:
            self.Ng = 2

        self.ilo = self.Ng
        self.ihi = self.Ng + self.Nx - 1

        self.grid_size = self.Nx + 2 * self.Ng
        self.k = np.arange( self.Nx + 2 * self.Ng )
        self.dx = (self.b-self.a)/float(self.Nx)

        # xs[ilo:ihi+1] contains the physical range of xs
        self.x = a + \
                (np.arange( self.Nx + 2 * self.Ng ) - self.Ng + 0.5 ) * self.dx
        self.physical = slice(self.ilo,self.ihi+1)

        # conserved quantities
        # rho , rhov , energy
        self.U = np.zeros((3,self.grid_size))

        # primitives
        # density, velocity, pressure
        self.W = np.zeros((3,self.grid_size))

        # speed of sound
        self.cs = np.zeros(self.grid_size)
        self.gamma = 1.4

    def setSod( self ,  x0=0.5 , left_states=[1,0,1] , right_states=[0.125,0.0,0.1],  gamma=1.4 ):
        """
        x0 - Float , value of x position to be the center of the riemann problem
        left-states - density, velocity , pressure
        right-states - density , velocity , pressure
        gamma - thermodynamic gamma to use for the evolution of fluid
        """
        self.gamma = gamma
        for i in range(3):
            self.W[i,:] = np.where( self.x <= x0 , left_states[i] , right_states[i] )
        self.U = self.prim_to_cons( self.W )

    def setIsentropicWave( self , rho0 , p0 , alpha , f , *args ):
        initial_wave = f( self.x , *args )
        rho = rho0 * (1.0 + alpha * initial_wave)
        p = p0 * ( rho / rho0 ) ** self.gamma
        cs = self.get_sound_speed(rho ,p)
        v = (2. / (self.gamma-1.) ) * (cs - self.get_sound_speed(rho0,p0))

        self.U[0,:] = rho
        self.U[1,:] = rho * v
        self.U[2,:] = ( p / (self.gamma-1.))+(rho*v*v/2.)
        self.W[0,:] = rho
        self.W[1,:] = v
        self.W[2,:] = p

    def setSmoothWave( self ):
        self.W[0,:] = np.sin(2 * np.pi * self.x)+2.0
        self.W[1,:] = 0.0
        self.W[2,:] = 1.0
        self.U[0,:] = self.W[0,:] # set initial density
        self.U[1,:] = self.W[0,:] * self.W[1,:]
        self.U[2,:] = 0.5 * self.W[0,:] * self.W[1,:]**2 + \
                self.W[2,:] / (self.gamma - 1.0)

    def lambdaP( self  , v , cs ):
        return v+cs

    def lambdaM( self , v , cs ):
        return v-cs

    def get_sound_speed(self, r , p):
        """ get relativistic sound speed """
        check_if_negative_pressures( p )
        return np.sqrt(self.gamma * p / r )

    def update_conservative_variables_RK3(self,dt):
        U0 = self.U
        udot = self.LU()
        self.U[:,self.physical] = U0[:,self.physical] + dt * udot[:,:]
        U1 = self.U
        udot = self.LU()
        self.U[:,self.physical] = 3./4. * U0[:,self.physical] + 1./4. * U1[:,self.physical] + \
        1./4. * dt * udot
        U2 = self.U
        udot = self.LU()
        self.U[:,self.physical] = 1./3. * U0[:,self.physical] + 2./3. * U2[:,self.physical] + \
        2./3. * dt * udot

    def update_conservative_variables_forward_euler( self , dt ):
        """
        dt - Float , time step
        """
        udot = self.LU()
        self.U[:,self.physical] += dt * udot

    def update_primitive_variables(self):
        """ update the member variable W """
        self.W[:,:] = self.cons_to_prim( self.U )

    def cons_to_prim( self , U ):
        """ perform a recovery of the primitive variables """
        W = np.zeros(U.shape)
        r = U[0,:]
        rv = U[1,:]
        E = U[2,:]
        W[0,:] = r
        W[1,:] = rv / r
        W[2,:] = ( self.gamma - 1.0 ) * ( E - 0.5 * rv * rv / r )
        return W

    def prim_to_cons( self , W ):
        """ compute conserved variables """
        U = np.zeros((3,self.grid_size))
        r = W[0,:]
        v = W[1,:]
        p = W[2,:]

        U[0,:] = r
        U[1,:] = r * v
        U[2,:] = p / ( self.gamma - 1. ) + 0.5 * r * v * v
        return U

    def fill_BCs( self ):
        if self.bc == "periodic":
            self.U[ : , 0 : self.ilo ] = self.U[ : , self.ihi - self.Ng + 1 : \
                                           self.ihi + 1 ]
            self.U[ : , self.ihi + 1 : ] = self.U[ : , self.ilo : \
                                             self.ilo + self.Ng ]
        if self.bc == "outflow":
            for i in range(3):
                self.U[ i , 0 : self.ilo ] = self.U[ i , self.ilo ]
                self.U[ i , self.ihi + 1 : ] = self.U[ i , self.ihi ]

    def evolve(self, tfinal):
        self.tfinal=tfinal
        while self.t < tfinal: # while time less than tfinal
            if (self.W[2,:]<0).any():
                print("Negative pressure:")
                print(self.W[2,:])
            self.cs = self.get_sound_speed( self.W[0,:] , self.W[2,:] )
            dt = self.get_dt()
            if self.t+dt > tfinal: # if we're about to overshoot,
                dt = tfinal - self.t # don't
            if self.time_order == 1:
                self.update_conservative_variables_forward_euler( dt )
            elif self.time_order != 1:
                self.update_conservative_variables_RK3( dt )
            self.update_primitive_variables()
            self.fill_BCs()
            self.t += dt # increment time

    def Physical_Fluxes( self , W ):
        """ compute fluxes for each cell using primitives
        rhov , rhov^2+P , (E+P)v
        """
        flux = np.zeros(W.shape)
        r = W[0,:]
        v = W[1,:]
        p = W[2,:]
        E = (p / (self.gamma - 1.0)) + 0.5 * r * v * v

        flux[0,:] = r * v
        flux[1,:] = r * v * v + p
        flux[2,:] = (E + p) * v
        return flux

    def HLLE_Flux( self , UL, UR , FL , FR , am , ap ):
        FHLL = np.zeros((3,self.Nx+1))
        Flux_Difference = np.zeros((3,self.Nx))
        FHLL[:,:] = ( ap * FL + am * FR - ap * am * ( UR - UL ) ) / ( ap + am )
        Flux_Difference = -( FHLL[:,1:] - FHLL[:,:-1] ) / self.dx
        return Flux_Difference

    def get_dt(self ):
        dt = self.cfl * self.dx / np.max([ np.max( np.fabs( self.lambdaP( self.W[1,:] , self.cs ) ) ) ,\
                                           np.max( np.fabs( self.lambdaM( self.W[1,:] , self.cs ) ) ) ])
        return(dt)

    def Reconstruct_States(self, U=None , theta=1.5 ):
        """ do a tvd reconstruction using generalized minmod slope limiter """
        # TODO: this is slow because it loops over Nx
        # even though we have 2 ghost cells on either side, the size of
        # interpolated states is Nx-2 and not -4 because we will need only one of the
        # ghost cells on either side to reconstruct a flux at every interface
        if U is None:
            U=self.U
        UL = np.zeros((3,self.Nx+1))
        UR = np.zeros((3,self.Nx+1))

        if self.spatial_order==1:
            # Godunov piecewise constant
            UL = U[:,self.ilo-1:self.ihi+1]
            UR = U[:,self.ilo:self.ihi+2]

        elif self.spatial_order!=1:
            # piecewise linear with minmod slope limiter
            UIL = np.zeros((3,self.Nx+1))
            UIR = np.zeros((3,self.Nx+1))
            for i in range( self.Nx+1 ):
                UIL[:,i] = U[:,i+1] + 0.5 *\
                minmod( theta * ( U[:,i+1] - U[:,i]),\
                0.5 * ( U[:,i+2] - U[:,i] ) ,\
                theta * (U[:,i+2] - U[:,i+1]))

                UIR[:,i] = U[:,i+2] - 0.5 *\
                minmod( theta * ( U[:,i+2]-U[:,i+1]),\
                0.5 * ( U[:,i+3]-U[:,i+1] ) ,\
                theta * (U[:,i+3]-U[:,i+2]))

            UL = UIL
            UR = UIR
        return UL, UR

    def alphaP( self , WL , csL , WR , csR ):
        ap = np.maximum( 0 , self.lambdaP( WL[1,:] , csL ) )
        ap = np.maximum( ap , self.lambdaP( WR[1,:] , csR ) )
        return ap

    def alphaM( self , WL , csL , WR , csR ):
        am = np.maximum( 0 , -self.lambdaM( WL[1,:] , csL ) )
        am = np.maximum( am , -self.lambdaM( WR[1,:] , csR ) )
        return am

    def LU(self):
        # using dirichlet boundary conditions by not updating ghost cells
        LU = np.zeros((3,self.Nx))
        ap = np.empty( self.Nx+1 )
        am = np.empty( self.Nx+1 )

        UL , UR = self.Reconstruct_States( self.U )

        WL = self.cons_to_prim( UL )
        WR = self.cons_to_prim( UR )
        csL = self.get_sound_speed( WL[0,:] , WL[2,:] )
        csR = self.get_sound_speed( WR[0,:] , WR[2,:] )

        ap = self.alphaP( WL , csL , WR , csR )
        am = self.alphaM( WL , csL , WR , csR )

        FL = self.Physical_Fluxes( WL )
        FR = self.Physical_Fluxes( WR )

        # LU[:,:] = np.where(self.x[self.ilo-1:self.ihi+1]/self.t < am , FL , \
        # np.where(np.logical_and(self.x[self.ilo-1:self.ihi+1]/self.t<ap,self.x[self.ilo-1:self.ihi+1]/self.t>am),\
        # self.HLLE_Flux( UL , UR , FL , FR , am , ap ),FR))

        LU = self.HLLE_Flux( UL , UR , FL , FR , am , ap )

        return LU

    def plot(self,title="",color='k'):
        labels = [r'$\rho$',r'$v$',r'$P$']
        fig, axes = plt.subplots(nrows=3,ncols=1, sharex=True)
        if title == "":
            title = "Sod Shock Tube Simulation"
        title += " to time t={}".format(self.tfinal)

        fig.suptitle( title ,fontsize=16)
        for i, axis in enumerate(axes):
            axis.plot( self.x , self.W[i,:], label=labels[i])
            axis.set_title(labels[i])
            axis.grid(True)
            axis.legend()
        plt.xlabel('x')
        return (axes)


if __name__=="__main__":
    tfinal = 0.1
    order = 'high'
    cfl = 0.3
    if order == 'low':
        e = EulerSolver( Nx=400 , a=0.0 , b=1.0 , cfl=cfl, time_order=1 , spatial_order=1 , bc='outflow' )
    if order == 'high':
        e = EulerSolver( Nx=400 , a=0.0 , b=1.0 , cfl=cfl, time_order=3 , spatial_order=2 , bc='outflow' )
    e.setSod()
    # e.setSmoothWave()
    title="Isentropic Wave"
    # rho0 = 1.0; p0 = 0.6; alpha = 0.2; x0=0.5; sigma=0.4
    # e.setIsentropicWave(rho0,p0,alpha,f,x0,sigma)
    # save initial configuration
    winit = e.W.copy()
    # evolve to final time
    e.evolve( tfinal )
    # plot the euler solver and capture resultant axes
    axes = e.plot(title=title)
    # add the initial configurations on each subplot
    init_labels = ['Initial Density','Initial Velocity','Initial Pressure']
    for i, axis in enumerate(axes):
        axis.plot( e.x , winit[i,:], label=init_labels[i],linestyle='dashed',alpha=0.7)
        axis.legend()


    # make animation
    # fig = plt.figure()
    # ax = plt.axes(xlim=(0, 1), ylim=(-2, 2))
    # line, = ax.plot([], [], lw=2)
    # anim = animation.FuncAnimation(fig, animate, init_func=initialize_animation,
    #                                frames=200, interval=20, blit=True)
    # anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    plt.show()
