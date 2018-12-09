import numpy as np
import matplotlib.pyplot as plt
from utils import *
from numerical_flux import HLLE_Flux
from eigenvalues import lambdaP, lambdaM
from alphas import alphaM , alphaP
from sound_speed import get_sound_speed
from Reconstruct_States import State_Reconstructor
from newton import newton, fp , dfp 

class EulerSolver:
    def __init__(self, Nx=10 ,  a=0.0 , b=1.0 ,cfl=0.5, spatial_order=1, time_order=1, bc='outflow',gamma=1.4):
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

        self.state_reconstructor = State_Reconstructor(self.time_order , self.spatial_order)

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
        self.gamma = gamma

    def setSod( self ,  x0=0.5 , params=None ):
        """
        x0 - Float , value of x position to be the center of the riemann problem
        left-states - density, velocity , pressure
        right-states - density , velocity , pressure
        gamma - thermodynamic gamma to use for the evolution of fluid
        """
        if params is None:
            print("Need to set sod parameters")
        self.gamma = params['gamma']
        self.W[0,:] = np.where( self.x <= x0 , params['rho_l'] , params['rho_r'] )
        self.W[1,:] = np.where( self.x <= x0 , params['v_l'] , params['v_r'] )
        self.W[2,:] = np.where( self.x <= x0 , params['p_l'] , params['p_r'] )

        self.U[:,:] = self.prim_to_cons( self.W )
    def setIsentropicWave( self , rho0 , p0 , alpha , f , *args ):
        initial_wave = f( self.x , *args )
        rho = rho0 * (1.0 + alpha * initial_wave)
        p = p0 * ( rho / rho0 ) ** self.gamma
        cs = get_sound_speed(rho ,p,self.gamma)
        v = (2. / (self.gamma-1.) ) * (cs - get_sound_speed(rho0,p0),self.gamma)

        self.U[0,:] = rho
        self.U[1,:] = rho * v
        self.U[2,:] = ( p / (self.gamma-1.))+(rho*v*v/2.)
        self.W[0,:] = rho
        self.W[1,:] = v
        self.W[2,:] = p
        self.U[:,:] = self.prim_to_cons( self.W )
    def setSmoothWave( self ):
        self.W[0,:] = np.sin(2 * np.pi * self.x)+2.0
        self.W[1,:] = 0.0
        self.W[2,:] = 1.0
        self.U[:,:] = self.prim_to_cons( self.W )

    def update_conservative_variables_RK3(self,dt):
        Un = self.U
        U1 = np.zeros(self.U.shape)
        U2 = np.zeros(self.U.shape)

        update = dt * self.LU( Un )
        U1[:,self.physical] = Un[:,self.physical] + update
        self.fill_BCs(U1)
        update = dt * self.LU( U1 )
        U2[:,self.physical] = (1./4.)*( 3. * Un[:,self.physical] + U1[:,self.physical] + update )
        self.fill_BCs(U2)
        update = dt * self.LU( U2 )
        self.U[:,self.physical] = \
        (1./3.)* (Un[:,self.physical] + 2. * U2[:,self.physical] + 2. * update )
        self.fill_BCs(self.U)

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
        D = U[0,:]
        S = U[1,:]
        tau = U[2,:]
        ps = np.zeros(U.shape[1])
        p0 = 5
        for i in range(U.shape[1]):
            proot = newton( fp, dfp , p0 , D[i] , S[i],tau[i], self.gamma)
            ps[i] = proot

        vs = S / ( tau + ps + D )
        lors = lorentz_factor( vs )
        rs = D / lors

        if (ps<0).any():
            print("Negative pressure")
            ps = np.where(ps<0,np.fabs(S-tau-D),ps)
            print(ps)
        if ( vs >= 1).any():
            print( "v greater than c")
            print(vs[vs>1])
            vs = np.where( vs>1 , 1e-6 , vs )

        W[2,:] = ps[:]
        W[1,:] = vs[:]
        W[0,:] = rs[:]

        return W

    def prim_to_cons( self , W ):
        """ compute relativistic conserved variables """
        U = np.zeros((3,self.grid_size))
        r = W[0,:]
        v = W[1,:]
        p = W[2,:]
        e = specific_internal_energy( r , p , self.gamma )
        h = specific_enthalpy( r , p , e )
        lorentz = lorentz_factor(v)
        U[0,:] = lorentz * r  # D
        U[1,:] = r * v * h * lorentz * lorentz # Sx
        U[2,:] = r * h * lorentz * lorentz - p - lorentz * r # tau
        return U

    def fill_BCs( self , U=None ):
        if U is None:
            U = self.U
        if self.bc == "periodic":
            U[ : , 0 : self.ilo ] = U[ : , self.ihi - self.Ng + 1 : \
                                           self.ihi + 1 ]
            U[ : , self.ihi + 1 : ] = U[ : , self.ilo : \
                                             self.ilo + self.Ng ]
        if self.bc == "outflow":
            for i in range(3):
                U[ i , 0 : self.ilo ] = U[ i , self.ilo ]
                U[ i , self.ihi + 1 : ] = U[ i , self.ihi ]

    def evolve(self, tfinal):
        self.tfinal=tfinal
        while self.t < tfinal: # while time less than tfinal
            self.cs = get_sound_speed( self.W[0,:] , self.W[2,:], self.gamma )
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
        flux = np.zeros((3,self.Nx+1))
        r = W[0,:]
        v = W[1,:]
        p = W[2,:]
        lorentz = lorentz_factor( v )
        e = specific_internal_energy( r , p , gamma=self.gamma )
        h = specific_enthalpy( r , p , e )
        flux[0,:] = r * v * lorentz
        flux[1,:] = r * h * lorentz * lorentz * v * v + p
        flux[2,:] = r * h * lorentz * lorentz * v - r * v * lorentz
        return flux

    def get_dt(self ):
        dt = self.cfl * self.dx / np.max([ np.max( np.fabs( lambdaP( self.W[1,:] , self.cs ) ) ) ,\
                                           np.max( np.fabs( lambdaM( self.W[1,:] , self.cs ) ) ) ])
        return(dt)

    def LU(self , U=None):
        # using dirichlet boundary conditions by not updating ghost cells
        if U is None:
            U = self.U
        LU = np.zeros((3,self.Nx))
        ap = np.empty( self.Nx+1 )
        am = np.empty( self.Nx+1 )

        UL , UR = self.state_reconstructor.Reconstruct_States( U=self.U, theta=1.5 )

        WL = self.cons_to_prim( UL )
        WR = self.cons_to_prim( UR )
        csL = get_sound_speed( WL[0,:] , WL[2,:] ,self.gamma)
        csR = get_sound_speed( WR[0,:] , WR[2,:] ,self.gamma)

        ap = alphaP( WL , csL , WR , csR )
        am = alphaM( WL , csL , WR , csR )

        FL = self.Physical_Fluxes( WL )
        FR = self.Physical_Fluxes( WR )

        # LU[:,:] = np.where(self.x[self.ilo-1:self.ihi+1]/self.t < am , FL , \
        # np.where(np.logical_and(self.x[self.ilo-1:self.ihi+1]/self.t<ap,self.x[self.ilo-1:self.ihi+1]/self.t>am),\
        # self.HLLE_Flux( UL , UR , FL , FR , am , ap ),FR))

        LU = HLLE_Flux( UL , UR , FL , FR , am , ap ) / self.dx

        return LU

    def plot(self,title="",color='k'):
        labels = [r'$\rho$',r'$v$',r'$P$']
        fig, axes = plt.subplots(nrows=3,ncols=1, sharex=True)
        if title == "":
            title = "Sod Shock Tube Simulation"
        title += " to time t={}".format(self.tfinal)

        fig.suptitle( title ,fontsize=16)
        for i, axis in enumerate(axes):
            axis.plot( self.x[self.physical] , self.W[i,self.physical], label=labels[i])
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

    plt.show()
