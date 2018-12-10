import numpy as np
import matplotlib.pyplot as plt
from utils import *
from numerical_flux import HLLE_Flux
from eigenvalues import lambdaP, lambdaM
from alphas import alphaM , alphaP
from sound_speed import get_sound_speed
from Reconstruct_States import State_Reconstructor
from cons_to_prim import cons_to_prim
from prim_to_cons import prim_to_cons
from initial_conditions import  Initial_Conditions

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
            self.title = "First Order in Space"
        else:
            self.Ng = 2
            self.title = "Second Order in Space"
        if self.time_order == 1:
            self.title += " First Order in Time"
        elif self.time_order != 1:
            self.title += " Third Order in Time"

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

        # initial conditions object
        self.IC_manager = Initial_Conditions( self.W , self.U , self.x )

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

    def update_conservative_variables_forward_euler( self , dt ):
        """
        dt - Float , time step
        """
        udot = self.LU()
        self.U[:,self.physical] += dt * udot

    def update_primitive_variables(self):
        """ update the member variable W """
        self.W[:,:] = cons_to_prim( self.U, self.gamma )

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
        self.title += " to time t={}".format(self.tfinal)
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
        if U is None:
            U = self.U
        LU = np.zeros((3,self.Nx))
        ap = np.empty( self.Nx+1 )
        am = np.empty( self.Nx+1 )

        theta = 1.5
        UL , UR = self.state_reconstructor.Reconstruct_States( U=self.U,theta=theta )
        W = cons_to_prim( U ,self.gamma)
        WL , WR = self.state_reconstructor.Reconstruct_States( W , theta = theta )
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
        fig.suptitle( self.title ,fontsize=16)
        for i, axis in enumerate(axes):
            axis.plot( self.x[self.physical] , self.W[i,self.physical], label=labels[i])
            axis.set_title(labels[i])
            axis.grid(True)
            axis.legend()
        plt.xlabel('x')
        return (axes)
