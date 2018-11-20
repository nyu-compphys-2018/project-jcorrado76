import numpy as np
import matplotlib.pyplot as plt
class EulerSolver:
    def __init__(self, Nx=10 , a=0.0 , b=1.0 ,cfl=0.5, spatial_order=1,\
            time_order=1):
        self.Nx = Nx
        self.a = a
        self.b = b
        self.Ng = 1
        self.t = 0.0
        self.cfl=cfl

        self.grid_size = self.Nx + 2 * self.Ng
        self.k = np.arange(Nx) + 0.5
        self.dx = (self.b-self.a)/float(self.Nx)
        self.x = a + self.dx * ( self.k )
        # indices of points grid

        # conserved quantities
        # rho
        # rhov
        # energy
        self.U = np.zeros((3,Nx))

        # flux vector
        # rhov
        # rhov^2+P
        # (E+P)v
        self.F = np.zeros((3,Nx))

        # primitives
        # density
        # velocity
        # pressure
        self.W = np.zeros((3,Nx))

        # speed of sound
        self.cs = np.zeros(Nx)
        self.gamma = 0

        # order in space
        self.spatial_order = spatial_order
        # order in time
        self.time_order=time_order

        # interpolated state values close to the boundaries
        self.UIL = np.zeros((3,self.Nx-2))
        self.UIR = np.zeros((3,self.Nx-2))
        # interpolated primitives values close to the boundaries
        self.WIL = np.zeros((3,self.Nx-2))
        self.WIR = np.zeros((3,self.Nx-2))

    def setSod( self ,  x0=0.5 , left_states=[1,0,1] , \
            right_states=[0.1,0.0,0.125],  gamma=1.4 ):
        """
        x0 - Float , value of x position to be the center of the riemann problem
        left-states - List , values of pressure velocity and pressure for the
        left hand states
        right-states - List , values of pressure velocity and pressure for the
        right-hand states
        gamma - thermodynamic gamma to use for the evolution of fluid
        """

        self.gamma = gamma
        for i in range(self.Nx):
            if self.x[i] <= x0:
                self.W[0,i] = left_states[0] # set density
                self.W[1,i] = left_states[1] # set velocity
                self.W[2,i] = left_states[2] # set pressure
            else:
                self.W[0,i] = right_states[0] # set density
                self.W[1,i] = right_states[1] # set velocity
                self.W[2,i] = right_states[2] # set pressure
        self.U[0,:] = self.W[0,:] # set initial density
        self.U[1,:] = self.W[0,:] * self.W[1,:]
        self.U[2,:] = 0.5 * self.W[0,:] * self.W[1,:]**2 + \
                self.W[2,:] / (self.gamma - 1.0)

    def setInitialWave( self , rho0 , p0 , alpha , f , *args ):
        initial_wave = f( self.x , *args )
        rho = isentropic_initial_rho( rho0 , alpha , initial_wave )
        p = isentropic_initial_P( rho , rho0 , p0 , self.gamma )
        self.cs = cs( p , rho , self.gamma )
        self.v = isentropic_initial_v( rho , rho0 , p ,p0 , self.gamma)

        self.U[0,:] = rho
        self.U[1,:] = rho * self.v
        self.U[2,:] = total_energy( p , rho , self.v , self.gamma)

    def update_sound_speed(self):
        self.cs = np.sqrt( self.gamma * self.W[2,:] / self.W[0,:] )

    def update_conservative_variables_RK3(self,dt):
        U0 = self.U
        udot = self.LU()
        self.U = U0 + dt * udot
        self.update_primitive_variables()
        U1 = self.U
        udot = self.LU()
        self.U = 3./4. * U0 + 1./4. * U1 + 1./4. * dt * udot
        self.update_primitive_variables()
        U2 = self.U
        udot = self.LU()
        self.U = 1./3. * U0 + 2./3. * U2 + 2./3. * dt * udot

    def update_conservative_variables_forward_euler( self , dt , udot ):
        """
        dt - Float , time step
        udot - 3 by Nx array of the updates for the conservative variables
        """
        self.U += dt * udot

    def update_primitive_variables(self):
        self.W[0,:] = self.U[0,:]
        self.W[1,:] = self.U[1,:] / self.U[0,:]
        self.W[2,:] = ( self.gamma - 1.0 ) *\
                ( self.U[2,:] - 0.5 * self.U[1,:]**2 / self.U[0,:] )

    def evolve(self, tfinal):
        self.tfinal=tfinal
        while self.t < tfinal: # while time less than tfinal
            dt = self.get_dt()
            if self.t+dt > tfinal: # if we're about to overshoot,
                dt = tfinal - self.t # don't
            # if order in time is 1
            if self.time_order == 1:
                udot = self.LU()
                # use forward euler
                self.update_conservative_variables_forward_euler( dt , udot )
                self.update_primitive_variables()
            # if higher order in time
            elif self.time_order != 1:
                # use RK3
                self.update_conservative_variables_RK3(dt)
                self.update_primitive_variables()
            self.t += dt # increment time

    def get_dt(self ):
        self.cs = np.sqrt( self.gamma * self.W[2,:] / self.W[0,:] )
        dt = self.cfl * self.dx / \
                np.max([ np.max( np.fabs( self.W[1,:] + self.cs ) ) ,\
                         np.max( np.fabs( self.W[1,:] - self.cs ) )])
        return(dt)

    def LU(self):
        # using dirichlet boundary conditions by not updating ghost cells
        # low order in space
        ap = np.empty( self.Nx - 1 )
        am = np.empty( self.Nx - 1 )
        if self.spatial_order == 1:
            self.update_sound_speed()
            for i in range( self.Nx - 1 ):
                ap[i] = max( 0 , +(self.W[1,i] + self.cs[i] ),\
                        +(self.W[1,i+1] + self.cs[i+1] ) )
                am[i] = max( 0 , -(self.W[1,i] - self.cs[i] ),\
                        -(self.W[1,i+1] - self.cs[i+1] ) )

            self.F[0,:] = self.W[0,:] * self.W[1,:]
            self.F[1,:] = self.W[0,:] * self.W[1,:]**2 + self.W[2,:]
            self.F[2,:] = ( (self.W[2,:] / (self.gamma - 1.0) + \
                0.5*self.W[0,:]*self.W[1,:]**2) + self.W[2,:]) * self.W[1,:]
            LU = self.getLU( self.U , ap , am )
        # high order in space
        elif self.spatial_order != 1:
            theta=1.5
            for i in range( self.Nx-2 ):
                self.UIL[:,i] = self.U[:,i+1] + 0.5 *\
                minmod( theta * ( self.U[:,i+1]-self.U[:,i]),\
                0.5 * ( self.U[:,i+2]-self.U[:,i] ) ,\
                theta * (self.U[:,i+2]-self.U[:,i+1]))

                self.UIR[:,i] = self.U[:,i+1] - 0.5 *\
                minmod( theta * ( self.U[:,i+1]-self.U[:,i]),\
                0.5 * ( self.U[:,i+2]-self.U[:,i] ) ,\
                theta * (self.U[:,i+2]-self.U[:,i+1]))

            # need to update interpolated values of primitives at boundaries
            self.WIL[0,:] = self.UIL[0,:]
            self.WIL[1,:] = self.UIL[1,:] / self.UIL[0,:]
            self.WIL[2,:] = ( self.gamma - 1.0 ) *\
                    ( self.UIL[2,:] - 0.5 * self.UIL[1,:]**2 / self.UIL[0,:] )
            self.WIR[0,:] = self.UIR[0,:]
            self.WIR[1,:] = self.UIR[1,:] / self.UIR[0,:]
            self.WIR[2,:] = ( self.gamma - 1.0 ) *\
                    ( self.UIR[2,:] - 0.5 * self.UIR[1,:]**2 / self.UIR[0,:] )

            self.csL = (self.gamma * self.WIL[2,:]/self.WIL[0,:])**0.5
            self.csR = (self.gamma * self.WIR[2,:]/self.WIR[0,:])**0.5
            # 2 ghost cells on either side
            # update ap and am at k * x + 0.5
            for i in range(1,self.Nx-2):
                ap[i] = max(0, +(self.WIL[1,i-1]+self.csL[i-1]), +(self.WIR[1,i]+self.csR[i]) )
                am[i] = max(0, -(self.WIL[1,i-1]-self.csL[i-1]), -(self.WIR[1,i]-self.csR[i]) )

            #Flux for the left (similar to rL,vL,PL) and right states (similar to rR,vR,pR) at boundaries
            F1L = self.WIL[0,:]*self.WIL[1,:]
            F2L = self.WIL[0,:]*self.WIL[1,:]**2 + self.WIL[2,:]
            F3L = ( (self.WIL[2,:]/(self.gamma - 1) +\
                    0.5*self.WIL[0,:]*self.WIL[1,:]**2) + self.WIL[2,:] )\
            *self.WIL[1,:]

            F1R = self.WIR[0,:]*self.WIR[1,:]
            F2R = self.WIR[0,:]*self.WIR[1,:]**2 + self.WIR[2,:]
            F3R = ( (self.WIR[2,:]/(self.gamma - 1) +\
                    0.5*self.WIR[0,:]*self.WIR[1,:]**2) + self.WIR[2,:] )\
            *self.WIR[1,:]

            #Right and left states at the boundaries
            UIL[0,:]=self.WIL[0,:]
            UIL[1,:]=self.WIL[0,:]*self.WIL[1,:]
            UIL[2,:]=self.WIL[2,:]/(self.gamma-1)\
                    +0.5*self.WIL[0,:]*self.WIL[1,:]**2
            UIR[0,:]=self.WIR[0,:]
            UIR[1,:]=self.WIR[0,:]*self.WIR[1,:]
            UIR[2,:]=self.WIR[2,:]/(self.gamma-1)\
                    +0.5*self.WIR[0,:]*self.WIR[1,:]**2

            #Calculating FHLL at the boundaries 1+0.5, 2+0.5, ... Nx-3+0.5
            FHLL1 = ( ap[1:-1]*F1L[:-1] + am[1:-1]*F1R[1:] - ap[1:-1]*am[1:-1]*(u1R[1:] - u1L[:-1]) ) / (ap[1:-1] + am[1:-1])
            FHLL2 = ( ap[1:-1]*F2L[:-1] + am[1:-1]*F2R[1:] - ap[1:-1]*am[1:-1]*(u2R[1:] - u2L[:-1]) ) / (ap[1:-1] + am[1:-1])
            FHLL3 = ( ap[1:-1]*F3L[:-1] + am[1:-1]*F3R[1:] - ap[1:-1]*am[1:-1]*(u3R[1:] - u3L[:-1]) ) / (ap[1:-1] + am[1:-1])

            LU = np.zeros((3,self.Nx))

            #The first and last 2 states are fixed, change is in cells 2,3,...Nx-3
            #from fluxes at boundaries of each cell (1+0.5, 2+0.5) , (2+0.5 , 3+0.5) ... (Nx-2+0.5, Nx-3+0.5) respectively
            LU[0,2:-2]=-(FHLL1[1:]-FHLL1[:-1])/self.dx
            LU[1,2:-2]=-(FHLL2[1:]-FHLL2[:-1])/self.dx
            LU[2,2:-2]=-(FHLL3[1:]-FHLL3[:-1])/self.dx
        return LU

    def getLU( self , U , ap , am ):
        FL = self.F[:,:-1]
        FR = self.F[:,1:]
        UL = U[:,:-1]
        UR = U[:,1:]
        FHLL = ( ap * FL + am * FR - ap * am * ( UR - UL )) / (ap+am)
        LU = np.zeros((3,self.Nx))
        LU[:,1:-1] = -(FHLL[:,1:]-FHLL[:,:-1]) / self.dx
        return LU

    def plot(self,title="",color='k'):
        labels = ['Density','Velocity','Pressure']
        fig, axes = plt.subplots(nrows=3,ncols=1, sharex=True)
        if title == "":
            title = "Sod Shock Tube Simulation"
        title += " time t={}".format(self.tfinal)

        fig.suptitle( title ,fontsize=16)
        for i, axis in enumerate(axes):
            axis.plot( self.x , self.W[i,:], label=labels[i])
            axis.set_xlabel("x")
            axis.set_title(labels[i])
            axis.grid(True)
            axis.legend()
        return (axes)

def minmod( x , y , z ):
    return( 1./4. * np.fabs( np.sign(x) + np.sign(y)) * \
            (np.sign(x) + np.sign(z)) * \
            min(np.minimum(np.minimum(np.fabs(x),np.fabs(y)),np.fabs(z))))

if __name__=="__main__":
    e = EulerSolver( 1000 , 0.0 , 1.0 , 0.5, time_order=2,spatial_order=1 )
    e.setSod()
    e.evolve(0.1)
    e.plot()
    plt.show()
