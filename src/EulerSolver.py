import numpy as np
import matplotlib.pyplot as plt
from eulerExact import riemann
from utils import *

c=3e8

def plot_convergence(order='low'):
    """
    this function plots the convergence rate of then
    scheme
    """
    Ns = [ 32 , 64 , 128 , 256 , 512 ]
    t=0.1
    x0=0.5
    a=0.0
    b=1.0
    rhol = 1.0; vl = 0.0; pl = 1.0
    rhor = 0.1; vr = 0.0; pr = 0.125
    gamma=1.4

    rerr = []
    verr = []
    perr = []

    L1 = np.zeros( len(Ns) )
    for i in range(len(Ns)):
        deltaX = (b-a)/Ns[i]
        if order=='low':
            e = EulerSolver(Nx = Ns[i],time_order=1,spatial_order=1)
        else:
            e = EulerSolver(Nx = Ns[i],time_order=2,spatial_order=2)
        e.setSod()
        e.evolve(t)
        xexact , rexact , vexact , pexact = riemann( a , b , x0 , Ns[i] , t ,rhol ,vl , pl , rhor , vr , pr , gamma )
        rerr.append( compute_l1_error( e.W[0,:] , rexact , deltaX ))
        verr.append( compute_l1_error( e.W[1,:] , vexact , deltaX ))
        perr.append( compute_l1_error( e.W[2,:] ,   pexact , deltaX ))

    print("L1 Errors on Density",rerr)
    print("L1 Errors on Velocity",verr)
    print("L1 Errors on Pressure",perr)

    logns = np.log( Ns )
    logrs = np.log( rerr )
    logvs = np.log( verr )
    logps = np.log( perr )

    rb , rm = fit_line( logns , logrs )
    vb , vm = fit_line( logns , logvs )
    pb , pm = fit_line( logns , logps )

    print("Density slope:" , rm)
    print("Velocity slope:" , vm)
    print("Pressure slope:" , pm)

    fig = plt.figure()

    plt.plot( logns, logrs , color='red' , marker='o',label="L1 Error on Density")
    plt.plot( logns, logvs , color='green' , marker='o',label="L1 Error on Velocity")
    plt.plot( logns, logps , color='blue' , marker='o',label="L1 Error on Pressure")

    plt.plot( logns , line(logns , rb , rm ) , linestyle='dashed', color='red' , label='L1 Error Fit on Density')
    plt.plot( logns , line(logns , vb , vm ) , linestyle='dashed',color='green' , label='L1 Error Fit on Velocity')
    plt.plot( logns , line(logns , pb , pm ) , linestyle='dashed',color='blue' , label='L1 Error Fit on Pressure')

    plt.title("Convergence Plot for Sod Shock Tube")
    plt.xlabel("$log(N)$")
    plt.ylabel("$log(L_1)$")
    plt.legend()

class EulerSolver:
    def __init__(self, Nx=10 , a=0.0 , b=1.0 ,cfl=0.5, spatial_order=1, time_order=1):
        self.Nx = Nx
        self.a = a
        self.b = b
        self.Ng = 1
        self.t = 0.0
        self.cfl=cfl

        self.grid_size = self.Nx + 2 * self.Ng
        self.k = np.arange(Nx) + 0.5
        self.dx = (self.b-self.a)/float(self.Nx)
        # x values of bin centers
        self.x = a + self.dx * ( self.k )
        # indices of points grid

        # conserved quantities
        # rho , rhov , energy
        self.U = np.zeros((3,Nx))

        # primitives
        # density, velocity, pressure
        self.W = np.zeros((3,Nx))

        # speed of sound
        self.cs = np.zeros(Nx)
        self.gamma = 1.4

        # order in space
        self.spatial_order = spatial_order
        # order in time
        self.time_order=time_order

    def setSod( self ,  x0=0.5 , left_states=[1,0,1] , right_states=[0.1,0.0,0.125],  gamma=1.4 ):
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
        self.cs = self.get_sound_speed(rho ,p)
        v = (2. / (self.gamma-1.) ) * (self.cs - self.get_sound_speed(rho0,p0))

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
        return (v+cs)/(1+v*cs/c**2)

    def lambdaM( self , v , cs ):
        return (v-cs)/(1-v*cs/c**2)

    def lorentz( self ):
        """ relativistic lorentz factor """
        return 1./np.sqrt(1.-self.W[1,:]*self.W[1,:]/c**2)

def get_sound_speed(self, r , p):
        # let me know if encountering negative pressures
        if isinstance(p,np.ndarray):
            if (p < 0.0).any():
                print("Warning, negative pressure encountered when computing sound speed")
        else:
            if p < 0.0:
                print("Negative pressure encountered when computing sound speed")
        specific_internal_energy = self.specific_internal_energy( r , p )
        h = 1. + specific_internal_energy + (p/r)
        return np.sqrt( ((self.gamma-1.)/h)*(specific_internal_energy + (p / r )) )



    def get_sound_speed(self, r , p):
        # let me know if encountering negative pressures
        if isinstance(p,np.ndarray):
            if (p < 0.0).any():
                print("Warning, negative pressure encountered when computing sound speed")
        else:
            if p < 0.0:
                print("Negative pressure encountered when computing sound speed")
        return np.sqrt( self.gamma * p / r )

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

    # TODO: make this function call cons_to_prim
    def update_primitive_variables(self):
        """ update the member variable W """
        self.W[0,:] = self.U[0,:]
        self.W[1,:] = self.U[1,:] / self.U[0,:]
        self.W[2,:] = ( self.gamma - 1.0 ) *\
                ( self.U[2,:] - 0.5 * self.U[1,:]**2 / self.U[0,:] )

    # need to perform newton raphson and recover relativistic primitives
    def cons_to_prim( self , U ):
        """ perform a recovery of the primitive variables """
        W = np.zeros(U.shape)
        W[0,:] = U[0,:]
        W[1,:] = U[1,:] / U[0,:]
        W[2,:] = ( self.gamma - 1.0 ) *\
                ( U[2,:] - 0.5 * U[1,:]**2 / U[0,:] )
        return W

    # TODO: make sure relativistic conserved variables are correct D,S,tau
    def prim_to_cons( self , W ):
        """ compute relativistic conserved variables """
        U = np.zeros((3,self.Nx))
        specific_internal_energy = W[2,:]/(W[0,:]*(self.gamma-1.))
        h = 1. + specific_internal_energy + (W[2,:]/W[0,:])
        U[0,:] = self.lorentz()*W[0,:] # set initial density
        U[1,:] = W[0,:] * W[1,:] * h * self.lorentz()**2
        U[2,:] = W[0,:]*h*self.lorentz()**2-W[2,:]-self.lorentz()*W[0,:]
        return U

    def evolve(self, tfinal):
        self.tfinal=tfinal
        while self.t < tfinal: # while time less than tfinal
            dt = self.get_dt()
            if self.t+dt > tfinal: # if we're about to overshoot,
                dt = tfinal - self.t # don't
            if self.time_order == 1:
                udot = self.LU()
                # use forward euler
                self.update_conservative_variables_forward_euler( dt , udot )
                self.update_primitive_variables()
            elif self.time_order != 1:
                # use RK3
                self.update_conservative_variables_RK3(dt)
                self.update_primitive_variables()
            self.t += dt # increment time

    # TODO: make relativistic physical fluxes
    def Euler_Flux( self , W ):
        """ compute fluxes for each cell using primitives
        rhov , rhov^2+P , (E+P)v
        """
        flux = np.zeros((3,self.Nx))
        flux[0,:] = self.W[0,:] * self.W[1,:]
        flux[1,:] = self.W[0,:] * self.W[1,:]**2 + self.W[2,:]
        flux[2,:] = ( (self.W[2,:] / (self.gamma - 1.0) + \
            0.5*self.W[0,:]*self.W[1,:]**2) + self.W[2,:]) * self.W[1,:]
        return flux

    # TODO: need obvious indices here
    def HLLE_Flux( self , UL, UR , FL , FR , am , ap ):
        # need Nx + 1 fluxes because Nx +1 interfaces
        # thats why everything on rhs has len 997 ( we have N=1000 )
        # i added 1 to start and stop of FL and FR because that array actually
        # has same size as U
        # when we take FHLL[1:] and FHLL[:-1], that'll yield len 996
        # which is the Nx - 2 Ng we expected to update all physical cells
        FHLL = np.zeros((3,self.Nx-3))
        LU = np.zeros((3,self.Nx))
        for i in range(3):
            FHLL[i,:] = ( ap[1:-1]*FL[i,1:-2] + am[1:-1]*FR[i,2:-1] - ap[1:-1]*am[1:-1]*( UR[i,1:] -  UL[i,:-1]) ) / (ap[1:-1] + am[1:-1])
            LU[i,2:-2] = -( FHLL[i,1:]-FHLL[i,:-1])/self.dx
        return LU

    def get_dt(self ):
        # TODO: implement new eigenvalues ( \lambda_{\pm}=(v\pm c_s)/(1\pm v c_s ) )
        self.cs = self.get_sound_speed( self.W[0,:],self.W[2,:])

        dt = self.cfl * self.dx / \
                np.max([ np.max( np.fabs( self.lambdaP(self.W[1,:] , self.cs) ) ) ,\
                         np.min( np.fabs( self.lambdaM(self.W[1,:] , self.cs) ) )])
        return(dt)

    #TODO: implement obvious indices
    def Reconstruct_States(self, theta=1.5 ):
        """ do a tvd reconstruction using generalized minmod slope limiter """
        # TODO: this is slow because it loops over Nx
        # even though we have 2 ghost cells on either side, the size of
        # interpolated states is Nx-2 and not -4 because we will need only one of the
        # ghost cells on either side to reconstruct a flux at every interface
        UIL = np.zeros((3,self.Nx-2))
        UIR = np.zeros((3,self.Nx-2))
        for i in range( self.Nx-2 ):
            UIL[:,i] = self.U[:,i+1] + 0.5 *\
            minmod( theta * ( self.U[:,i+1]-self.U[:,i]),\
            0.5 * ( self.U[:,i+2]-self.U[:,i] ) ,\
            theta * (self.U[:,i+2]-self.U[:,i+1]))

            UIR[:,i] = self.U[:,i+1] - 0.5 *\
            minmod( theta * ( self.U[:,i+1]-self.U[:,i]),\
            0.5 * ( self.U[:,i+2]-self.U[:,i] ) ,\
            theta * (self.U[:,i+2]-self.U[:,i+1]))
        return UIL, UIR

    def LU(self):
        # using dirichlet boundary conditions by not updating ghost cells
        LU = np.zeros((3,self.Nx))
        ap = np.empty( self.Nx - 1 )
        am = np.empty( self.Nx - 1 )
        if self.spatial_order == 1:
            cs = self.get_sound_speed(self.W[0,:],self.W[2,:])
            indices = np.arange( self.Nx - 1 )
            ap = np.maximum( 0 , self.lambdaP(self.W[1,indices] , cs[indices]) )
            ap = np.maximum( ap , self.lambdaP(self.W[1,indices+1] , cs[indices+1]) )
            am = np.maximum( 0 , -self.lambdaM(self.W[1,indices] , cs[indices]) )
            am = np.maximum( am , -self.lambdaM(self.W[1,indices+1] , cs[indices+1]))

            F = self.Euler_Flux( self.W )
            LU[:,:] = self.getLU( self.U , F , ap , am )

        elif self.spatial_order != 1:
            LU = np.zeros((3,self.Nx))
            UIL, UIR = self.Reconstruct_States()
            WIL = self.cons_to_prim( UIL )
            WIR = self.cons_to_prim( UIR )
            csL = self.get_sound_speed( WIL[0,:], WIL[2,:] )
            csR = self.get_sound_speed( WIR[0,:], WIR[2,:] )
            for i in range(1,self.Nx-2):
                ap[i] = max(0, self.lambdaP(WIL[1,i-1] , csL[i-1]), self.lambdaP(WIR[1,i] , csR[i]) )
                am[i] = max(0, -self.lambdaM(WIL[1,i-1] , csL[i-1]), -self.lambdaM(WIR[1,i] , csR[i]) )
            # compute physical fluxes
            FL = self.Euler_Flux( WIL )
            FR = self.Euler_Flux( WIR )
            LU[:,:] = self.HLLE_Flux( UIL, UIR , FL , FR , am , ap )

        return LU

    def getLU( self , U , F , ap , am ):
        UL = U[:,:-1]
        UR = U[:,1:]
        FL = F[:,:-1]
        FR = F[:,1:]
        FHLL = ( ap * FL + am * FR - ap * am * ( UR - UL )) / (ap+am)
        LU = np.zeros((3,self.Nx))
        LU[:,1:-1] = -(FHLL[:,1:]-FHLL[:,:-1]) / self.dx
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
    # final time
    t = 0.2
    # initialize euler solver object
    e = EulerSolver( Nx=400 , a=0.0 , b=1.0 , cfl=0.3, time_order=2,spatial_order=2 )
    # set initial conditions
    e.setSod()
    # e.setSmoothWave()
    title="Sod"
    # rho0 = 1.0; p0 = 0.6; alpha = 0.2; x0=0.5; sigma=0.4
    # e.setIsentropicWave(rho0,p0,alpha,f,x0,sigma)
    # save initial configuration
    winit = e.W.copy()
    # evolve to final time
    e.evolve( t )
    # plot the euler solver and capture resultant axes
    axes = e.plot(title=title)
    # add the initial configurations on each subplot
    init_labels = ['Initial Density','Initial Velocity','Initial Pressure']
    for i, axis in enumerate(axes):
        axis.plot( e.x , winit[i,:], label=init_labels[i],linestyle='dashed',alpha=0.7)
        axis.legend()

    # do convergence plot
    # plot_convergence(order='high')
    plt.show()
