import numpy as np
import matplotlib.pyplot as plt
from eulerExact import riemann
from numpy.polynomial.polynomial import polyfit

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

        # flux vector
        # rhov , rhov^2+P , (E+P)v
        self.F = np.zeros((3,Nx))

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

        # interpolated state values close to the boundaries
        self.UIL = np.zeros((3,self.Nx-2))
        self.UIR = np.zeros((3,self.Nx-2))
        # interpolated primitives values close to the boundaries
        self.WIL = np.zeros((3,self.Nx-2))
        self.WIR = np.zeros((3,self.Nx-2))

    def setSod( self ,  x0=0.5 , left_states=[1,0,1] , right_states=[0.1,0.0,0.125],  gamma=1.4 ):
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

    def update_primitive_variables(self):
        self.W[0,:] = self.U[0,:]
        self.W[1,:] = self.U[1,:] / self.U[0,:]
        self.W[2,:] = ( self.gamma - 1.0 ) *\
                ( self.U[2,:] - 0.5 * self.U[1,:]**2 / self.U[0,:] )

    def cons_to_prim( self , U ):
        """ perform a recovery of the primitive variables """
        W = np.zeros((3,self.Nx))
        W[0,:] = self.U[0,:]
        W[1,:] = self.U[1,:] / self.U[0,:]
        W[2,:] = ( self.gamma - 1.0 ) *\
                ( self.U[2,:] - 0.5 * self.U[1,:]**2 / self.U[0,:] )
        return W

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

    def Euler_Flux( self , W ):
        """ compute fluxes for each cell using primitives """
        flux = np.zeros((3,self.Nx))
        flux[0,:] = self.W[0,:] * self.W[1,:]
        flux[1,:] = self.W[0,:] * self.W[1,:]**2 + self.W[2,:]
        flux[2,:] = ( (self.W[2,:] / (self.gamma - 1.0) + \
            0.5*self.W[0,:]*self.W[1,:]**2) + self.W[2,:]) * self.W[1,:]
        return flux

    def get_dt(self ):
        # TODO: implement new eigenvalues ( \lambda_{\pm}=(v\pm c_s)/(1\pm v c_s ) )
        self.cs = self.get_sound_speed( self.W[0,:],self.W[2,:])

        dt = self.cfl * self.dx / \
                np.max([ np.max( np.fabs( self.W[1,:] + self.cs ) ) ,\
                         np.max( np.fabs( self.W[1,:] - self.cs ) )])
        return(dt)

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
        # low order in space
        ap = np.empty( self.Nx - 1 )
        am = np.empty( self.Nx - 1 )
        if self.spatial_order == 1:
            self.get_sound_speed(self.W[0,:],self.W[2,:])

            indices = np.arange( self.Nx - 1 )
            ap = np.maximum( 0 , self.W[1,indices]+self.cs[indices] )
            ap = np.maximum( ap , self.W[1,indices+1] + self.cs[indices+1] )
            am = np.maximum( 0 , -self.W[1,indices] + self.cs[indices] )
            am = np.maximum( am , -self.W[1,indices+1] + self.cs[indices+1])

            self.F[:,:] = self.Euler_Flux( self.W )
            LU = self.getLU( self.U , ap , am )
        elif self.spatial_order != 1:
            FHLL = np.zeros((3,self.Nx-3))
            LU = np.zeros((3,self.Nx))
            UIL, UIR = self.Reconstruct_States()

            WIL = self.cons_to_prim( UIL )
            WIR = self.cons_to_prim( UIR )

            csL = self.get_sound_speed( WIL[0,:], WIL[2,:] )
            csR = self.get_sound_speed( WIR[0,:], WIR[2,:] )

            for i in range(1,self.Nx-2):
                ap[i] = max(0, +(WIL[1,i-1] + csL[i-1]), + (WIR[1,i] + csR[i]) )
                am[i] = max(0, -(WIL[1,i-1] - csL[i-1]), - (WIR[1,i] - csR[i]) )

            FL = self.Euler_Flux( WIL )
            FR = self.Euler_Flux( WIR )

            # need Nx + 1 fluxes because Nx +1 interfaces
            # thats why everything on rhs has len 997 ( we have N=1000 )
            # i added 1 to start and stop of FL and FR because that array actually
            # has same size as U
            # when we take FHLL[1:] and FHLL[:-1], that'll yield len 996
            # which is the Nx - 2 Ng we expected to update all physical cells
            for i in range(3):
                FHLL[i,:] = ( ap[1:-1]*FL[i,1:-2] + am[1:-1]*FR[i,2:-1] - ap[1:-1]*am[1:-1]*( UIR[i,1:] -  UIL[i,:-1]) ) / (ap[1:-1] + am[1:-1])
                LU[i,2:-2] = -( FHLL[i,1:]-FHLL[i,:-1])/self.dx

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
        labels = [r'$\rho$',r'$v$',r'$P$']
        fig, axes = plt.subplots(nrows=3,ncols=1, sharex=True)
        if title == "":
            title = "Sod Shock Tube Simulation"
        title += " time t={}".format(self.tfinal)

        fig.suptitle( title ,fontsize=16)
        for i, axis in enumerate(axes):
            axis.plot( self.x , self.W[i,:], label=labels[i])
            axis.set_title(labels[i])
            axis.grid(True)
            axis.legend()
        plt.xlabel('x')
        return (axes)


def minmod( x , y , z ):
    return( 1./4. * np.fabs( np.sign(x) + np.sign(y)) * \
            (np.sign(x) + np.sign(z)) * \
            min(np.minimum(np.minimum(np.fabs(x),np.fabs(y)),np.fabs(z))))

def compute_l1_error( numerical , exact , deltaX ):
    diff = np.abs( numerical - exact )
    summ = diff.sum()
    return deltaX * summ

def fit_line( xdata , ydata ):
    # fit data to a line
    b, m = polyfit( xdata , ydata , 1 )
    return b , m

def line( x , b , m ):
    return m * x + b

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
    deltaX = ()
    L1 = np.zeros( len(Ns) )
    for i in range(len(Ns)):
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

    logns = np.log( Ns )
    logrs = -np.log( rerr )
    logvs = -np.log( verr )
    logps = -np.log( perr )

    rb , rm = fit_line( logns , logrs )
    vb , vm = fit_line( logns , logvs )
    pb , pm = fit_line( logns , logps )

    print("Density slope:" , rm)
    print("Velocity slope:" , vm)
    print("Pressure slope:" , pm)

    fig = plt.figure()
    plt.plot( logns , line(logrs) , color='red' , label='L1 Error on Density')
    plt.plot( logns , line(logvs) , color='green' , label='L1 Error on Velocity')
    plt.plot( logns , line(logps) , color='blue' , label='L1 Error on Pressure')

    plt.title("Convergence Plot for Sod Shock Tube")
    plt.xlabel("$log(N)$")
    plt.ylabel("$log(L_1)$")
    plt.legend()

def f(x,x0,sigma):
    return(np.where(abs(x-x0)<sigma , (1-((x-x0)/sigma)**2)**2 , 0))

if __name__=="__main__":
    t = 0.1
    e = EulerSolver( 400 , 0.0 , 1.0 , 0.5, time_order=2,spatial_order=1 )
    # e.setSod()
    # e.setSmoothWave()
    rho0 = 1.0; p0 = 0.6; alpha = 0.2; x0=0.5; sigma=0.4
    e.setIsentropicWave(rho0,p0,alpha,f,x0,sigma)
    winit = e.W.copy()
    e.evolve( t )
    axes = e.plot()
    # add the initial configurations on each subplot
    init_labels = ['Initial Density','Initial Velocity','Initial Pressure']
    for i, axis in enumerate(axes):
        axis.plot( e.x , winit[i,:], label=init_labels[i])
        axis.legend()

    # do convergence plot
    plot_convergence(order='high')
    plt.show()
