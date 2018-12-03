from utils import *
import numpy as np
import matplotlib.pyplot as plt
import pdb
class EulerSolver:
    def __init__(self, Nx=10 , a=0.0 , b=1.0 ,gamma=1.4,cfl=0.5):
        self.Nx = Nx
        self.a = a
        self.b = b
        self.gamma = gamma
        self.Ng = 1
        self.t = 0.0
        self.tfinal=0.0
        self.cfl=cfl
        self.grid_size = self.Nx + 2 * self.Ng
        self.dx = (b-a)/float(Nx)
        # rho rhov energy
        self.U = np.zeros((3,Nx))
        # rhov rhov^2+P (E+P)v
        self.F = np.zeros((3,Nx))
        self.v = np.zeros(Nx)
        self.pressure = np.zeros(Nx)
        self.k = np.arange(Nx) + 0.5
        self.spatial_order = 1
        self.time_order=1
        self.cs = np.zeros(Nx)
        # self.ap = np.zeros(Nx)
        # self.am = np.zeros(Nx)
        self.x = a + self.dx * ( self.k )

    def setInitConditions(self, pl=1.0 ,rhol=1.0 , vl=0.0 , pr=0.125 , rhor = 0.1 ,vr = 0.0 ):
        midway_point = int(np.floor(self.Nx/2.))
        self.U[0][:midway_point] = rhol
        self.U[1][:midway_point] = rhol * vl
        self.U[2][:midway_point] = total_energy( pl , rhol , vl , self.gamma )

        self.U[0][midway_point:] = rhor
        self.U[1][midway_point:] = rhor * vr
        self.U[2][midway_point:] = total_energy( pr, rhor , vr , self.gamma )

        self.v[:midway_point] = vl
        self.v[midway_point:] = vr

        self.pressure[:midway_point] = pl
        self.pressure[midway_point:] = pr

        self.update_sound_speed()
        self.update_alphas()

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
        self.cs = cs(self.pressure , self.U[0] , self.gamma)

    def update_velocity(self):
        self.v = self.U[1][:] / self.U[0][:]

    def update_pressure( self ):
        self.pressure = compute_pressure(self.U[0][:],\
                self.U[1][:],self.U[2][:],self.gamma)

    def update_alphas(self):
        def LP( v , sound_speed ):
            return v + sound_speed
        def LM( v , sound_speed ):
            return v - sound_speed
        v = self.v
        self.ap = np.maximum(np.maximum( 0 , LP( v[:-1] , self.cs[:-1] ) ), LP( v[1:] ,self.cs[1:] ) )
        self.am = np.maximum(np.maximum( 0 , -LM( v[:-1] , self.cs[:-1] )) , -LM( v[1:] ,self.cs[1:] ) )
        # for i in range( self.Nx - 1):
            # self.ap[i] = max( 0 , LP( v[i] , self.cs[i] ) , LP( v[i+1] ,self.cs[i+1] ) )
            # self.am[i] = max( 0 , -LM( v[i] , self.cs[i] ) , -LM( v[i+1] ,self.cs[i+1] ) )

    def update_time_RK3(self,dt):
        U1 = np.zeros(self.U.shape)
        U2 = np.zeros(self.U.shape)
        U1 = self.U + dt * self.LU( self.U )
        U2 = (3./4.)*self.U+(1./4.)*U1+(1./4.)*dt*self.LU(U1)
        self.U = (1./3.)*self.U +(2./3.)*U2+(2./3.)*dt*self.LU(U2)

    def forward_euler(self,dt):
        udot = self.LU(self.U)
        # update states
        self.U += dt*udot # update with time integration

    def evolve(self, tfinal):
        # pdb.set_trace()
        self.tfinal = tfinal
        while self.t < tfinal: # while time less than tfinal
            self.update_pressure() # sound speed in each region depends on pressure
            self.update_sound_speed() # or else the alphas are zero on first step
            self.update_alphas() # dt depends on alphas
            dt = self.get_dt()
            if (self.time_order==1):
                self.forward_euler(dt)
            elif ( self.time_order == 3):
                self.update_time_RK3(dt)
            else:
                print("Invalid time order")

            self.t += dt # increment time

    def get_dt(self ):
        dt = self.cfl * self.dx / max(np.maximum(self.am,self.ap))
        if self.t+dt > self.tfinal: # if we're about to overshoot,
            dt = self.tfinal - self.t # don't
        return(dt)

    def LU(self,U):
        self.update_pressure()
        self.update_velocity()

        F = np.zeros((3,self.Nx))
        self.F[0][:] = U[1][:]
        self.F[1][:] = (U[1][:] * U[1][:] / U[0][:]) + self.pressure
        self.F[2][:] = (U[2][:] + self.pressure ) * self.v


        self.update_sound_speed()
        self.update_alphas()
        FL = self.F[:,:-1]
        FR = self.F[:,1:]
        LU = np.zeros((3,self.Nx))

        # using dirichlet boundary conditions by not updating ghost cells
        # low order in space
        if self.spatial_order == 1:
            UL = U[:,:-1]
            UR = U[:,1:]
            FHLL = ( self.ap * FL + self.am * FR - self.ap * self.am * ( UR - UL ))  / ( self.ap + self.am )
            # one ghost cell
            LU[:,1:-1] = -( FHLL[:,1:] - FHLL[:,:-1] )/self.dx
        # high order in space
        elif self.spatial_order == 2:
            cl = U.copy()[:,:-1]
            cr = U.copy()[:,1:]
            theta = 1.0
            for i in range(1,self.Nx-2):
                cl[:,i] = U[:,i] + 0.5 * minmod( theta * ( U[:,i] - U[:,i-1] ),\
                        0.5 * ( U[:,i+1] - U[:,i-1] ) , \
                        theta * ( U[:,i+1] - U[:,i] ) )
                cr[:,i] = U[:,i+1] - 0.5 * minmod( theta * ( U[:,i+1] - U[:,i] ) , \
                        0.5 * ( U[:,i+2] - U[:,i] ) , \
                        theta * ( U[:,i+2] - U[:,i+1] ) )
            FHLL = ( self.ap * FL + self.am * FR - self.ap * self.am * (\
                cr - cl ))  / ( self.ap + self.am )
            # two ghost cells
            LU[:,2:-2] = -( FHLL[:,3:] - FHLL[:,:-3] ) / self.dx
            print("LU:",LU[0])
        else:
            print("Invalid order in space")

        return LU

    def plot(self,title="",color='k'):
        labels = ['Pressure' , 'Density' , 'Velocity']
        primitives = [self.pressure , self.U[0] , self.v ]
        fig, axes = plt.subplots(nrows=3,ncols=1, sharex=True)
        if title == "":
            title = "Sod Shock Tube Simulation"
        title += " time t={}".format(self.tfinal)

        fig.suptitle( title ,fontsize=16)
        for i, axis in enumerate(axes):
            axis.plot( self.x , primitives[i], label=labels[i])
            axis.set_xlabel("x")
            axis.set_title(labels[i])
            axis.grid(True)
            axis.legend()
        return (axes)

    def write_to_file():
        ''' write final solution to file'''


if __name__=="__main__":
    e = EulerSolver(Nx=100 , a=0.0 , b=1.0 ,gamma=1.4,cfl=0.5)
    e.setInitConditions()
    e.evolve( 0.1 )
    e.plot()
