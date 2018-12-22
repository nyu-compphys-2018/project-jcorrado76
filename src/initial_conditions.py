from prim_to_cons import prim_to_cons
from numpy import where,sin,pi
from sound_speed import get_sound_speed

class Initial_Conditions( object ):
    def __init__( self , W , U , x ):
        self.W = W
        self.U = U
        self.x = x
    def setSod( self , params ,  x0=0.5 ):
        """
        x0 - Float , value of x position to be the center of the riemann problem
        left-states - density, velocity , pressure
        right-states - density , velocity , pressure
        gamma - thermodynamic gamma to use for the evolution of fluid
        """
        gamma = params['gamma']
        self.W[0,:] = where( self.x <= x0 , params['rho_l'] , params['rho_r'] )
        self.W[1,:] = where( self.x <= x0 , params['v_l'] , params['v_r'] )
        self.W[2,:] = where( self.x <= x0 , params['p_l'] , params['p_r'] )

        self.U[:,:] = prim_to_cons( self.W , gamma )
    def setIsentropicWave( self , params , f , *args ):
        alpha = params['alpha'];gamma=params['gamma']
        p0=params['p0'];rho0=params['rho0']
        initial_wave = f( self.x , *args )
        rho = rho0 * (1.0 + alpha * initial_wave)
        p = p0 * ( rho / rho0 ) ** gamma
        cs = get_sound_speed(rho ,p,gamma)
        v = (2. / (gamma-1.) ) * (cs - get_sound_speed(rho0,p0,gamma))

        self.U[0,:] = rho
        self.U[1,:] = rho * v
        self.U[2,:] = ( p / (gamma-1.))+(rho*v*v/2.)
        self.W[0,:] = rho
        self.W[1,:] = v
        self.W[2,:] = p
        self.U[:,:] = prim_to_cons( self.W , gamma )
    def setSmoothWave( self ):
        self.W[0,:] = sin(2 * pi * self.x)+2.0
        self.W[1,:] = 0.0
        self.W[2,:] = 1.0
        self.U[:,:] = prim_to_cons( self.W )
