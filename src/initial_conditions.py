

class Initial_Conditions( object ):
    def __init__(self, W , U ,  params, method="" ):
        self.W = W
        self.U = U
        if method == 'sod':
            print("Hello world")



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
