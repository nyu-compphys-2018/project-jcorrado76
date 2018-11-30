import numpy as np
class IC_Manager( object ):
    def __init__(self , grid=None , type='tophat' ):
        self.grid = grid

    def set_ICs( self , type='tophat' ):
        if type == 'tophat':
            self.grid.U[ : , np.logical_and( self.grid.xs >= 0.333, \
                    self.grid.xs <= 0.666 ) ] = 1.0

        elif type == "sine":
            # initialize all state variables to 1.0
            self.grid.U[ : , : ] = 1.0

            # indices are the middle third of self.grid
            index = np.logical_and( self.grid.xs >= 0.333, \
                    self.grid.xs <= 0.666 )

            self.grid.U[ : , index ] += \
                    0.5*np.sin(2.0*np.pi*(self.grid.xs[index]-0.333)/0.333)

        elif type == "rarefaction":
            self.grid.U[:] = 1.0
            self.grid.U[ :, self.grid.xs > 0.5 ]  = 2.0

        elif type == "sod":
            self.setSod()

    def setSod( self , left_states = {'rho_l':1.0,'u_l':0.0,'p_l':1.0} ,right_states = {'rho_r':0.125,'u_r':0.0,'p_r':0.1} , gamma =1.4 ):
        # this assumes ideal gas eos to compute E_l
        x0 = (self.grid.xmax+self.grid.xmin)/2.
        print(x0)
        # initial density
        self.grid.U[0,:] = np.where( self.grid.xs < x0 , left_states['rho_l'],right_states['rho_r'] )
        # initial momentum density
        rhov_l = left_states['rho_l'] * left_states['u_l']
        rhov_r = right_states['rho_r'] * right_states['u_r']
        self.grid.U[1,:] = np.where( self.grid.xs < x0 , rhov_l , rhov_r )
        # initial energy density
        E_l = left_states['p_l']/(gamma-1.)+0.5*left_states['rho_l']*left_states['u_l']*left_states['u_l']
        E_r = right_states['p_r']/(gamma-1.)+0.5*right_states['rho_r']*right_states['u_r']*right_states['u_r']
        self.grid.U[2,:] = np.where( self.grid.xs < x0 , E_l, E_r )

        self.grid.W[0,:] = np.where( self.grid.xs < x0 , left_states['rho_l'],right_states['rho_r'])
        self.grid.W[1,:] = np.where( self.grid.xs < x0 , left_states['u_l'],right_states['u_r'])
        self.grid.W[2,:] = np.where( self.grid.xs < x0 , left_states['p_l'],right_states['p_r'])
