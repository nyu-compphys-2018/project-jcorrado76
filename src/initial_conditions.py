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
            self.grid.U
