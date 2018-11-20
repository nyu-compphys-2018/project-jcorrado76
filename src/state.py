

class State( object ):
    """ a simple object to hold a primitive variable state"""
    def __init__( self , p=1.0 , u=0.0 , rho=1.0 ):
        self.p = p
        self.u = u
        self.rho = rho

    def __str__(self):
        return "rho: {}\nu: {}\np: {}"\
                .format( self.rho , self.u , self.p )

