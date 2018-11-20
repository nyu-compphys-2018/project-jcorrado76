import numpy as np
import sys


class Grid1d( object ):
    def __init__( self , N , Ng , xmin=0.0 , xmax=1.0 , bc='outflow' ):
        self.N = N
        self.Ng = Ng

        self.xmin = xmin
        self.xmax = xmax

        self.bc = bc

        self.ilo = Ng
        self.ihi = Ng + N - 1

        self.dx = (xmax - xmin) / N
        self.xs = xmin + \
                (np.arange( N + 2 * Ng ) - Ng + 0.5 ) * self.dx

        self.U = np.zeros( N + 2 * Ng , dtype=np.float64 )

    def get_scratch_array( self ):
        """ return a scratch array dimensioned for our grid """

        return( np.zeros( (self.N + 2 * self.Ng ) , dtype=np.float64 ))

    def fill_BCs( self ):
        """ fill ghostcells """

        if self.bc == "periodic":
            # left boundary
            self.U[ 0 : self.ilo ] = self.U[ self.ihi - self.Ng + 1 :\
                                             self.ihi + 1 ]

            # right boundary 
            self.U[ self.ihi + 1 : ] = self.U[ self.ilo : \
                                        self.ilo + self.Ng ]

        elif self.bc == "outflow":
            # left boundary 
            self.U[ 0 : self.ilo ] = self.U[ self.ilo ]

            # right boundary 
            self.U[ self.ihi + 1 : ] = self.U[ self.ihi ]

        else:
            sys.exit("invalid BC")



