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
        self.grid_size = N + 2 * Ng
        # these xs are the x values at centers of cells
        self.xs = xmin + \
                (np.arange( N + 2 * Ng ) - Ng + 0.5 ) * self.dx
        self.U = np.zeros( self.grid_size , dtype=np.float64 )
        self.physical = slice(self.ilo,self.ihi+1)

    def get_scratch_array( self ):
        """ return a scratch array dimensioned for our grid """
        return( np.zeros( (self.grid_size ) , dtype=np.float64 ))

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
class Grid1d_Euler( Grid1d ):
    def __init__( self , N , Ng , xmin=0.0 , xmax=1.0 , bc='outflow' ):
        # call constructor of parent
        super().__init__(N , Ng , xmin=xmin , xmax=xmax , bc=bc)
        # indices of row for each variable
        self.WRHO = 0
        self.URHO = 0
        self.WV   = 1
        self.URHOV= 1
        self.WP   = 2
        self.UENER= 2
        # number of variables to evolve
        self.NVAR = 3
        # conserved variables, 3 rows, many columns
        self.U = np.zeros( (self.NVAR , self.grid_size ) , dtype=np.float64 )
        # primitive variables
        self.W = np.zeros( (self.NVAR , self.grid_size ) , dtype=np.float64 )
        self.cons_vars = [r'$\rho$',r'$v$',r'$E$']

    def get_scratch_array( self ):
        return( np.zeros( (self.NVAR, self.grid_size ) , dtype=np.float64 ) )

    def fill_BCs( self ):
        """ fill ghostcells """
        # see if 1d burgers
        try:
            NVAR = self.U.shape[0]
        # if not, make sure BCs know
        except:
            NVAR = 1

        if self.bc == "periodic":
            # if 1 variable
            if NVAR == 1:
                # left boundary
                self.U[ 0 : self.ilo ] = self.U [ self.ihi - self.Ng + 1 : \
                                        self.ihi + 1 ]
                # right boundary
                self.U[ self.ihi + 1 : ] = self.U [ self.ilo : \
                                          self.ilo + self.Ng ]

            # if not, multiple variables
            else:
                # left boundary
                self.U[ : , 0 : self.ilo ] = self.U[ : , self.ihi - self.Ng + 1 : \
                                               self.ihi + 1 ]

                # right boundary
                self.U[ : , self.ihi + 1 : ] = self.U[ : , self.ilo : \
                                                 self.ilo + self.Ng ]

        elif self.bc == "outflow":
            # if 1 variable
            if NVAR == 1:
                self.U[ 0 : self.ilo ] = self.U[ self.ilo ]
                self.U[ self.ihi + 1 : ] = self.U[ self.ihi ]
            else:
                # left boundary
                for i in range(NVAR):
                    self.U[ i , 0 : self.ilo ] = self.U[ i , self.ilo ]

                # right boundary
                for i in range(NVAR):
                    self.U[ i , self.ihi + 1 : ] = self.U[ i , self.ihi ]

        else:
            sys.exit("invalid BC")
