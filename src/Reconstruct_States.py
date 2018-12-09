import numpy as np
def minmod( x , y , z ):
    prefactor = (1./4.) * np.abs( np.sign( x ) + np.sign ( y ) )
    minimum = np.minimum( np.abs( x ) , np.abs ( y ) )
    minimum = np.minimum( minimum , np.abs( z ) )
    return( prefactor * ( np.sign( x ) + np.sign( z ) ) * minimum )

class State_Reconstructor( object ):
    def __init__(self , time_order=1,spatial_order=1,method=""):
        self.time_order=time_order
        self.spatial_order=spatial_order
        self.method=method
        if self.spatial_order==1:
            self.Ng=1
        else:
            self.Ng=2
        self.ilo=self.Ng

    def Reconstruct_States( self , U=None , theta=1.5 ):
        self.Nx = U.shape[1]-2*self.Ng
        self.ihi = self.Ng + self.Nx - 1
        self.U = U
        UL = np.zeros((3, self.Nx+1 ) )
        UR = np.zeros((3, self.Nx+1 ) )
        if self.spatial_order==1:
            UL , UR = self.Piecewise_Constant()
        else:
            UL , UR = self.Piecewise_Linear( theta=theta )
        return UL , UR

    def Piecewise_Constant( self ):
        # Godunov piecewise constant
        UL = self.U[:,self.ilo-1:self.ihi+1]
        UR = self.U[:,self.ilo:self.ihi+2]
        return UL , UR

    def Piecewise_Linear( self , theta ):
        UIL = np.zeros((3,self.Nx+1))
        UIR = np.zeros((3,self.Nx+1))
        for i in range( self.Nx+1 ):
            UIL[:,i] = self.U[:,i+1] + 0.5 *\
            minmod( theta * ( self.U[:,i+1] - self.U[:,i]),\
            0.5 * ( self.U[:,i+2] - self.U[:,i] ) ,\
            theta * (self.U[:,i+2] - self.U[:,i+1]))

            UIR[:,i] = self.U[:,i+2] - 0.5 *\
            minmod( theta * ( self.U[:,i+2]-self.U[:,i+1]),\
            0.5 * ( self.U[:,i+3]-self.U[:,i+1] ) ,\
            theta * (self.U[:,i+3]-self.U[:,i+2]))

        UL = UIL
        UR = UIR
        return UL, UR

