#include "EulerSolver.h"
double compute_conserved_rho( float rho ){
    return rho;
}
double compute_rhov( float rho , float v ){
    return rho * v ;
}
double compute_E( float rho , float v , float p , float gamma ){
    return( p / ( gamma - 1. ) + 0.5 * rho * v * v );
}


void EulerSolver::prim_to_cons( float rho , float v , float p , float& u_rho , float& rhov , float& E ){
    



}


    //def prim_to_cons( self , W ):
        //""" compute conserved variables """
        //U = np.zeros((3,self.grid_size))
        //r = W[0,:]
        //v = W[1,:]
        //p = W[2,:]

        //U[0,:] = r
        //U[1,:] = r * v
        //U[2,:] = p / ( self.gamma - 1. ) + 0.5 * r * v * v
        //return U
