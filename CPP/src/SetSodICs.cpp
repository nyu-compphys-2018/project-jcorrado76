#include "EulerSolver.h"
void EulerSolver::setSod( float x0=0.5 , float rho_=1. , float v_l=0. , float p_l=1. , float rho_r=0.125 , float v_r=0.0 , float p_r=0.1 , float gam=1.4 ){
    gamma = gam;
    for ( int i = 0 ; i < grid_size ; i++ ){
        if ( x[i] <= x0 ){
            rho[i] = rho_l;
            v[i] = v_l;
            p[i] = p_l;

        }else{
            rho[i] = rho_r;
            v[i] = v_r;
            p[i] = p_r;
        }







    //def setSod( self ,  x0=0.5 , left_states=[1,0,1] , right_states=[0.125,0.0,0.1],  gamma=1.4 ):
        //"""
        //x0 - Float , value of x position to be the center of the riemann problem
        //left-states - density, velocity , pressure
        //right-states - density , velocity , pressure
        //gamma - thermodynamic gamma to use for the evolution of fluid
        //"""
        //self.gamma = gamma
        //for i in range(3):
            //self.W[i,:] = np.where( self.x <= x0 , left_states[i] , right_states[i] )
        //self.U = self.prim_to_cons( self.W )
