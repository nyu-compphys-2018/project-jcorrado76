#include "EulerSolver.h"


void EulerSolver::alphaP( const std::vector<float> &rho_l , const std::vector<float> &v_l , const std::vector<float> &p_l , 
                            const std::vector<float> &rho_r , const std::vector<float> &v_r , const std::vector<float> &p_r  , 
                            const std::vector<float> &csL , const std::vector<float> &csR , 
                            std::vector<float> &ap ){

    float zero = 0.0;
    float left_ap;
    for ( int  i = 0 ; i < Nx+1 ; i++ ){
        left_ap = std::max( zero , lambdaP( v_l[i] , csL[i] ) );
        ap[i] = std::max( left_ap , lambdaP( v_r[i] , csR[i] ) );
    }
}



void EulerSolver::alphaM( const std::vector<float> &rho_l , const std::vector<float> &v_l , const std::vector<float> &p_l , 
                            const std::vector<float> &rho_r , const std::vector<float> &v_r , const std::vector<float> &p_r , 
                            const std::vector<float> &csL , const std::vector<float> &csR , 
                            std::vector<float> &am ){

    float left_am;
    float zero = 0.0;
    for ( int  i = 0 ; i < Nx+1 ; i++ ){
        left_am = std::max(zero , -lambdaM( v_l[i] , csL[i] ) );
        am[i] = std::max( left_am , -lambdaM( v_r[i] , csR[i] ) );
    }
}

