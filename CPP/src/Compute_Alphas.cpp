#include "EulerSolver.h"


double EulerSolver::alphaP( const std::vector<float> &rho_l , const std::vector<float> &v_l , const std::vector<float> &p_l , 
                            const std::vector<float> &rho_r , const std::vector<float> &v_r , const std::vector<float> &p_r  , 
                            const std::vector<float> &csL , const std::vector<float> &csR , 
                            std::vector<float> &ap ){

    float left_ap;
    for ( int  i = 0 ; i < Nx+1 ; i++ ){
        left_ap = std::max(0 , lambdaP( v_l[i] , csL[i] ) );
        ap[i] = std::max( left_ap , lambdaP( v_r[i] , csR[i] ) );
    }
}



double EulerSolver::alphaM( const std::vector<float> &rho_l , const std::vector<float> &v_l , const std::vector<float> &p_l , 
                            const std::vector<float> &rho_r , const std::vector<float> &v_r , const std::vector<float> &p_r , 
                            const std::vector<float> &csL , const std::vector<float> &csR , 
                            std::vector<float> &am ){

    float left_am;
    for ( int  i = 0 ; i < Nx+1 ; i++ ){
        left_am = std::max(0 , -lambdaM( v_l[i] , csL[i] ) );
        am[i] = std::max( left_am , -lambdaM( v_r[i] , csR[i] ) );
    }
}

