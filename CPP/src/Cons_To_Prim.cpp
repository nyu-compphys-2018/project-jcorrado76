#include "EulerSolver.h"

void EulerSolver::cons_to_prim( const std::vector<float> rho , const std::vector<float> rhov , const std::vector<float> E ,
                                   std::vector<float> &vs , std::vector<float> &ps ){

    for ( int i = 0 ; i < grid_size ; i++ ){
        vs[i] = rhov[i] / rho[i] ;
        ps[i] = ( gamma - 1. ) * ( E[i] - 0.5 * rhov[i] * rhov[i] / rho[i] );
    }
}


