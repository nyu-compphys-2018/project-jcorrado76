#include "EulerSolver.h"

double compute_rhov( float rho , float v ){
    return rho * v ;
}
double compute_E( float rho , float v , float p , float gamma ){
    return( p / ( gamma - 1. ) + 0.5 * rho * v * v );
}

void EulerSolver::prim_to_cons( const std::vector<float> &rho , const std::vector<float> &v , const std::vector<float> &p , std::vector<float>& rhov , std::vector<float>& E ){
    for ( int i = 0 ; i < grid_size ; i++ ){
        rhov[i] = compute_rhov( rho[i] , v[i] );
        E[i] = compute_E( rho[i] , v[i] , p[i] , gamma );
    }
}
