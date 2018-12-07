#include <algorithm>
#include "EulerSolver.h"

float EulerSolver::get_dt(){
    float maxP;
    float maxM;
    float lambdaPlus[grid_size];
    float lambdaMinus[grid_size];

    for ( int  i = 0 ; i < grid_size ; i++ ) {
        lambdaPlus[i] = lambdaP( v[i] , cs[i] );
        lambdaMinus[i] = lambdaM( v[i] , cs[i] );
    }
    maxP = *std::max_element( std::begin( lambdaPlus ) , std::end( lambdaPlus ) );
    maxM = *std::max_element( std::begin( lambdaMinus ) , std::end( lambdaMinus ) );
    float maxEigVal = std::max( maxP , maxM );

    float dt = cfl * dx / maxEigVal;
    return dt;
}

