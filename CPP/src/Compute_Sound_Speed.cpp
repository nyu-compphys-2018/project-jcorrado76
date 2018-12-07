#include <cmath>
#include "EulerSolver.h"

float* EulerSolver::get_sound_speed( float* rho , float* p ){
    float cslocal[grid_size];
    for ( int i = 0 ; i < grid_size ; i++ ){
        cslocal[i] = sqrt( gamma * p[i] / rho[i] );
    }
    return cslocal;
}
