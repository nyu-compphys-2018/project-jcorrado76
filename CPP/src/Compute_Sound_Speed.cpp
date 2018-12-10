#include <cmath>
#include "EulerSolver.h"

void EulerSolver::get_sound_speed( const std::vector<float> &rho , const std::vector<float> &p , std::vector<float> &cslocal ){
    for ( int i = 0 ; i < cslocal.size() ; i++ ){
        cslocal[i] = (sqrt( gamma * p[i] / rho[i] ));
    }
}
