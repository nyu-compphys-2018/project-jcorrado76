#include <cmath>
#include "EulerSolver.h"

void EulerSolver::get_sound_speed( const std::vector<float> &rho , const std::vector<float> &p , std::vector<float> &csLocal ){
    for ( int i = 0 ; i < csLocal.size() ; i++ ){
        cslocal[i] = (sqrt( gamma * p[i] / rho[i] ));
    }
}
