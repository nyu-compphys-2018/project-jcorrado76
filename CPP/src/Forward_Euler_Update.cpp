#include "EulerSolver.h"

void EulerSolver::Forward_Euler_Update( float dt ){
    udot = LU();
    self.U[:,sself.physical] += dt * udot;
}
