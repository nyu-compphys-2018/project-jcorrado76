#include "EulerSolver.h"

void EulerSolver::Forward_Euler_Update( float dt ){
    // HLL updates. same size as number physical cells
    std::vector<float> fhlle_1(Nx);
    std::vector<float> fhlle_2(Nx);
    std::vector<float> fhlle_3(Nx);
    LU( rho , rhov , energy , fhlle_1 , fhlle_2 , fhlle_3 );
    for ( int i = 0; i < Nx ; i++ ){
        rho[ ilo + i ] += fhlle_1[i] * dt;
        rhov[ ilo + i ] += fhlle_2[i] * dt;
        energy[ ilo + i ] += fhlle_3[i] * dt;
    }
}
