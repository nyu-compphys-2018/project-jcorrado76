#include "EulerSolver.h"



void EulerSolver::Reconstruct_States( const std::vector<float> rhos , const std::vector<float> rhovs , const std::vector<float> Es , 
                                      std::vector<float> &rho_Ls , std::vector<float> &rhov_Ls , std::vector<float> &E_Ls
                                      std::vector<float> &rho_Rs , std::vector<float> &rhov_Rs , std::vector<float> &E_Rs,   float theta ){
    if ( spatial_order == 1 ){
        // piecewise-constant left and right states 
        for ( int i = 0 ; i < Nx + 1 ; i++ ){
            rho_Ls[i] = rhos[ ilo + i - 1 ];
            rhov_Ls[i] = rhovs[ ilo + i - 1 ];
            E_Ls[i] = Es[ ilo + i - 1 ];

            rho_Rs[i] = rhos[ ilo + i ];
            rhov_Rs[i] = rhovs[ ilo + i ];
            E_Rs[i] = Es[ ilo + i ];
        }else{
        }
}
