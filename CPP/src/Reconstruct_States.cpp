#include "EulerSolver.h"



void EulerSolver::Reconstruct_States( const std::vector<float> &rhos , const std::vector<float> &rhovs , const std::vector<float> &Es , 
                                      std::vector<float> &rho_Ls , std::vector<float> &rhov_Ls , std::vector<float> &E_Ls,
                                      std::vector<float> &rho_Rs , std::vector<float> &rhov_Rs , std::vector<float> &E_Rs,   const float theta ){
    if ( spatial_order == 1 ){
        // piecewise-constant left and right states 
        for ( int i = 0 ; i < Nx + 1 ; i++ ){
            rho_Ls[i] = rhos[ ilo + i - 1 ];
            rhov_Ls[i] = rhovs[ ilo + i - 1 ];
            E_Ls[i] = Es[ ilo + i - 1 ];

            rho_Rs[i] = rhos[ ilo + i ];
            rhov_Rs[i] = rhovs[ ilo + i ];
            E_Rs[i] = Es[ ilo + i ];

            if ( std::isnan(rhos[ilo+i-1]) ){
                std::cout << "nan in rhoLs at " << ilo+i-1 << ": " << rho_Ls[ilo+i-1] << std::endl;
            }
            if ( std::isnan(rhovs[ilo+i-1]) ){
                std::cout << "nan in rhovLs at " << ilo+i-1 << ": " << rhov_Ls[ilo+i-1] << std::endl;
            }
            if ( std::isnan(Es[ilo+i-1]) ){
                std::cout << "nan in ELs at " << ilo+i-1 << ": " << E_Ls[ilo+i-1] << std::endl;
            }
            if ( std::isnan(rhos[ilo+i]) ){
                std::cout << "nan in rhoRs at " << ilo+i << ": " << rho_Rs[ilo+i-1] << std::endl;
            }
            if ( std::isnan(rhovs[ilo+i]) ){
                std::cout << "nan in rhovRs at " << ilo+i << ": " << rhov_Rs[ilo+i-1] << std::endl;
            }
            if ( std::isnan(Es[ilo+i]) ){
                std::cout << "nan in ERs at " << ilo+i  << ": " << E_Rs[ilo+i-1] << std::endl;
            }
            
        }
    }else{}
}
