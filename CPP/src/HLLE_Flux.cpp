#include "EulerSolver.h"

void EulerSolver::HLLE_Flux( const std::vector<float>& rho_l , const std::vector<float>& rhov_l , const std::vector<float>& E_l , 
                             const std::vector<float>& rho_r , const std::vector<float>& rhov_r , const std::vector<float>& E_r ,
                             const std::vector<float>& f1_l , const std::vector<float>& f2_l , const std::vector<float>& f3_l, 
                             const std::vector<float>& f1_r , const std::vector<float>& f2_r , const std::vector<float>& f3_r, 
                             const std::vector<float>& am , const std::vector<float>& ap,
                             std::vector<float>& hlle1 , std::vector<float>& hlle2 , std::vector<float>& hlle3){

    std::vector<float> fhll1(Nx+1);
    for ( int  i = 0 ; i < Nx+1 ; i++ ){
        fhll1 = ( ap[i] * f1_l[i] + am[i] * f1_2[i] - ap[i] * am[i] * ( rho_r[i] - rho_l[i] ) ) / ( ap[i]+am[i] );
        fhll2 = ( ap[i] * f2_l[i] + am[i] * f2_2[i] - ap[i] * am[i] * ( rhov_r[i] - rhov_l[i] ) ) / ( ap[i]+am[i] );
        fhll3 = ( ap[i] * f3_l[i] + am[i] * f3_2[i] - ap[i] * am[i] * ( E_r[i] - E_l[i] ) ) / ( ap[i]+am[i] );
    }

    for ( int i = 0 ; i < Nx ; i++ ){
        hlle1[i] = -(fhll1[i+1]-fhll1[i]) / dx;
        hlle2[i] = -(fhll2[i+2]-fhll2[i]) / dx;
        hlle3[i] = -(fhll3[i+3]-fhll3[i]) / dx;
    }
}
