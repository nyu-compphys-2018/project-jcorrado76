#include "EulerSolver.h"

void EulerSolver::LU( const std::vector<float> rho , const std::vector<float> rhov , const std::vector<float> E , 
    std::vector<float> &rho_update, std::vector<float> &rhov_update, std::vector<float> &energy_density_update){

    std::vector<float> alphaP(grid_size);
    std::vector<float> alphaM(grid_size);

    std::vector<float> rho_Ls , rhov_Ls , E_Ls; // left conserved variables
    std::vector<float> rho_Rs , rhov_Rs , E_Rs; // right conserved variables
    std::vector<float> v_Ls , p_Ls, cs_Ls; // left primitives
    std::vector<float> v_Rs , p_Rs, cs_Rs; // right primitives

    // Reconstruct left and right states 
    Reconstruct_States(  &rho ,  &rhov ,  &energy , 
                         &rho_Ls ,  &rhov_Ls ,  &E_Ls
                         &rho_Rs ,  &rhov_Rs ,  &E_Rs );


    // convert left states to primitive 
    cons_to_prim( &rho_Ls , &rhov_Ls , &E_Ls , &v_Ls , &p_Ls );
    cons_to_prim( &rho_Rs , &rhov_Rs , &E_Rs , &v_Rs , &p_Rs );


    // get sound speeds 
    get_sound_speed( &rho_Ls , &p_Ls , &cs_Ls );
    get_sound_speed( &rho_Rs , &p_Rs , &cs_Rs );

    // compute alphas
    //
    //
    //
    //
    //
    //
    // compute physical fluxes as function of primitives
    //
    //
    //
    //
    // compute HLLE flux as a function of UL, UR, FL, FR , ap, am
    //

}
