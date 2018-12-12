#include "EulerSolver.h"

void EulerSolver::LU( const std::vector<float> rho , const std::vector<float> rhov , const std::vector<float> E , 
    std::vector<float> &rho_update, std::vector<float> &rhov_update, std::vector<float> &energy_density_update){

    std::vector<float> ap(Nx+1);
    std::vector<float> am(Nx+1);

    // left conserved variables
    std::vector<float> rho_Ls(Nx+1);
    std::vector<float> rhov_Ls(Nx+1);
    std::vector<float> E_Ls(Nx+1);     
    // right conserved variables
    std::vector<float> rho_Rs(Nx+1);
    std::vector<float> rhov_Rs(Nx+1);
    std::vector<float> E_Rs(Nx+1); 
    //left primitives 
    std::vector<float> v_Ls(Nx+1);
    std::vector<float> p_Ls(Nx+1);
    std::vector<float> cs_Ls(Nx+1); 
    // right primitives 
    std::vector<float> v_Rs(Nx+1);
    std::vector<float> p_Rs(Nx+1);
    std::vector<float> cs_Rs(Nx+1); 
    // left fluxes
    std::vector<float> fl_1(Nx+1);
    std::vector<float> fl_2(Nx+1);
    std::vector<float> fl_3(Nx+1);
    // right fluxes
    std::vector<float> fr_1(Nx+1);
    std::vector<float> fr_2(Nx+1);
    std::vector<float> fr_3(Nx+1);

    for ( int i = 0 ; i < grid_size ; i++ ){
        if ( std::isnan(rho[i]) || std::isnan( rhov[i]) || std::isnan( E[i] ) ){
            std::cout << "Detected nan in time averaged states" << std::endl;
            return;
        }
    }

    // Reconstruct left and right states 
    Reconstruct_States(  rho ,  rhov ,  energy , 
                         rho_Ls ,  rhov_Ls ,  E_Ls ,
                         rho_Rs ,  rhov_Rs ,  E_Rs );


    for ( int i = 0 ; i < Nx+1 ; i++ ){
        if ( std::isnan(rho_Ls[i]) || std::isnan( rho_Rs[i]) || std::isnan( rhov_Ls[i] ) ){
            std::cout << "Detected nan in reconstructed left and right states" << std::endl;
            return;
        }
    }

    // convert left states to primitive 
    cons_to_prim( rho_Ls , rhov_Ls , E_Ls , v_Ls , p_Ls );
    cons_to_prim( rho_Rs , rhov_Rs , E_Rs , v_Rs , p_Rs );


    // get sound speeds 
    get_sound_speed( rho_Ls , p_Ls , cs_Ls );
    get_sound_speed( rho_Rs , p_Rs , cs_Rs );

    alphaP( rho_Ls ,   v_Ls ,   p_Ls , 
             rho_Rs ,   v_Rs ,   p_Rs  , 
             cs_Ls ,   cs_Rs , 
             ap );

    alphaM( rho_Ls ,   v_Ls ,   p_Ls , 
             rho_Rs ,   v_Rs ,   p_Rs  , 
             cs_Ls ,   cs_Rs , 
             am );

    // compute left physical fluxes
    Physical_Flux( rho_Ls , v_Ls , p_Ls, fl_1 , fl_2 , fl_3 );
    // compute right physical fluxes
    Physical_Flux( rho_Rs , v_Rs , p_Rs, fr_1 , fr_2 , fr_3 );

    HLLE_Flux(  rho_Ls ,  rhov_Ls ,  E_Ls , 
                rho_Rs ,  rhov_Rs ,  E_Rs ,
                fl_1 ,  fl_2,  fl_3, 
                fr_1 ,  fr_2 ,  fr_3, 
                am ,  ap,
                rho_update , rhov_update, energy_density_update);

    // checking for nans in the update
    for ( int i = 0 ; i < Nx+1 ; i++ ){
        if ( std::isnan(rho_update[i]) || std::isnan( rhov_update[i]) || std::isnan( energy_density_update[i] ) ){
            std::cout << "Detected nan in update" << std::endl;
            return;
        }
    }
}
