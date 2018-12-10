#include "EulerSolver.h"

void EulerSolver::setSod( float x0, float rho_l , float v_l , float p_l, float rho_r , float v_r , float p_r , float gam ){
    std::cout << "Setting sod initial conditions" << std::endl;
    gamma = gam;
    std::cout << "gamma: " << gamma << std::endl;
    std::cout << "grid size: " << grid_size << std::endl;
    std::cout << "x0: " << x0<< std::endl;
    for ( int i = 0 ; i < grid_size ; i++ ){
        if ( x[i] <= x0 ){
            rho[i] = rho_l;
            v[i] = v_l;
            p[i] = p_l;

        }else{
            rho[i] = rho_r;
            v[i] = v_r;
            p[i] = p_r;
        }
    }
    prim_to_cons( rho , v , p , rhov , energy );

    std::string filename = "initial_conditions.csv";
    write_to_file( filename );

}
