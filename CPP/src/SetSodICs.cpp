#include "EulerSolver.h"

void EulerSolver::setSod( float x0, float rho_l , float v_l , float p_l, float rho_r , float v_r , float p_r , float gam ){
    std::cout << "Setting sod initial conditions" << std::endl;
    gamma = gam;
    std::cout << "gamma: " << gamma << std::endl;
    std::cout << "grid size: " << grid_size << std::endl;
    std::cout << "x0: " << x0<< std::endl;
    for ( int i = 0 ; i < grid_size ; i++ ){
        rho[i] = ( (x[i] <= x0) ? rho_l : rho_r );
        v[i] = ( (x[i] <= x0) ? v_l : v_r );
        p[i] = ( (x[i] <= x0) ? p_l : p_r );
    }
    prim_to_cons( rho , v , p , rhov , energy );

    std::string filename = "initial_conditions.csv";
    write_to_file( filename );

}
