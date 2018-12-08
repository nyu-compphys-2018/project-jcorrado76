#include "EulerSolver.h"
void EulerSolver::setSod( float x0, float rho_l , float v_l , float p_l, float rho_r , float v_r , float p_r , float gam ){
    std::cout << "Setting sod initial conditions" << std::endl;
    gamma = gam;
    std::cout << "gamma: " << gamma << std::endl;
    std::cout << "grid size: " << grid_size << std::endl;
    std::cout << "x0: " << x0<< std::endl;
    for ( int i = 0 ; i < grid_size ; i++ ){
        std::cout << "i: " << i << std::endl;
        std::cout << "x[i] " << x[i] << std::endl;
        std::cout << "x0: " << x0 << std::endl;
        if ( x[i] <= x0 ){
            std::cout << "Filling left hand states" << std::endl;
            rho[i] = rho_l;
            v[i] = v_l;
            p[i] = p_l;

        }else{
            std::cout << "Filling right hand states" << std::endl;
            rho[i] = rho_r;
            v[i] = v_r;
            p[i] = p_r;
        }
    }
    //Now you need to call prim to cons to initialize conserved variables 
}
