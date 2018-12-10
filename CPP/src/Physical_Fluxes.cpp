#include "EulerSolver.h"

float compute_rho_flux( float rho , float v ){
    return rho * v ;
}
float compute_rhov_flux( float rho , float v , float p ){
    return rho * v * v + p;
}
float compute_energy_flux( float v , float p , float E ){
    return (E + p)*v;
}
float compute_total_energy_density( float rho , float v , float p , float gamma ){
    return(( p / ( gamma - 1.0 ) ) + 0.5 * rho * v * v);
}

void EulerSolver::Physical_Flux( const std::vector<float> rho , const std::vector<float> v , const std::vector<float> p,
        std::vector<float> &rho_flux , std::vector<float> &rhov_flux , std::vector<float> &total_energy_density_flux){
    float E;
    for ( int i = 0 ; i < grid_size ; i++ ){
        rho_flux[i] = compute_rho_flux( rho[i] , v[i] );
        rhov_flux[i] = compute_rhov_flux( rho[i] , v[i] , p[i]);
        E = compute_total_energy_density( rho[i] , v[i] , p[i] , gamma );
        total_energy_density_flux[i] = compute_energy_flux( v[i], p[i], E );
    }
}

