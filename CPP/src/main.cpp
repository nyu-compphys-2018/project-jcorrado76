#include "EulerSolver.h"
#include <iostream>



int main(){
    float tfinal = 0.1;
    int Nx = 400;
    float xmin=0.0; float xmax=1.0;
    float cfl = 0.5;
    int spatial_order=1; int time_order=1;
    std::string bc="outflow";
    EulerSolver e( Nx , xmin, xmax,cfl, spatial_order, time_order, bc);
    e.setSod();
    e.evolve( tfinal );
    return 1;
}
