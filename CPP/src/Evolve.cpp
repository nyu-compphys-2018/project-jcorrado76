#include "EulerSolver.h"

void EulerSolver::evolve( float tend ){
    tfinal = tend;
    float t = 0.0;
    while(t < tfinal){
        std::cout << "Time: " << t << std::endl;
        get_sound_speed( rho , p , cs);
        float dt = get_dt();
        if ( t + dt > tfinal ){ // avoid overshooting 
            dt = tfinal - t;
        }
        if ( time_order == 1){
            Forward_Euler_Update( dt );
        }else{
            RK3_Update( dt );
        }
        cons_to_prim( rho , rhov , energy , v , p );
        fill_BCs();
        t += dt;
    }
    write_to_file();
}
