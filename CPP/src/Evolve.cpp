

void EulerSolver::evolve( float tend ){
    tfinal = tend;
    while t < tfinal{
        cs = get_sound_speed( rho , p );
        dt = get_dt();
        if ( t + dt > tfinal ){ // avoid overshooting 
            dt = tfinal - t;
        }
        if ( time_order == 1){
            Forward_Euler_Update();
        }else{
            RK3_Update();
        }
        cons_to_prim();
        fill_BCs();
        t += dt;
    }
}
