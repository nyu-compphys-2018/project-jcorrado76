#include <string>

class EulerSolver{
    private:
        int Nx; // number of physical cells 
        float xmin , xmax; // x values for physical boundaries 
        float cfl; // CFL number 
        float t; // current time 
        int spatial_order, time_order; // order in time and space
        std::string bc; // boundary condition
        if ( spatial_order == 1 ){
            int Ng = 1;
        }else{
            int Ng = 2;
        }
        int ilo = Ng; // first physical index
        int ihi = Ng + Nx -1; // last physical index 
        int grid_size = Nx + 2 * Ng; // total size of numerical grid, including guard cells
        float dx = (xmax - xmin) / Nx; // physical width of cells 
        float x [grid_size]; // array to store x values 
        // fill x values 
        for ( int i = 0 ; i < grid_size ; i++ ){
            x[i] = xmin + ( i - Ng + 0.5 ) * dx;
        }
        // conserved variables
        float rho[grid_size]; // mass density 
        float rhov[grid_size]; // momentum density 
        float energy[grid_size]; // total energy density
        // primitives 
        float v[grid_size]; // velocity 
        float p[grid_size]; // pressure 
        float cs[grid_size]; // speed of sound 
        float gamma = 1.4; // adiabatic index 
        float tfinal; // tfinal

    public:
        EulerSolver( int Nx=10 , float xmin=0.0 , float xmax=1.0 , float cfl=0.5 , int spatial_order=1 , int time_order=1 , std::string bc="outflow"):
            Nx(Nx) , xmin(xmin) , xmax(xmax) , cfl(cfl) , spatial_order(spatial_order) , time_order(time_order) , bc(bc) {}
        float lambdaP( float v , float cs );
        float lambdaM( float v , float cs );
        void setSod( float x0=0.5 , float rho_=1. , float v_l=0. , float p_l=1. , float rho_r=0.125 , float v_r=0.0 , float p_r=0.1 , float gam=1.4 );
        float get_sound_speed( float rho , float p );
        void Forward_Euler_Update( float dt );
        void RK3_Update( float dt );
        void cons_to_prim( float rho , float rhov , float E );
        void prim_to_cons( float rho , float v , float p );
        void fill_BCs();
        void evolve( float tend );
        void Physical_Flux( float rho , float v , float p );
        void HLLE_Flux( float rho_l , float rhov_l , float E_l , float rho_r , float rhov_r , float E_r ,
                        float f1_l , float f2_l , float f3_l,float f1_r , float f2_r , float f3_r, float am , float ap );
        float get_dt();
        void Reconstruct_States( float rho , float rhov , float E , float theta=1.5 );
        double alphaP( float rho_l , float v_l , float p_l , float rho_r , float v_r , float p_r  , float csL , float csR );
        double alphaM( float rho_l , float v_l , float p_l , float rho_r , float v_r , float p_r  , float csL , float csR );
        void LU( float rho , float rhov , float E );
};
