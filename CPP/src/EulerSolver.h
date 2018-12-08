#include <string>
#include <iostream>
#include <vector>

class EulerSolver{
    private:
        int Nx; // number of physical cells 
        float xmin; // x value for lower bound 
        float xmax; // x values for upper bound 
        float cfl; // CFL number 
        float t; // current time 
        int spatial_order; // order in space
        int time_order; // order in time
        std::string bc; // boundary condition
        int Ng; // number of guard cells 
        int ilo,ihi; // first,last physical index
        int grid_size; // total size of numerical grid, including guard cells
        float dx; // physical width of cells 
        std::vector<float> x; // array to store x values 
        // conserved variables
        std::vector<float> rho;
        std::vector<float> rhov;
        std::vector<float> energy;
        // primitives 
        std::vector<float> v;
        std::vector<float> p;
        std::vector<float> cs;
        float gamma; // adiabatic index 
        float tfinal; // tfinal

    public:
        EulerSolver( int Nx=10 , float xmin=0.0 , float xmax=1.0 , float cfl=0.5 , int spatial_order=1 , int time_order=1 , std::string bc="outflow")
            : Nx(Nx) , xmin(xmin) , xmax(xmax) , cfl(cfl) , spatial_order(spatial_order) , time_order(time_order) , bc(bc) , grid_size(Nx + 2 * Ng){
            std::cout << spatial_order << std::endl;
            if ( spatial_order == 1 ){
                Ng = 1;
            }else{
                Ng = 2;
            }
            ilo = Ng; // first physical index
            ihi = Ng + Nx -1; // last physical index 
            grid_size = Nx + 2 * Ng; // total size of numerical grid, including guard cells
            dx = (xmax - xmin) / Nx; // physical width of cells 
            for ( int i = 0 ; i < grid_size ; i++ ){ // fill x values 
                x.push_back(xmin + ( i - Ng + 0.5 ) * dx);
            }
            x.resize(grid_size); 
            rho.resize(grid_size);
            rhov.resize(grid_size);
            energy.resize(grid_size);
            v.resize(grid_size);
            p.resize(grid_size);
            cs.resize(grid_size);
        };
        float lambdaP( float v , float cs );
        float lambdaM( float v , float cs );
        void setSod( float x0=0.5 , float rho_=1. , float v_l=0. , float p_l=1. , float rho_r=0.125 , float v_r=0.0 , float p_r=0.1 , float gam=1.4 );
        std::vector<float> get_sound_speed( std::vector<float> rho , std::vector<float> p );
        void Forward_Euler_Update( float dt );
        void RK3_Update( float dt );
        void cons_to_prim( std::vector<float> rho , std::vector<float> rhov , std::vector<float> E );
        void prim_to_cons( float rho , float v , float p, float& u_rho , float& rhov , float& E );
        void fill_BCs();
        void evolve( float tend );
        void Physical_Flux( float rho , float v , float p );
        void HLLE_Flux( float rho_l , float rhov_l , float E_l , float rho_r , float rhov_r , float E_r ,
                        float f1_l , float f2_l , float f3_l,float f1_r , float f2_r , float f3_r, float am , float ap );
        float get_dt();
        double alphaP( float rho_l , float v_l , float p_l , float rho_r , float v_r , float p_r  , float csL , float csR );
        double alphaM( float rho_l , float v_l , float p_l , float rho_r , float v_r , float p_r  , float csL , float csR );
        void LU( float rho , float rhov , float E );
        void Reconstruct_States( const std::vector<float> rhos , const std::vector<float> rhovs , const std::vector<float> Es , 
                                      std::vector<float> &rho_Ls , std::vector<float> &rhov_Ls , std::vector<float> &E_Ls
                                      std::vector<float> &rho_Rs , std::vector<float> &rhov_Rs , std::vector<float> &E_Rs,   float theta=1.5 );
        void write_to_file( std::string fname="" );
};
