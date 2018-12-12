#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fenv.h> // this include to have functions to stop when nan is encountered

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
        void setSod( float x0=0.5 , float rho_L=1. , float v_l=0. , float p_l=1. , float rho_r=0.125 , float v_r=0.0 , float p_r=0.1 , float gam=1.4 );
        void get_sound_speed( const std::vector<float> &rho , const std::vector<float> &p , std::vector<float> &cslocal);
        void Forward_Euler_Update( float dt );
        void RK3_Update( float dt );
        void cons_to_prim( const std::vector<float> &rho , const std::vector<float> &rhov , const std::vector<float> &E ,std::vector<float> &vs , std::vector<float> &ps);
        void prim_to_cons( const std::vector<float> &rho , const std::vector<float> &v , const std::vector<float> &p, std::vector<float>& rhov , std::vector<float>& E );
        void fill_BCs();
        void evolve( float tend );
        void Physical_Flux( const std::vector<float> rho , const std::vector<float> v , const std::vector<float> p, 
                std::vector<float> &rho_flux , std::vector<float> &rhov_flux , std::vector<float> &total_energy_density_flux);

        void HLLE_Flux( const std::vector<float>& rho_l , const std::vector<float>& rhov_l , const std::vector<float>& E_l , 
                             const std::vector<float>& rho_r , const std::vector<float>& rhov_r , const std::vector<float>& E_r ,
                             const std::vector<float>& f1_l , const std::vector<float>& f2_l , const std::vector<float>& f3_l, 
                             const std::vector<float>& f1_r , const std::vector<float>& f2_r , const std::vector<float>& f3_r, 
                             const std::vector<float>& am , const std::vector<float>& ap,
                             std::vector<float>& hlle1 , std::vector<float>& hlle2 , std::vector<float>& hlle3);
        float get_dt();
        void alphaP( const std::vector<float> &rho_l , const std::vector<float> &v_l , const std::vector<float> &p_l , 
                            const std::vector<float> &rho_r , const std::vector<float> &v_r , const std::vector<float> &p_r  , 
                            const std::vector<float> &csL , const std::vector<float> &csR , 
                            std::vector<float> &ap );
        void alphaM( const std::vector<float> &rho_l , const std::vector<float> &v_l , const std::vector<float> &p_l , 
                            const std::vector<float> &rho_r , const std::vector<float> &v_r , const std::vector<float> &p_r  , 
                            const std::vector<float> &csL , const std::vector<float> &csR , 
                            std::vector<float> &am);
        
        void LU( const std::vector<float> rho , const std::vector<float> rhov , const std::vector<float> E , 
                 std::vector<float> &rho_update, std::vector<float> &rhov_update, std::vector<float> &energy_density_update);

        void Reconstruct_States(const std::vector<float> &rhos , const std::vector<float> &rhovs , const std::vector<float> &Es , 
                                      std::vector<float> &rho_Ls , std::vector<float> &rhov_Ls , std::vector<float> &E_Ls ,
                                      std::vector<float> &rho_Rs , std::vector<float> &rhov_Rs , std::vector<float> &E_Rs,   const float theta=1.5 );
        void write_to_file( std::string fname="outputFile.csv" );
};
