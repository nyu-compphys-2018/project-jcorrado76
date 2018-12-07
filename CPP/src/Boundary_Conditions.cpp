#include "EulerSolver.h"
void EulerSolver::fill_BCs(){
    if ( bc == "periodic" ){
        for ( int  i = 0 ; i < ilo ; i++ ){
            rho[i] = rho[ ihi - Ng + i + 1];
            rhov[i] = rhov[ ihi - Ng + i + 1];
            energy[i] = energy[ ihi - Ng + i + 1];
            rho[ ihi + i + 1 ] = rho[ ilo + i ];
            rhov[ ihi + i + 1 ] = rhov[ ilo + i ];
            energy[ ihi + i + 1 ] = energy[ ilo + i ];
        }
    }
    else if ( bc == "outflow" ){
        for ( int i = 0 ; i < ilo ; i++ ){
            rho[ i ] = rho[ ilo ];
            rhov[ i ] = rhov[ ilo ];
            energy[ i ] = energy[ ilo ];

            rho[ ihi + 1 + i ] = rho[ ihi ];
            rhov[ ihi + 1 + i ] = rhov[ ihi ];
            energy[ ihi + 1 + i ] = energy[ ihi ];
        }
    }
}
