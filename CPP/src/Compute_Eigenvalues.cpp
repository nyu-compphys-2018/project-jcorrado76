#include "EulerSolver.h"
float EulerSolver::lambdaP( float v , float cs ){
    return(v + cs);
}
float EulerSolver::lambdaM( float v , float cs ){
    return( v - cs );
}
