#include <iostream>
#include <fstream>
#include "EulerSolver.h"

void EulerSolver::write_to_file( std::string fname ){
    std::cout << "Writing primitives to file: " << fname << std::endl;
    std::ofstream outFile;
    outFile.open("outputFile.csv");
    for ( int i = ilo ; i < ihi+1 ; i++ ) {
        outFile << x[i] << "," << rho[i] << "," << v[i] << "," << p[i] << "\n";
    }
    outFile.close();
}
