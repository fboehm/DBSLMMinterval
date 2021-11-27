/*
 Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)
 Copyright (C) 2019  Sheng Yang and Xiang Zhou
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <armadillo>

#include "../include/dbslmm.hpp"

using namespace arma;
using namespace std;

int main(int argc, char * argv[])
{
  DBSLMM cDB;
  PARAM cPar;
  
  if (argc <= 1) {
    cDB.printHeader();
    return EXIT_SUCCESS;
  }
  if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h') {
    cDB.printHelp();
    return EXIT_SUCCESS;
  }
  /*for (int i = 0; i < argc; ++i)
    cout << argv[i] << "\n";
  */
  for(int i = 0; i < argc; ++i)
    cout << argv[i] << '\n';
  
  
  cDB.Assign(argc, argv, cPar);
    //initialize a arma::field to store outputs for var calcs!
    //double sigma2_s = cPar.h / (double)cPar.nsnp;
    // call BatchRun
  arma::field <arma::mat> ff = cDB.BatchRun(cPar);
  ff.save("out.bin", arma_binary);
  
  return EXIT_SUCCESS;
}
