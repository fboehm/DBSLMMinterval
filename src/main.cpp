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
#include <armadillo>
#include <sys/stat.h>
#include <sys/types.h>

#include "../include/dbslmm.hpp"
#include "../include/read_control.h"

using namespace std;
using namespace arma;

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
  cDB.Assign(argc, argv, cPar);
  //read param settings file here!
  std::vector<std::vector<std::string>> control = readControlFile(cPar.control_file);
  cout << "control length is: " << control.size() << endl;
  //initialize a arma::field to store outputs for var calcs!
  arma::field < arma::field < arma::mat> > results; 
  for (int i = 0; i < control.size(); ++i){
    std::vector < std::string> rr = control[i];
    //create cPar object here:
    cPar.b = rr[1];// chr is column 0; then do the other args in alphabetical order
    cPar.eff = rr[2];
    cPar.h = std::stod(rr[3]);
    cPar.l = rr[4];
    cPar.mafMax = std::stod(rr[5]);//https://www.programiz.com/cpp-programming/string-float-conversion
    cPar.n = std::stoi(rr[6]);
    cPar.nsnp = std::stoi(rr[7]);
    cPar.r = rr[8];
    cPar.s = rr[9];
    cPar.t = std::stoi(rr[10]);
    
    // call BatchRun
    results(i) = cDB.BatchRun(cPar);
    
  }
  //var calcs here! use contents of results field of field of matrices
  //1. assemble genome-wide matrices from "results"
  // results is a 22-long field where each entry is itself a 2d field containing matrices
  
  //2. input matrices to calc_asymptotic_variance
  
  //3. write diagonal of var to a csv file
  
  
  
  return EXIT_SUCCESS;
}
