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
  int nchr = control.size() / 2; // nchr is number of chromosomes in genome
  //control file always has a number of lines equal to twice the number of chromosomes
  //initialize a arma::field to store outputs for var calcs!
  arma::field <arma::mat > ff;
  arma::field < arma::mat> training(nchr, 5);
  arma::field < arma::mat> test(nchr, 5);
  for (int i = 0; i < 2 * nchr; ++i){
    std::vector < std::string> rr = control[i];
    //create cPar object here: COORDINATE WITH CONTROL FILE STRUCTURE
    // 
    cPar.b = rr[1];// chr is column 0; then do the other args in alphabetical order
    // b is the path to the block data files directory
    cPar.eff = rr[2]; //file path for outputting the snp effect estimates
    cPar.h = std::stod(rr[3]); //heritability
    cPar.l = rr[4]; //large effects summary file path
    cPar.mafMax = std::stod(rr[5]);//https://www.programiz.com/cpp-programming/string-float-conversion
    cPar.n = std::stoi(rr[6]); //sample size
    cPar.nsnp = std::stoi(rr[7]); //number of snps (genomewide)
    cPar.r = rr[8]; // bfile for reference data 
    cPar.s = rr[9]; //small effects summary file path
    cPar.t = std::stoi(rr[10]); //number of threads
    double sigma2_s = cPar.h / (double)cPar.nsnp;
    // call BatchRun
    ff = cDB.BatchRun(cPar);
    if (i < nchr){
      training.row(i) = assembleMatrices(ff);
    } else {
      test.row(i - nchr) = assembleMatrices(ff);
    }
  }
  //var calcs here! 
  //1. assemble genome-wide matrices from "results"
  // results is a 22-long field where each entry is itself a 1d field containing 5 matrices
  arma::field < arma::mat > mats_training = assembleMatrices(training);
  arma::field < arma::mat > mats_test = assembleMatrices(test);
    //2. input matrices to calc_asymptotic_variance
  arma::mat vv = calc_asymptotic_variance(mats_training(2), t(mats_training(1)), mats_training(0), 
                                          cPar.h, cPar.n, mats_test(4), mats_test(3));
  //3. write diagonal of var to a csv file
  arma::vec vd = diagvec(vv);
  vd.save("out.csv", csv_ascii);
  
  
  return EXIT_SUCCESS;
}
