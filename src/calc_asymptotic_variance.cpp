#define ARMA_DONT_USE_WRAPPER

#include <iostream>
#include <vector>
#include <armadillo>
#include <string>
#include <boost/algorithm/string.hpp> // split
#include <algorithm>

#include "../include/calc_asymptotic_variance.h"

using namespace std;
using namespace arma;

//' Calculate the asymptotic variance for the predicted y values
//' 
//' @param Sigma_ll Sigma_ll matrix for the whole genome
//' @param Sigma_ls Sigma_ls matrix for the whole genome
//' @param Sigma_ss Sigma_ss matrix for the whole genome
//' @param sigma2_s estimated value of sigma^2_s
//' @param n sample size
//' @param Xl_test genotypes matrix for large effect SNPs for test subjects
//' @param Xs_test genotypes matrix for small effect SNPs for test subjects
//' @return variance of predicted y values

arma::mat calc_asymptotic_variance(arma::mat Sigma_ll, 
                                   arma::mat Sigma_ls, 
                                   arma::mat Sigma_ss, 
                                   double sigma2_s, 
                                   unsigned int n,
                                   arma::mat Xl_test, 
                                   arma::mat Xs_test){
  arma::mat Ainv = calc_A_inverse(Sigma_ss, sigma2_s, n);
  arma::mat var_bl = calc_var_betal(Sigma_ll, 
                                    Sigma_ls, 
                                    Sigma_ss, 
                                    Ainv, 
                                    n);
  arma::mat var_bs = calc_var_betas(Sigma_ss, 
                                    Sigma_ls,
                                    Ainv,
                                    sigma2_s,
                                    n,
                                    var_bl);
  arma::mat result = arma::trans(Xl_test) * var_bl * Xl_test + arma::trans(Xs_test) * var_bs * Xs_test;
  return(result);
}

//' Calculate A inverse matrix
//' 
//' @details (sigma^{-2}n^{-1} I_ms + Sigma_ss) = A.
//' @param Sigma_ss Sigma_ss matrix for a single LD block
//' @param sigma2_s estimate of sigma^2_s
//' @param n sample size
//' @return A inverse matrix

arma::mat calc_A_inverse(arma::mat Sigma_ss, 
                         double sigma2_s, 
                         unsigned int n){
  //determine m_s
  unsigned int m_s = Sigma_ss.n_rows;
  // calculate result
  arma::mat result = arma::inv(arma::eye(m_s, m_s) / (n * sigma2_s) + Sigma_ss);
  return result;
}



//' Calculate variance of coefficient estimator for large effects
//' 
//' @param Sigma_ll Sigma_ll constructed for one LD block
//' @param Sigma_ls Sigma_ls constructed for one LD block
//' @param Sigma_ss Sigma_ss constructed for one LD block
//' @param A_inverse inverse of (sigma^{-2}n^{-1} I_ms + Sigma_ss)
//' @param n sample size
//' @return covariance matrix

arma::mat calc_var_betal(arma::mat Sigma_ll, 
                      arma::mat Sigma_ls, 
                      arma::mat Sigma_ss,
                      arma::mat A_inverse,
                      unsigned int n){
  //calculate second matrix
  arma::mat big = Sigma_ll - Sigma_ls * A_inverse * arma::trans(Sigma_ls);
  //invert and divide by n
  arma::mat result = arma::inv(big) / n;
  return result;
}

//' Calculate variance of coefficient estimator for small effects
//' 
//' @param Sigma_ss Sigma_ss matrix 
//' @param Sigma_ls Sigma_ls matrix 
//' @param A_inverse A inverse matrix 
//' @param sigma2_s estimated value of sigma^2_s
//' @param n sample size
//' @param var_bl variance of beta hat l
//' @return covariance matrix
  
arma::mat calc_var_betas(arma::mat Sigma_ss, 
                         arma::mat Sigma_ls,
                         arma::mat A_inverse,
                         double sigma2_s,
                         unsigned int n,
                         arma::mat var_bl){
  arma::mat small = arma::trans(Sigma_ls) - Sigma_ss * A_inverse * arma::trans(Sigma_ls);
  arma::mat term2 = small * var_bl * arma::trans(small);
  arma::mat term1 = Sigma_ss - Sigma_ss * A_inverse * Sigma_ss;
  arma::mat result = sigma2_s * sigma2_s * n * (term1 + term2);
  return result;
}


//' Construct a block diagonal matrix from a collection of matrices
//' 
//' @param x a field of matrices, possibly of different sizes
//' @return a block diagonal matrix
//' @reference https://stackoverflow.com/questions/29198893/block-diagonal-matrix-armadillo

arma::mat BlockDiag( arma::field<arma::mat> x ) {
  
  unsigned int len = x.n_rows;
  cout << "number of rows in field is " << len << endl;
  int drow = 0;
  int dcol = 0;
  arma::ivec rvec(len);
  arma::ivec cvec(len);
  //get dimensions of each matrix in the field
  for(unsigned int i=0; i< len; i++) {
    rvec(i) = x(i, 0).n_rows ; 
    drow += rvec(i) ;
    cvec(i) = x(i, 0).n_cols ; 
    dcol += cvec(i);
  }
  cout << "number of rows in big block diagonal is " << drow << endl;
  cout << "number of columns in big block diagonal is " << dcol << endl;
  
  //initialize matrix to be returned
  arma::mat X(drow, dcol, fill::zeros);
  int idx_row = 0;
  int idx_col = 0;
  // place matrices at correct places
  for(unsigned int i=0; i< len; i++) {
    X.submat(idx_row, 
             idx_col, 
             idx_row + rvec(i) - 1, 
             idx_col + cvec(i) - 1) = x(i, 0) ;
    idx_row = idx_row + rvec(i) ;
    idx_col = idx_col + cvec(i);
  }
  return(X);
}


