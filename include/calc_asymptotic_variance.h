#define ARMA_DONT_USE_WRAPPER

#include <iostream>
#include <vector>
#include <armadillo>
#include <string>
#include <boost/algorithm/string.hpp> // split
#include <tuple> // std::tuple, std::get, std::tie, std::ignore
#include <algorithm> 


arma::mat calc_asymptotic_variance(arma::mat Sigma_ll, 
                                   arma::mat Sigma_ls, 
                                   arma::mat Sigma_ss, 
                                   double sigma2_s, 
                                   unsigned int n,
                                   arma::mat Xl_test, 
                                   arma::mat Xs_test);

arma::mat calc_var_betal(arma::mat Sigma_ll, 
                         arma::mat Sigma_ls, 
                         arma::mat Sigma_ss,
                         arma::mat A_inverse,
                         unsigned int n);

  arma::mat calc_var_betas(arma::mat Sigma_ss, 
                           arma::mat Sigma_ls,
                           arma::mat A_inverse,
                           double sigma2_s,
                           unsigned int n,
                           arma::mat var_bl);

arma::mat calc_A_inverse(arma::mat Sigma_ss, 
                         double sigma2_s, 
                         unsigned int n);


arma::mat BlockDiag( arma::field<mat> x );
