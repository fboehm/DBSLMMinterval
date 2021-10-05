#define ARMA_DONT_USE_WRAPPER

#include <iostream>
#include <vector>
#include <armadillo>
#include <string>
#include <boost/algorithm/string.hpp> // split
#include <tuple> // std::tuple, std::get, std::tie, std::ignore
#include <algorithm>

#include "../include/calc_asymptotic_variance.h"

using namespace std;


//' Calculate the asymptotic variance for the predicted y values
//' 
//' @param Xl_training genotypes matrix for large effects markers in training cohort
//' @param Xs_training genotypes matrix for small effects markers in training cohort
//' @param Xl_test genotypes matrix for large effects in test cohort
//' @param Xs_test genotypes matrix for small effects in test cohort
//' @param sigma2_s estimated value of sigma^2_s
//' @param y_training trait values vector
//' @return variance of predicted y values

arma::mat calc_asymptotic_variance(arma::mat Xl_training, 
                                   arma::mat Xs_training, 
                                   arma::mat Xl_test, 
                                   arma::mat Xs_test,
                                   double sigma2_s,
                                   arma::vec y_training){
  arma::mat Hinv = calc_Hinv(Xs_training, sigma2_s);
  arma::mat var_bl = calc_var_betal(Xl_training, 
                                    Hinv, 
                                    y_training);
  arma::mat var_bs = calc_var_betas(Xl_training, 
                                    Xs_training, 
                                    Hinv, 
                                    sigma2_s,
                                    y_training, 
                                    var_bl);
  arma::mat result = Xl_test.t() * var_bl * Xl_test + Xs_test.t() * var_bs * Xs_test;
  return(result);
}


//' Calculate variance of coefficient estimator for large effects
//' 
//' @return covariance matrix

arma::mat calc_var_betal(arma::mat Sigma_ll, 
                      arma::mat Sigma_ls, 
                      arma::mat Sigma_ss, 
                      double sigma2_s, 
                      unsigned int n){
  //set m_s
  unsigned int m_s = Sigma_ss.n_rows;
  arma::mat small = arma::eye(m_s, m_s) / (n * sigma2_s) + Sigma_ss;
  arma::mat big = Sigma_ll - Sigma_ls * arma::inv(small) * arma::trans(Sigma_ls);
  arma::mat result = arma::inv(big) / n;
  return result;
}

//' Calculate variance of coefficient estimator for small effects
//' 
//' @param Xs matrix of genotypes for small effect markers
//' @param Xl matrix of genotypes for large effect markers
//' @param Hinv inverse of H matrix
//' @param sigma2_s estimated value of sigma^2_s
//' @param y trait values vector
//' @param var_bl variance of beta hat l
//' @return covariance matrix
  
arma::mat calc_var_betas(arma::mat Xl, 
                         arma::mat Xs,
                         arma::mat Hinv,
                         double sigma2_s,
                         arma::vec y,
                         arma::mat var_bl){
  //var y
  arma::mat vy = y * y.t(); 
  // 
  arma::mat arg1 = sigma2_s * sigma2_s * Xs.t() * Hinv * vy * Hinv * Xs;
  arma::mat result = arg1 + sigma2_s * sigma2_s * Xs.t() * Hinv * Xl * var_bl * Xl.t() * Hinv * Xs;
  return result;
}

//' Read pheno data from plink fam file
//' 
//' @param file_path path to plink fam file
//' @param col_number column number in fam file. Default value is 6
//' @return tuple of two string vectors. First is the ids and second is the phenotype, as a string.
//' @references https://techoverflow.net/2020/01/30/how-to-read-tsv-tab-separated-values-in-c/ https://gist.github.com/jbwashington/8b53a7f561322e827c59

std::tuple<vector<string>, vector<string> > read_pheno(std::string file_path, 
                                                       int col_number){
  ifstream fin(file_path);
  string line;
  vector<string> id, pheno;
  while (getline(fin, line)) {
    vector<string> parts;
    boost::algorithm::split(parts, line, boost::algorithm::is_any_of("\t"));
    pheno.push_back(parts[col_number - 1]); // - 1 since indexing starts with zero
    id.push_back(parts[0]);
  }
  //put id and pheno into a single object
  auto result = std::make_tuple (id, pheno);
  return result;
}

//' Convert a string vector containing doubles, as strings, into a numeric vector
//' 
//' @param string_vector a vector containing numeric entries as strings
//' @return numeric vector

std::vector<double> convert_string_vector_to_double_vector(vector<string> string_vector){
  std::vector<double> double_vector(string_vector.size());
  std::transform(string_vector.begin(), string_vector.end(), double_vector.begin(), [](const std::string& val)
  {
    return std::stod(val);
  });
  return(double_vector);
}

//http://arma.sourceforge.net/docs.html#conv_to conv_to for converting between 
// std::vector and arma::vec

//' Mean-center a vector
//' 
//' @param vector
//' @return mean-centered vector

arma::vec center_vector(arma::vec vector){
  arma::vec result = vector - arma::mean(vector);
  return result;
}
