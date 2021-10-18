#include <tuple> // std::tuple, std::get, std::tie, std::ignore
#include <armadillo>
#include <math.h>       /* floor */

using namespace std;

//' Calculate qplus alpha and qminus alpha for a vector of LOO residuals
//' 
//' @param residuals vector of leave-one out residuals
//' @param alpha
//' @return a tuple of length two, with each element containing exactly one entry, an entry from residuals vector, corresponding to the empirical quantiles of interest
//' @reference see page 4 of Barber et al 2020 preprint "Predictive inference with the jackknife+"

std::tuple<double, double> calc_q(arma::vec residuals, 
                                  double alpha){
  // sort residuals
  std::vector<double> resid_std = arma::conv_to<std::vector<double> >::from(residuals);
  std::sort(resid_std.begin(), resid_std.end()); //essentially overwrites resid_std with the sorted vector, smallest to largest
  //resid_std is now ordered from least to greatest
  int lower_index = floor(alpha * (1 + resid_std.size()));
  int upper_index = ceil((1 - alpha) * (1 + resid_std.size()));
  //store results as a tuple
  std::tuple<double, double> result = std::make_tuple (resid_std[lower_index - 1], resid_std[upper_index - 1]); //-1 due to indexing starting with zero
  return(result);
}


