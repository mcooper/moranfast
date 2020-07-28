#include <Rcpp.h>
using namespace Rcpp;

// This function is useful for calculating Moran's I 
// Without holding an entire distance matrix in memory

// [[Rcpp::export]]
double distanceCalculate(double x1, double y1, double x2, double y2)
{
	double x = x1 - x2; //calculating number to square in next step
	double y = y1 - y2;
	double dist;

	dist = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
	dist = sqrt(dist);                  
  
  if (dist > 0){
    dist = 1/dist;
  }

	return dist;
}

// [[Rcpp::export]]
NumericVector normalize(NumericVector x)
{
  double x_bar = mean(x);
  NumericVector x_norm = x - x_bar;

  return x_norm;
}

// [[Rcpp::export]]
double calc_moran(NumericVector x, NumericVector c1, NumericVector c2)
{
  // Easy variables to calculate
  NumericVector x_norm = normalize(x);
  double N = x.length();
  double denom = sum(x_norm*x_norm);
 
  Rcout << "The value of N : " << N << "\n";
  Rcout << "The value of denom : " << denom << "\n";

  // Variables to calculate through iteration
  double W = 0;
  double num = 0;

  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j){
      double w_ij = distanceCalculate(c1[i], c2[i], c1[j], c2[j]);
      W += w_ij;

      double w_x1_x2 = w_ij * x_norm[i] * x_norm[j];
      num += w_x1_x2;
    }
  }
 
  Rcout << "The value of W : " << W << "\n";
  Rcout << "The value of num : " << num << "\n";

  double I = (N/W)*(num/denom);

  return(I);
}
