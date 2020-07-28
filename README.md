# moranfast
Code to calculate Moran's I in R quickly and memory efficiently

This package is still in development.  Right now it just does one thing (calculates Moran's I).

`moranfast` is an improvement over any other package I know of for calculating Moran's I for two reasons:

1. It is _memory efficient_, because it calculates the distance matrix on-the-fly.  It shouldnt take up any more memory than it takes to hold a dataframe of point observations in R.
2. It is _fast_, because it uses Rcpp. I found it calculated the Moran's I for 100,000 observations in under a minute.
