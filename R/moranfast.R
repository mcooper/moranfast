#' Calculate Moran's I quickly for point data
#' @param x a numeric vector, the value distributed over space
#' @param c1 a numeric vector, the first (x) value of a column of coordinates
#' @param c2 a numeric vector, the second (y) value of a column of coordinates
#' @param alternative a character sring specifying the alternative hypothesis that is tested against; must be one of "two.sided", "less", or "greater", or any unambiguous abbreviation of these.

moranfast <- function(x, c1, c2, alternative='two.sided'){
  res <- calc_moran(x, c1, c2)
  
  names(res) <- c('observed', 'expected', 'sd')
  res <- as.list(res)
    
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  pv <- pnorm(res$observed, mean = res$expected, sd = res$sd)
  if (alternative == "two.sided"){
    if (res$observed <= -1/(length(x) - 1)){
      pv <- 2 * pv
    }else{
      pv <- 2 * (1 - pv)
    }
  }
  if (alternative == "greater"){
    pv <- 1 - pv
  }
  
  res[['p.value']] <- pv
  
  return(res)
}


