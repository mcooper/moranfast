library(ape)
library(tidyverse)
library(Rcpp)

dat <- read.csv('~/mortalityblob/gbv/GBV_sel.csv') %>%
  filter(country=='AO')

mod <- glm(viol_phys ~ plos_age + woman_literate + is_married + plos_births + 
           plos_hhsize + plos_rural + husband_education_level + plos_husband_age + 
           drought_cat, data=dat, family='binomial')

dat$residuals <- mod$residuals

dmat <- as.matrix(dist(dat[ , c('longitude', 'latitude')]))
dmat <- 1/dmat
dmat[is.infinite(dmat)] <- 0

Moran.I <- function (x, weight, scaled = FALSE, na.rm = FALSE, alternative = "two.sided")
{
    if (dim(weight)[1] != dim(weight)[2])
        stop("'weight' must be a square matrix")
    n <- length(x)
    if (dim(weight)[1] != n)
        stop("'weight' must have as many rows as observations in 'x'")
    ei <- -1/(n - 1)
    nas <- is.na(x)
    if (any(nas)) {
        if (na.rm) {
            x <- x[!nas]
            n <- length(x)
            weight <- weight[!nas, !nas]
        }
        else {
            warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
            return(list(observed = NA, expected = ei, sd = NA,
                p.value = NA))
        }
    }
    #ROWSUM <- rowSums(weight)
    #ROWSUM[ROWSUM == 0] <- 1
    #weight <- weight/ROWSUM
    s <- sum(weight)
    m <- mean(x)
    y <- x - m
    cv <- sum(weight * y %o% y)
    v <- sum(y^2)
    obs <- (n/s) * (cv/v)
    if (scaled) {
        i.max <- (n/s) * (sd(rowSums(weight) * y)/sqrt(v/(n -
            1)))
        obs <- obs/i.max
    }
    S1 <- 0.5 * sum((weight + t(weight))^2)
    S2 <- sum((apply(weight, 1, sum) + apply(weight, 2, sum))^2)
    s.sq <- s^2
    k <- (sum(y^4)/n)/(v/n)^2
    sdi <- sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq) -
        k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq))/((n -
        1) * (n - 2) * (n - 3) * s.sq) - 1/((n - 1)^2))
    alternative <- match.arg(alternative, c("two.sided", "less",
        "greater"))
    pv <- pnorm(obs, mean = ei, sd = sdi)
    if (alternative == "two.sided")
        pv <- if (obs <= ei)
            2 * pv
        else 2 * (1 - pv)
    if (alternative == "greater")
        pv <- 1 - pv
    list(observed = obs, expected = ei, sd = sdi, p.value = pv)
}


real <- Moran.I(dat$residuals, dmat)

Rcpp::sourceCpp('~/moran_fast/calc_moran.cpp')

n <- 50000
start <- Sys.time()
calc_moran(dat$residuals[1:n], dat$longitude[1:n], dat$latitude[1:n])
end <- Sys.time()
end - start

df <- data.frame(n=c(10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000))

for (i in 1:nrow(df)){
  start <- Sys.time()
  calc_moran(dat$residuals[1:df$n[i]], dat$longitude[1:df$n[i]], dat$latitude[1:df$n[i]])
  end <- Sys.time()
  df$time[i] <- as.numeric(end-start)
}





