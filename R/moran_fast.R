#' Calculate Moran's I quickly for point data
#' @param x a numeric vector, the value distributed over space
#' @param c1 a numeric vector, the first (x) value of a column of coordinates
#' @param c2 a numeric vector, the second (y) value of a column of coordinates
moran_fast <- function(x, c1, c2){
  res <- calc_moran(x, c1, c2)

  return(res)
}
# 
# Moran.I <- function (x, weight, scaled = FALSE, na.rm = FALSE, alternative = "two.sided")
# {
#     if (dim(weight)[1] != dim(weight)[2])
#         stop("'weight' must be a square matrix")
#     n <- length(x)
#     if (dim(weight)[1] != n)
#         stop("'weight' must have as many rows as observations in 'x'")
#     ei <- -1/(n - 1)
#     nas <- is.na(x)
#     if (any(nas)) {
#         if (na.rm) {
#             x <- x[!nas]
#             n <- length(x)
#             weight <- weight[!nas, !nas]
#         }
#         else {
#             warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
#             return(list(observed = NA, expected = ei, sd = NA,
#                 p.value = NA))
#         }
#     }
#     #ROWSUM <- rowSums(weight)
#     #ROWSUM[ROWSUM == 0] <- 1
#     #weight <- weight/ROWSUM
#     s <- sum(weight)
#     m <- mean(x)
#     y <- x - m
#     cv <- sum(weight * y %o% y)
#     v <- sum(y^2)
#     obs <- (n/s) * (cv/v)
#     if (scaled) {
#         i.max <- (n/s) * (sd(rowSums(weight) * y)/sqrt(v/(n -
#             1)))
#         obs <- obs/i.max
#     }
#     S1 <- 0.5 * sum((weight + t(weight))^2)
#     S2 <- sum((apply(weight, 1, sum) + apply(weight, 2, sum))^2)
#     s.sq <- s^2
#     k <- (sum(y^4)/n)/(v/n)^2
#     sdi <- sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq) -
#         k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq))/((n -
#         1) * (n - 2) * (n - 3) * s.sq) - 1/((n - 1)^2))
#     alternative <- match.arg(alternative, c("two.sided", "less",
#         "greater"))
#     pv <- pnorm(obs, mean = ei, sd = sdi)
#     if (alternative == "two.sided")
#         pv <- if (obs <= ei)
#             2 * pv
#         else 2 * (1 - pv)
#     if (alternative == "greater")
#         pv <- 1 - pv
#     list(observed = obs, expected = ei, sd = sdi, p.value = pv)
# }


