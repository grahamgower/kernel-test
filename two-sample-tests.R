# Copyright (c) 2018 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

require(parallel)

# Pairwise distance matrix.
euclidean.dmatrix <- function(a, b) {
  n <- length(a)
  x <- matrix(nrow=n, ncol=n)
  for (j in 1:n)
    x[,j] <- sqrt((a-a[j])**2 + (b-b[j])**2)
  return (x)
}

# Pairwise distance matrix (Great-circle distance, in km).
greatcircle.dmatrix <- function(lat, lon, degrees=TRUE) {

  # Trig functions use radians.
  if (degrees) {
    lat <- lat * pi / 180
    lon <- lon * pi / 180
  }
  n <- length(lat)
  x <- matrix(nrow=n, ncol=n)

  # Column-wise vector calculations.
  for (j in 1:n)
    x[,j] <- sin(lat)*sin(lat[j]) + cos(lat)*cos(lat[j]) * cos(lon-lon[j])

  # Avoid NaNs due to floating point imprecision.
  x[which(x < -1)] <- -1
  x[which(x > 1)] <- 1

  radius = 6371 # radius of Earth (km)
  d <- radius * acos(x)
  return (d)
}

# Maximum Mean Discrepancy between two distributions.
# Distance between two samples, with kernel matrix x.
# Gretton et al. 2012, "A kernel two-sample test", Journal of Machine Learning Research.
MMD <- function(x, sizes) {

  # This shouldn't happen, and checking makes things much slower.
  #if (any(eigen(x, symmetric=TRUE, only.values=TRUE)$values < .Machine$double.eps))
  #  stop("matrix is not positive definite---cannot be used as a kernel")

  n <- sizes[[1]]
  m <- sizes[[2]]
  Xii <- x[1:n,1:n]
  Xij <- x[1:n,n+1:m]
  Xjj <- x[n+1:m,n+1:m]
  Mii <- 2 * sum(Xii[upper.tri(Xii)]) / (n*n)
  Mij <- sum(Xij) / (n*m)
  Mjj <- 2 * sum(Xjj[upper.tri(Xjj)]) / (m*m)
  return (2*Mij -Mii -Mjj)
}

# kernel distance
kernel.dist <- function(x, sizes, sigma=NA) {
  fn <- function(s) MMD(exp(-x/s), sizes)
  if (is.na(sigma)) {
    # https://normaldeviate.wordpress.com/2012/07/14/modern-two-sample-tests/
    # Optimise over the scaling parameter sigma.  There is one minima, and
    # the Brent method works well.
    best <- optimise(fn, c(.Machine$double.eps, upper=max(x)))
    return (-best$objective)
  } else {
    return (-fn(sigma))
  }
}


# Kernel test for the distance between two distributions.
# The p-value is derived from a permutation test.
kernel.test <- function(x, sizes, R=1000, kernel="Gaussian", sigma=NA) {

  dname <- deparse(substitute(x))

  if (!isSymmetric(x))
    stop("matrix is not symmetric")
  if (any(x < 0))
    stop("matrix contains negative entries---cannot be a proper metric")
  # The diagonal entries should be zero, but calculations with finite precision
  # can yield non-zero numbers. Allow some wiggle room here.
  if (any(diag(x) > 0.001))
    stop("matrix has non-zero trace---cannot be a proper metric")

  ker <- match.arg(kernel, c("Laplacian", "Gaussian"))
  if (ker == "Gaussian") {
    x <- x*x
    sigma <- sigma*sigma
  }

  D <- kernel.dist(x, sizes, sigma)

  fn <- function(k) {
    shuffle <- sample(nrow(x))
    x2 <- x[shuffle,shuffle] # shuffle rows and columns
    D2 <- kernel.dist(x2, sizes, sigma)
    return (D2)
  }
  D.perms <- unlist(mclapply(1:R, fn))
  p <- sum(D.perms > D) / R

  method <- paste("Two-sample (",ker,") kernel test", sep="")
  alternative <- "samples come from different distributions"
  names(D) = "T"

  ht = list(p.value=p, statistic=D, alternative=alternative,
            data.name=dname,
            method=method,
            permutations=D.perms)
  class(ht) <- "htest"
  return (ht)
}

# Cramer-von Mises distance.
# Anderson, 1962, doi://10.1214/aoms/1177704477
cvm.dist <- function(a, b) {
  e.a <- ecdf(a)
  e.b <- ecdf(b)
  sqdiff <- function (x) (e.a(x)-e.b(x))^2
  T2 <- sum(sqdiff(a))+sum(sqdiff(b))
  #T2 <- T2 * length(a)*length(b)/(length(a)+length(b))^2
  return (T2)
}

# Cramer-von Mises two-sample test.
# The p-value is derived from a permutation test.
# This is *not* the same as the Cramer test of Baringhaus and Franz (2004)
# as implemented in cramer.test from the cramer package.
cvm.test <- function(a, b, R=1000) {

  dname <- paste(deparse(substitute(a)), "and", deparse(substitute(b)))

  C <- cvm.dist(a, b)
  ab <- c(a, b)
  n <- length(ab)
  m <- length(a)
  fn <- function(k) {
    ab2 <- ab[sample(n)]
    a2 <- ab2[1:m]
    b2 <- ab2[(m+1):n]
    C2 <- cvm.dist(a2, b2)
    return (C2)
  }
  C.perms <- unlist(mclapply(1:R, fn))
  p <- sum(C.perms > C) / R

  method <- "Two-sample Cramer-von Mises test"
  alternative <- "samples come from different distributions"
  names(C) = "T"

  ht <- list(p.value=p, statistic=C, alternative=alternative,
            data.name=dname,
            method=method,
            permutations=C.perms)
  class(ht) <- "htest"
  return (ht)
}
