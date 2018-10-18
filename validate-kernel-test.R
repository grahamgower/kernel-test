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

library(dplyr)
library(MASS)
library(energy)
library(cramer)

source('two-sample-tests.R')


# EU
mu1 <- c(55,20) # pop. centre
s1_m <- matrix(c(20,-10,-10,60),ncol=2) # cov. matrix
s1_f <- matrix(c(10,10,10,15),ncol=2) # cov. matrix (alternative hypothesis)

# AM
mu2 <- c(65,220) # pop. centre
s2_m <- matrix(c(50,-115,-115,1150),ncol=2) # cov. matrix
s2_f <- matrix(c(30,-115,-115,500),ncol=2) # cov. matrix (alternative hypothesis)

# observed bison sample numbers
pops <- c(45,132)
sexes <- c(43,134)

set.seed(31415)
dsimlist <- list()
j <- 1
nsims <- 1000
nperms = 1000

##
## simulate locations for ndraws populations
for (hypothesis in c("Ha", "H0")) {

  if (hypothesis == "H0") {
    s1_f <- s1_m
    s2_f <- s2_m
  }

  for (k in 1:nsims) {

    # draw new sample numbers
    pop1 <- sum(rbinom(sum(pops), 1, pops[1]/sum(pops)))
    pop2 <- sum(pops) - pop1
    sex1p1 <- sum(rbinom(pop1, 1, sexes[1]/sum(sexes)))
    sex2p1 <- pop1 - sex1p1
    sex1p2 <- sum(rbinom(pop2, 1, sexes[1]/sum(sexes)))
    sex2p2 <- pop2 - sex1p2

    dsim <- data.frame(rbind(
      mvrnorm(sex1p1, mu1, s1_f),
      mvrnorm(sex2p1, mu1, s1_m),
      mvrnorm(sex1p2, mu2, s2_f),
      mvrnorm(sex2p2, mu2, s2_m)))
    colnames(dsim) <- c("lat","lon")
    dsim$sex <- factor(c(rep("F",sex1p1),rep("M",sex2p1),rep("F",sex1p2),rep("M",sex2p2)))
    dsim <- dsim[order(dsim$sex),] # order by sex

    dsimlist[[j]] <- dsim
    j <- j+1
  }
}


nah <- rep(NA, 11*nsims)
res <- tibble(test = as.character(nah),
              h = as.character(nah),
              k = as.integer(nah),
              p = as.numeric(nah))
i <- 1
j <- 1

for (hypothesis in c("Ha", "H0")) {
  print(hypothesis)
  for (k in 1:nsims) {
    print(k)

    dsim <- dsimlist[[j]]
    j <- j+1

    gcm <- greatcircle.dmatrix(dsim$lat, dsim$lon)
    eum <- euclidean.dmatrix(dsim$lat, dsim$lon)

    # tests for differences in spatial distribution between sexes

    # energy distance test (Szekely and Rizzo, 2004)
    e1 <- eqdist.etest(gcm, table(dsim$sex), distance=T, R=nperms)
    res[i,] <- list(test="Energy distance ($d_{gc}$)", h=hypothesis, k=k, p=e1$p.value)
    i <- i + 1

    e2 <- eqdist.etest(eum, table(dsim$sex), distance=T, R=nperms)
    res[i,] <- list(test="Energy distance ($d_{E}$)", h=hypothesis, k=k, p=e2$p.value)
    i <- i + 1

    # Cramer test (Baringhaus and Franz, 2004)
    cr <- cramer.test(as.matrix(subset(dsim,sex=="M")[c("lat","lon")]),
                as.matrix(subset(dsim,sex=="F")[c("lat","lon")]),
                sim="permutation", replicates=nperms)
    res[i,] <- list(test="Cram\\'{e}r test", h=hypothesis, k=k, p=cr$p.value)
    i <- i + 1

    # univariate Kolmogorov-Smirnov test, latitude
    ks.lat <- ks.test(subset(dsim,sex=="M")$lat, subset(dsim,sex=="F")$lat)
    res[i,] <- list(test="Kolmogorov-Smirnov (lat)", h=hypothesis, k=k, p=ks.lat$p.value)
    i <- i + 1

    # univariate Kolmogorov-Smirnov test, longitude
    ks.lon <- ks.test(subset(dsim,sex=="M")$lon, subset(dsim,sex=="F")$lon)
    res[i,] <- list(test="Kolmogorov-Smirnov (lon)", h=hypothesis, k=k, p=ks.lon$p.value)
    i <- i + 1

    # univariate Cramer-von Mises test, latitude
    cvm.lat <- cvm.test(subset(dsim,sex=="M")$lat, subset(dsim,sex=="F")$lat, R=nperms)
    res[i,] <- list(test="Cram\\'{e}r-von Mises (lat)", h=hypothesis, k=k, p=cvm.lat$p.value)
    i <- i + 1

    # univariate Cramer-von Mises test, longitude
    cvm.lon <- cvm.test(subset(dsim,sex=="M")$lon, subset(dsim,sex=="F")$lon, R=nperms)
    res[i,] <- list(test="Cram\\'{e}r-von Mises (lon)", h=hypothesis, k=k, p=cvm.lon$p.value)
    i <- i + 1

    # kernel test, Laplacian kernel
    kern1 <- kernel.test(gcm, table(dsim$sex), kernel="Laplacian", R=nperms)
    res[i,] <- list(test="kernel test ($k_L,d_{gc}$)", h=hypothesis, k=k, p=kern1$p.value)
    i <- i + 1

    # kernel test, Gaussian kernel
    kern2 <- kernel.test(gcm, table(dsim$sex), kernel="Gaussian", R=nperms)
    res[i,] <- list(test="kernel test ($k_G,d_{gc}$)", h=hypothesis, k=k, p=kern2$p.value)
    i <- i + 1

    # kernel test, Laplacian kernel
    kern3 <- kernel.test(eum, table(dsim$sex), kernel="Laplacian", R=nperms)
    res[i,] <- list(test="kernel test ($k_L,d_{E}$)", h=hypothesis, k=k, p=kern3$p.value)
    i <- i + 1

    # kernel test, Gaussian kernel
    kern4 <- kernel.test(eum, table(dsim$sex), kernel="Gaussian", R=nperms)
    res[i,] <- list(test="kernel test ($k_G,d_{E}$)", h=hypothesis, k=k, p=kern4$p.value)
    i <- i + 1

    # kernel test, Laplacian kernel, fixed bandwidth
    kern5 <- kernel.test(gcm, table(dsim$sex), kernel="Laplacian", R=nperms, sigma=median(gcm))
    res[i,] <- list(test="kernel test ($k_L,d_{gc},\\sigma=median$)", h=hypothesis, k=k, p=kern5$p.value)
    i <- i + 1

    # kernel test, Gaussian kernel, fixed bandwidth
    kern6 <- kernel.test(gcm, table(dsim$sex), kernel="Gaussian", R=nperms, sigma=median(gcm))
    res[i,] <- list(test="kernel test ($k_G,d_{gc},\\sigma=median$)", h=hypothesis, k=k, p=kern6$p.value)
    i <- i + 1

    # kernel test, Laplacian kernel, fixed bandwidth
    kern7 <- kernel.test(eum, table(dsim$sex), kernel="Laplacian", R=nperms, sigma=median(eum))
    res[i,] <- list(test="kernel test ($k_L,d_{E},\\sigma=median$)", h=hypothesis, k=k, p=kern7$p.value)
    i <- i + 1

    # kernel test, Gaussian kernel, fixed bandwidth
    kern8 <- kernel.test(eum, table(dsim$sex), kernel="Gaussian", R=nperms, sigma=median(eum))
    res[i,] <- list(test="kernel test ($k_G,d_{E},\\sigma=median$)", h=hypothesis, k=k, p=kern8$p.value)
    i <- i + 1
  }

}

warnings()

res2 <- summarise(group_by(res, test),
                 type1 = sum(h=="H0" & p<0.05)/nsims,
                 power = sum(h=="Ha" & p<0.05)/nsims)

# for copy/paste into LaTeX
write.table(res2, row.names=FALSE, quote=FALSE, sep=" & ", eol=" \\\\\n")


#kk <- kernel.test(gcm, table(dsim$sex), R=10000)
#hist(kk$perms, 100)
#abline(v=kk$D, col="red")

#require(ggplot2)
#ggplot(dsim, aes(x=lon, y=lat, colour=sex, shape=sex)) +
#  geom_jitter(width=.2, height=.2, size=2) +
#  scale_shape_manual(values=c(1,4)) +
#  scale_color_manual(values=c("red", "blue")) +
#  geom_density_2d()
