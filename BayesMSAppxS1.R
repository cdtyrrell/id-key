# Tyrrell—Applications in Plant Sciences 2018 #(#)—Data Supplement S1
#
# A method for implementing continuous characters in 
# interactive identification keys and estimating the 
# degree of belief in the taxonomic annotation
# (c) Christopher D. Tyrrell, 2018

###### Methods Example ######
spp.A <- dnorm(x = 3, mean = 4, sd = 1)
spp.B <- dnorm(x = 3, mean = 6, sd = 1)

P.sppA <- spp.A / (spp.A + spp.B)
P.sppB <- spp.B / (spp.A + spp.B)

###### EXAMPLE 1: Invasive vs. native North American Myriophyllum ######

# Appendix 1. Data from Moody & Les 2007
appendix1 <- cbind(c(9:39), 
                   c(1,2,4,2,2,3,12,14,17,14,18,8,14,8,4,0,0,0,1,rep(0,12)),
                   c(rep(0,14),1,1,2,6,9,14,16,22,23,29,26,26,17,12,5,6,1),
                   c(rep(0,7),1,1,2,10,10,11,11,19,27,30,32,45,26,28,11,21,7,6,7,8,3,3,2,1))

# Expand frequency table (Appendix 1) to get raw data
sibiricum <- rep(appendix1[,1], appendix1[,2])
spicatum <- rep(appendix1[,1], appendix1[,3])
hybrid <- rep(appendix1[,1], appendix1[,4])

# Normality tests
shapiro.test(sibiricum)
shapiro.test(spicatum)
shapiro.test(hybrid)

## Figure 1 ##
# Generate data for plotting
lfseg.min <- 9  # min(appendix1[,1])
lfseg.max <- 39  # max(appendix1[,1])
x <- seq(from = lfseg.min, to = lfseg.max, length.out = 301)  # 301 makes nice 0.1 intervals
y.sibiricum <- dnorm(x, mean = mean(sibiricum), sd = sd(sibiricum))
y.spicatum <- dnorm(x, mean = mean(spicatum), sd = sd(spicatum))
y.hybrid <- dnorm(x, mean = mean(hybrid), sd = sd(hybrid))
y <- y.sibiricum + y.spicatum + y.hybrid
plot(c(lfseg.min, lfseg.max), c(0, max(c(y.sibiricum, y.spicatum, y.hybrid))), type = "n",
     xlab = "No. of Leaf Segments", ylab = "Probability Density")
lines(x, y.sibiricum, lwd = 4, lty = 1)  # solid black
lines(x, y.spicatum, lwd = 4, lty = 2)  # dashed black
lines(x, y.hybrid, lwd = 4, lty = 1, col = "grey")  # solid gray

## What is the probability for 23 leaf segments?
# calculate proportions
rP.sibiricum <- y.sibiricum / y
rP.spicatum <- y.spicatum / y
rP.hybrid <- y.hybrid / y
# sample proportions at x = 23 (index 141)
rP.sibiricum[which(x == 23)]
rP.spicatum[which(x == 23)]
rP.hybrid[which(x == 23)]


## Figure 2 ##
# Generate absolute likelihoods for each taxon across leaf segment
# states, given a measurement error (precision)
measerr <- 1 # assume a count precision of +/- 1 segment

P.sibiricum <- P.spicatum <- P.hybrid <- NULL
for(i in x) {
  # I am using a cheater method of integration here becuase I know I
  # have normal distributions. With more complex PDFs one should use
  # the integrate() function or similar.
  tmp <- pnorm(c(i - measerr, i + measerr), mean = mean(sibiricum), sd = sd(sibiricum))
  P.sibiricum <- c(P.sibiricum, tmp[2] - tmp[1])
  tmp <- pnorm(c(i - measerr, i + measerr), mean = mean(spicatum), sd = sd(spicatum))
  P.spicatum <- c(P.spicatum, tmp[2] - tmp[1])
  tmp <- pnorm(c(i - measerr, i + measerr), mean = mean(hybrid), sd = sd(hybrid))
  P.hybrid <- c(P.hybrid, tmp[2] - tmp[1])
}
rm(tmp)

# Calcuate Jaynes evidence
y.sibiricum <- 10 * log10(P.sibiricum / (P.spicatum + P.hybrid))
y.spicatum <- 10 * log10(P.spicatum / (P.sibiricum + P.hybrid))
y.hybrid <- 10 * log10(P.hybrid / (P.sibiricum + P.spicatum))
y <- c(y.sibiricum, y.spicatum, y.hybrid)

plot(c(lfseg.min, lfseg.max), c(min(y), max(y)), 
     type = "n", xlab = "No. of Leaf Segments", ylab = "Evidence (dB)")
lines(x, y.sibiricum, lwd = 3)  # solid black
lines(x, y.spicatum, lty = 2, lwd = 3)  # dashed black
lines(x, y.hybrid, col = "grey", lwd = 3)  # solid grey
abline(h = 0, lty = 3)  # dotted

## Where are the "zones of ambiguity"?
#x[which(y.sibiricum < 1 & y.sibiricum > -1)]  # not run because these are
#x[which(y.spicatum < 1 & y.spicatum > -1)]  # repeated in y.hybrid below
x[which(y.hybrid < 1 & y.hybrid > -1)]


###### EXAMPLE 2: Vegetative Rhipidocladum bamboos in Mexico ######

# Appendix 2.
appendix2.means <- matrix(c(7.17, 34.14, 4.00, 150.0, 10.60, 24.71, 5.18, 43.05),
                          nrow = 4, ncol = 2, byrow = TRUE)
appendix2.sd <- matrix(c(1.53, 29.3, 0.50, 25.0, 4.94, 10.0, 1.48, 28.5),
                       nrow = 4, ncol = 2, byrow = TRUE)
appendix2.range <- list(branches = matrix(c(45, 60, 100, 200, 20, 40, 60, 80),
                                          nrow = 4, ncol = 2, byrow = TRUE),
                        widths = matrix(c(2, 4, 6, 11, 6, 17, 3, 9),
                                        nrow = 4, ncol = 2, byrow = TRUE))
# rownames(appendix2.means) <- rownames(appendix2.sd) <- c("R. bartlettii", "R. martinezii", "R. pittieri", "R. racemiflorum")	
# colnames(appendix2.means) <- colnames(appendix2.sd) <- c("Leaf width (mm)", "No. branches")

## Figure 3 ##
layout(matrix(1:2, nrow = 2, ncol = 1))
par(mar = c(4, 4, 2, 4))  # looks better

# first character
x1 <- seq(from = 0, to = 20, length.out = 201)  # 201 gives nice intervals
y1.bart <- dnorm(x1, mean = appendix2.means[1,1], sd = appendix2.sd[1,1])
y1.mart <- dnorm(x1, mean = appendix2.means[2,1], sd = appendix2.sd[2,1])
y1.pitt <- dnorm(x1, mean = appendix2.means[3,1], sd = appendix2.sd[3,1])
y1.race <- dnorm(x1, mean = appendix2.means[4,1], sd = appendix2.sd[4,1])
y1 <- c(y1.bart, y1.mart, y1.pitt, y1.race)
plot(c(min(x1), max(x1)), c(0, max(y1)),
     xlab = "Leaf Width (mm)", ylab = "Probability Density", type = "n")
lines(x1, y1.bart, lwd = 4, lty = 1)  # solid black
lines(x1, y1.mart, lwd = 4, lty = 4)  # dot-dash black
lines(x1, y1.pitt, lwd = 4, lty = 2)  # dashed black
lines(x1, y1.race, lwd = 4, lty = 1, col = "grey")  # solid grey

# second character
x2 <- 1:201
y2.bart <- dnorm(x2, mean = appendix2.means[1,2], sd = appendix2.sd[1,2])
y2.mart <- dnorm(x2, mean = appendix2.means[2,2], sd = appendix2.sd[2,2])
y2.pitt <- dnorm(x2, mean = appendix2.means[3,2], sd = appendix2.sd[3,2])
y2.race <- dnorm(x2, mean = appendix2.means[4,2], sd = appendix2.sd[4,2])
y2 <- c(y2.bart, y2.mart, y2.pitt, y2.race)
plot(c(min(x2), max(x2)), c(0, max(y2)),
     xlab = "No. Branches", ylab = "Probability Density", type = "n")
lines(x2, y2.bart, lwd = 4, lty = 1)  # solid black
lines(x2, y2.mart, lwd = 4, lty = 4)  # dot-dash black
lines(x2, y2.pitt, lwd = 4, lty = 2)  # dashed black
lines(x2, y2.race, lwd = 4, lty = 1, col = "grey")  # solid grey

# Supposed unknown specimen example
prior <- 1/4  # all species have same initial probabilities

lw <- which(x1 == 10)  # leaf width = 10 mm at index 101
y1.sum <- sum(y1.bart[lw], y1.mart[lw], y1.pitt[lw], y1.race[lw])
bart.lw <- y1.bart[lw] / y1.sum
bart.lw * prior
mart.lw <- y1.mart[lw] / y1.sum
mart.lw * prior
pitt.lw <- y1.pitt[lw] / y1.sum
pitt.lw  * prior
race.lw <- y1.race[lw] / y1.sum
race.lw  * prior

br <- which(x2 == 100)  # no. branches = 100 at index 100
y2.sum <- sum(y2.bart[br], y2.mart[br], y2.pitt[br], y2.race[br])
bart.br <- y2.bart[br] / y2.sum
bart.br * prior
mart.br <- y2.mart[br] / y2.sum
mart.br * prior
pitt.br <- y2.pitt[br] / y2.sum
pitt.br * prior
race.br <- y2.race[br] / y2.sum
race.br * prior

## Figure 4 ##
# Calculate absolute likelihoods for each taxon across each character,
# given some measurement precision
measerr <- c(0.5, 1)  # assume a precision of +/- 0.5 mm width or +/- 1 branch
P.bamboo <- list(P.bart = matrix(NA, nrow = 201, ncol = 201),
                 P.mart = matrix(NA, nrow = 201, ncol = 201),
                 P.pitt = matrix(NA, nrow = 201, ncol = 201),
                 P.race = matrix(NA, nrow = 201, ncol = 201))
for (taxon in 1:4) {
  rcount <- 0
  for(i in x1) {
    ccount <- 0
    rcount <- rcount + 1
    for (j in x2) {
      ccount <- ccount + 1
      # Cheater integration again
      tmp1 <- pnorm(c(i - measerr[1], i + measerr[1]), mean = appendix2.means[taxon, 1], sd = appendix2.sd[taxon, 1])
      tmp2 <- pnorm(c(j - measerr[2], j + measerr[2]), mean = appendix2.means[taxon, 2], sd = appendix2.sd[taxon, 2])
      P.bamboo[[taxon]][rcount, ccount] <- (tmp1[2] - tmp1[1]) * (tmp2[2] - tmp2[1])
    }
  }
}
rm(list = c("tmp1", "tmp2"))

# Calculate Jaynes evidence
E.bart <- 10 * log10(P.bamboo$P.bart / (P.bamboo$P.mart + P.bamboo$P.pitt + P.bamboo$P.race))
E.mart <- 10 * log10(P.bamboo$P.mart / (P.bamboo$P.bart + P.bamboo$P.pitt + P.bamboo$P.race))
E.pitt <- 10 * log10(P.bamboo$P.pitt / (P.bamboo$P.bart + P.bamboo$P.mart + P.bamboo$P.race))
E.race <- 10 * log10(P.bamboo$P.race / (P.bamboo$P.bart + P.bamboo$P.mart + P.bamboo$P.pitt))

# Need to "downsample" evidence to allow gridlines on plots
sampleseq <- seq(from = 1, to = 202, by = 10)
samplegrid <- expand.grid(sampleseq, sampleseq)
e.bart <- matrix(E.bart[as.matrix(samplegrid)], nrow = 21, ncol = 21)
e.mart <- matrix(E.mart[as.matrix(samplegrid)], nrow = 21, ncol = 21)
e.pitt <- matrix(E.pitt[as.matrix(samplegrid)], nrow = 21, ncol = 21)
e.race <- matrix(E.race[as.matrix(samplegrid)], nrow = 21, ncol = 21)

layout(matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2))
par(mai = c(0, 0.1, 0, 0))  # shrink margins; looks better
# Fig. 4A
persp(x = sampleseq, y = sampleseq, z = e.bart, 
      theta = 330, phi = 20, col = "grey", border = "white", xlab = "No. Branches", 
      ylab = "Leaf Width (mm)", zlab = "Evidence (dB)",
      zlim = c(-226, 141) ) -> transmtx
xy <- trans3d(x = c(1, 1, 202), 
              y = c(202, 1, 1), 
              z = 0, transmtx)
lines(xy$x, xy$y, lty = 2)
legend("topleft", legend = "", title = "A", bty = "n", cex = 2)
# Fig. 4B
persp(x = sampleseq, y = sampleseq, z = e.mart, 
      theta = 330, phi = 20, col = "grey", border = "white", xlab = "No. Branches", 
      ylab = "Leaf Width (mm)", zlab = "Evidence (dB)",
      zlim = c(-226, 141))
lines(xy$x, xy$y, lty = 2)
legend("topleft", legend = "", title = "B", bty = "n", cex = 2)
# Fig. 4C
persp(x = sampleseq, y = sampleseq, z = e.pitt, 
      theta = 330, phi = 20, col = "grey", border = "white", xlab = "No. Branches", 
      ylab = "Leaf Width (mm)", zlab = "Evidence (dB)",
      zlim = c(-226, 141))
lines(xy$x, xy$y, lty = 2)
legend("topleft", legend = "", title = "C", bty = "n", cex = 2)
# Fig. 4D
persp(x = sampleseq, y = sampleseq, z = e.race, 
      theta = 330, phi = 20, col = "grey", border = "white", xlab = "No. Branches", 
      ylab = "Leaf Width (mm)", zlab = "Evidence (dB)",
      zlim = c(-226, 141))
lines(xy$x, xy$y, lty = 2)
legend("topleft", legend = "", title = "D", bty = "n", cex = 2)

## Figure 5 ##
# clear values less than zero
E.bart.pos <- E.bart
E.mart.pos <- E.mart
E.pitt.pos <- E.pitt
E.race.pos <- E.race
E.bart.pos[E.bart.pos < 0] <- -1
E.mart.pos[E.mart.pos < 0] <- -1
E.pitt.pos[E.pitt.pos < 0] <- -1
E.race.pos[E.race.pos < 0] <- -1

dev.off() # reset par & layout
plot(c(min(x2), max(x2)), c(min(x1), max(x1)), type = "n", xlab = "No. of branches", ylab = "Leaf width (mm)")
contour(x2, x1, t(E.bart.pos), add=TRUE, lwd = 2, lty = 1)  # solid black
contour(x2, x1, t(E.mart.pos), add=TRUE, lwd = 2, lty = 4)  # dot-dash black
contour(x2, x1, t(E.pitt.pos), add=TRUE, lwd = 2, lty = 2)  # dashed black
contour(x2, x1, t(E.race.pos), add=TRUE, lwd = 2, col = "grey")  # solid grey

#ambigzone <- E.bart.pos + E.mart.pos + E.pitt.pos + E.race.pos
#pts <- which(ambigzone == -4, arr.ind = TRUE)