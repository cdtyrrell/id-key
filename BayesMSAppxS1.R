# Tyrrell—Applications in Plant Sciences 2018-D-18-00191 - Appendix S2
#
# A method for implementing continuous characters in 
# interactive identification keys and estimating the 
# degree of belief in the taxonomic annotation
# (c) Christopher D. Tyrrell, 2018

# APPS figure dimensions, inches
Fig1.5col <- 5.5
Fig2col <- 7.5
FigHeight <- 9
# Universal design color palette
bluegreen <- "#009E73"
reddish <- "#CC79A7"
skyblue <- "#56B4E9"
vermillion <- "#D55E00"

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

# Quantile-quantile plots to visualize possible deviation from normality (Appendix S1)
pdf(file = "Appendix_S1.pdf", width = Fig1.5col, height = FigHeight)
layout(matrix(1:3, nrow = 3, ncol = 1))
qqnorm(sibiricum, pch = 20, main = "", ylab = "Sample quantiles", xlab = "Theoretical quantiles")
qqline(sibiricum)
mtext("A", side = 2, line = 3, las = 1, at = 27)
qqnorm(spicatum, pch = 20, main = "", ylab = "Sample quantiles", xlab = "Theoretical quantiles")
qqline(spicatum)
mtext("B", side = 2, line = 3, las = 1, at = 40)
qqnorm(hybrid, pch = 20, main = "", ylab = "Sample quantiles", xlab = "Theoretical quantiles")
qqline(hybrid)
mtext("C", side = 2, line = 3, las = 1, at = 40)
dev.off()

# Normality tests (not overly useful, but can be run if uncommented)
#shapiro.test(sibiricum)
#shapiro.test(spicatum)
#shapiro.test(hybrid)

## Figure 1 ##
# Generate data for plotting
lfseg.min <- 9  # min(appendix1[,1])
lfseg.max <- 39  # max(appendix1[,1])
x <- seq(from = lfseg.min, to = lfseg.max, length.out = 301)  # 301 makes nice 0.1 intervals
absProb.sibiricum <- dnorm(x, mean = mean(sibiricum), sd = sd(sibiricum))  # generate normal probability densities (curves) across leaf segment states for each taxon
absProb.spicatum <- dnorm(x, mean = mean(spicatum), sd = sd(spicatum))
absProb.hybrid <- dnorm(x, mean = mean(hybrid), sd = sd(hybrid))
plotHeight <- max(c(absProb.sibiricum, absProb.spicatum, absProb.hybrid))
postscript(file = "Figure_1.eps", width = Fig2col, height = FigHeight / 3)
par(mar = c(4,4,1,1))  # more efficient for journal
plot(c(lfseg.min, lfseg.max), c(0, plotHeight), type = "n",
     xlab = "No. of leaf segments", ylab = "Probability density")
lines(x, absProb.sibiricum, lwd = 4, lty = 1, col = vermillion)  # solid orange
lines(x, absProb.spicatum, lwd = 4, lty = 1, col = bluegreen)  # solid green
lines(x, absProb.hybrid, lwd = 4, lty = 2, col = skyblue)  # dashed blue
legend("topleft", legend = c("M. sibiricum", "M. spicatum", "hybrid"), text.font = c(3,3,1), lwd = 2, lty = c(1, 1, 2), col = c(vermillion, bluegreen, skyblue), bty = "n")
dev.off()

## What is the probability for 23 leaf segments?
sumProbs <- absProb.sibiricum + absProb.spicatum + absProb.hybrid
# calculate proportions
relProb.sibiricum <- absProb.sibiricum / sumProbs
relProb.spicatum <- absProb.spicatum / sumProbs
relProb.hybrid <- absProb.hybrid / sumProbs
# sample proportions at x = 23 (index 141)
relProb.sibiricum[which(x == 23)]
relProb.spicatum[which(x == 23)]
relProb.hybrid[which(x == 23)]


## Figure 2 ##
# Generate absolute likelihoods for each taxon across leaf segment
# states, given a measurement error (precision)
measerr <- 1 # assume a count precision of +/- 1 segment

postProb.sibiricum <- postProb.spicatum <- postProb.hybrid <- NULL
for(i in x) {
  # I am using a cheater method of integration here becuase I know I
  # have normal distributions and R has built-in functions for them.
  # Essentially we calculate the area from the left tail to the 
  # measurement (i) plus error (measerr) and subtract that from the
  # area from the left tail to the measurement minus error to get
  # the area we want. Functions for other PDFs/PMFs can be substituted
  # into this method. For more complex PDFs, one could use create their
  # own function and use integrate() or similar.
  tmp <- pnorm(c(i - measerr, i + measerr), mean = mean(sibiricum), sd = sd(sibiricum)) # get both areas at same time
  postProb.sibiricum <- c(postProb.sibiricum, tmp[2] - tmp[1]) # subtract smaller from larger
  tmp <- pnorm(c(i - measerr, i + measerr), mean = mean(spicatum), sd = sd(spicatum))
  postProb.spicatum <- c(postProb.spicatum, tmp[2] - tmp[1])
  tmp <- pnorm(c(i - measerr, i + measerr), mean = mean(hybrid), sd = sd(hybrid))
  postProb.hybrid <- c(postProb.hybrid, tmp[2] - tmp[1])
}
rm(tmp)

# Calcuate Jaynes evidence
evi.sibiricum <- 10 * log10(postProb.sibiricum / (postProb.spicatum + postProb.hybrid))
evi.spicatum <- 10 * log10(postProb.spicatum / (postProb.sibiricum + postProb.hybrid))
evi.hybrid <- 10 * log10(postProb.hybrid / (postProb.sibiricum + postProb.spicatum))
evi.sum <- c(evi.sibiricum, evi.spicatum, evi.hybrid)

## Figure 2 ##
postscript(file = "Figure_2.eps", width = Fig2col, height = FigHeight / 3)
par(mar = c(4,4,1,1))
plot(c(lfseg.min, lfseg.max), c(min(evi.sum), max(evi.sum)), 
     type = "n", xlab = "No. of leaf segments", ylab = "Evidence (dB)")
lines(x, evi.sibiricum, lty = 1, col = vermillion, lwd = 3)  # solid orange
lines(x, evi.spicatum, lty = 1, col = bluegreen, lwd = 3)  # solid green
lines(x, evi.hybrid, lty = 2, col = skyblue, lwd = 3)  # dashed blue
abline(h = 0, lty = 3)  # dotted
legend("bottom", legend = c("M. sibiricum", "M. spicatum", "hybrid"), text.font = c(3,3,1), lwd = 2, lty = c(1, 1, 2), col = c(vermillion, bluegreen, skyblue), bty = "n")
dev.off()

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
# FYI:
# rownames(appendix2.means) <- rownames(appendix2.sd) <- c("R. bartlettii", "R. martinezii", "R. pittieri", "R. racemiflorum")	
# colnames(appendix2.means) <- colnames(appendix2.sd) <- c("Leaf width (mm)", "No. branches")

## Figure 3 ##
postscript(file = "Figure_3.eps", width = Fig2col, height = FigHeight * (2/3))
layout(matrix(1:2, nrow = 2, ncol = 1))
par(mar = c(4, 4, 1, 1))  # looks better
# first character
x1 <- seq(from = 0, to = 20, length.out = 201)  # 201 gives nice intervals
absP1.bart <- dnorm(x1, mean = appendix2.means[1,1], sd = appendix2.sd[1,1])
absP1.mart <- dnorm(x1, mean = appendix2.means[2,1], sd = appendix2.sd[2,1])
absP1.pitt <- dnorm(x1, mean = appendix2.means[3,1], sd = appendix2.sd[3,1])
absP1.race <- dnorm(x1, mean = appendix2.means[4,1], sd = appendix2.sd[4,1])
absP1.all <- c(absP1.bart, absP1.mart, absP1.pitt, absP1.race)
plot(c(min(x1), max(x1)), c(0, max(absP1.all)),
     xlab = "Leaf width (mm)", ylab = "Probability density", type = "n")
lines(x1, absP1.bart, lwd = 4, lty = 1, col = bluegreen)  # solid green
lines(x1, absP1.mart, lwd = 4, lty = 4, col = reddish)  # dot-dash purple
lines(x1, absP1.pitt, lwd = 4, lty = 2, col = skyblue)  # dashed blue
lines(x1, absP1.race, lwd = 4, lty = 1, col = vermillion)  # solid orange
legend("topright", legend = c("R. bartlettii", "R. martinezii", "R. pittieri", "R. racemiflorum"), text.font = 3, lwd = 2, lty = c(1,4,2,1), col = c(bluegreen, reddish, skyblue, vermillion), bty = "n")
mtext("A",2, las = 1, line = 3, at = 0.8)
# second character
x2 <- 1:201
absP2.bart <- dnorm(x2, mean = appendix2.means[1,2], sd = appendix2.sd[1,2])
absP2.mart <- dnorm(x2, mean = appendix2.means[2,2], sd = appendix2.sd[2,2])
absP2.pitt <- dnorm(x2, mean = appendix2.means[3,2], sd = appendix2.sd[3,2])
absP2.race <- dnorm(x2, mean = appendix2.means[4,2], sd = appendix2.sd[4,2])
absP2.all <- c(absP2.bart, absP2.mart, absP2.pitt, absP2.race)
plot(c(min(x2), max(x2)), c(0, max(absP2.all)),
     xlab = "No. of branches", ylab = "Probability density", type = "n")
lines(x2, absP2.bart, lwd = 4, lty = 1, col = bluegreen)  # solid green
lines(x2, absP2.mart, lwd = 4, lty = 4, col = reddish)  # dot-dash purple
lines(x2, absP2.pitt, lwd = 4, lty = 2, col = skyblue)  # dashed blue
lines(x2, absP2.race, lwd = 4, lty = 1, col = vermillion)  # solid orange
mtext("B",2, las = 1, line = 3, at = 0.04)
dev.off()

## Type specimen examples ##
prior <- 1/4  # all species have same initial probabilities
lf.widths.types <- c(7.5, 3, 6, 6.2)
num.branches.types <- c(23, 135, 25, NA)
probs.width <- matrix(NA, nrow = 4, ncol = length(lf.widths.types))
for(lcv in 1:length(lf.widths.types)) {  # lcv = loop control variable
  idx <- which(x1 == lf.widths.types[lcv])  # find the index (idx) where the type specimen's leaf width matches x1
  absP1.sum <- sum(absP1.bart[idx], absP1.mart[idx], absP1.pitt[idx], absP1.race[idx])
  probs.width[1,lcv] <- absP1.bart[idx] / absP1.sum
  probs.width[2,lcv] <- absP1.mart[idx] / absP1.sum
  probs.width[3,lcv] <- absP1.pitt[idx] / absP1.sum
  probs.width[4,lcv] <- absP1.race[idx] / absP1.sum
}
numtypes <- length(num.branches.types)
probs.branch <- matrix(NA, nrow = 4, ncol = numtypes) 
for(lcv in 1:numtypes-1) {  # -1 because of no data for R. racemiflorum
  idx <- which(x2 == num.branches.types[lcv])  # find the index of x2 where number of branches matches the type specimen's
  absP2.sum <- sum(absP2.bart[idx], absP2.mart[idx], absP2.pitt[idx], absP2.race[idx])
  probs.branch[1,lcv] <- absP2.bart[idx] / absP2.sum
  probs.branch[2,lcv] <- absP2.mart[idx] / absP2.sum
  probs.branch[3,lcv] <- absP2.pitt[idx] / absP2.sum
}
rownames(probs.width) <- rownames(probs.branch) <- c("R. bartlettii", "R. martinezii", "R. pittieri", "R. racemiflorum")
postProbs.types <- list(lf.width.only = prior * probs.width, num.branches.only = prior * probs.branch,
                    both = prior * probs.width * probs.branch)
## Table 5 ##
round(postProbs.types$lf.width.only, 3)
round(postProbs.types$num.branches.only, 3)
round(postProbs.types$both, 3)

## Bamboo Posterior Probabilities Across Character States ##
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
      # Cheater integration using the built-in R function
      tmp1 <- pnorm(c(i - measerr[1], i + measerr[1]), mean = appendix2.means[taxon, 1], sd = appendix2.sd[taxon, 1])
      tmp2 <- pnorm(c(j - measerr[2], j + measerr[2]), mean = appendix2.means[taxon, 2], sd = appendix2.sd[taxon, 2])
      # To use a model other than normal, one could substitute
      # the pnorm() function for another of R's built-in PDFs
      # such as plnorm() [log normal], ppois() [Poisson], etc.
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

## Figure 4 ##
# In order to color the perspective plots, must determine the z-axis center (average value) for each "tile"
bart.facets <- (E.bart[-1, -1] + E.bart[-1, -ncol(E.bart)] + E.bart[-nrow(E.bart), -1] + E.bart[-nrow(E.bart), -ncol(E.bart)]) / 4
mart.facets <- (E.mart[-1, -1] + E.mart[-1, -ncol(E.mart)] + E.mart[-nrow(E.mart), -1] + E.mart[-nrow(E.mart), -ncol(E.mart)]) / 4
pitt.facets <- (E.pitt[-1, -1] + E.pitt[-1, -ncol(E.pitt)] + E.pitt[-nrow(E.pitt), -1] + E.pitt[-nrow(E.pitt), -ncol(E.pitt)]) / 4
race.facets <- (E.race[-1, -1] + E.race[-1, -ncol(E.race)] + E.race[-nrow(E.race), -1] + E.race[-nrow(E.race), -ncol(E.race)]) / 4

postscript(file = "Figure_4.eps", width = Fig2col, height = FigHeight * (2/3))
# Layout 4 panels
layout(matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2))
par(mai = c(0, 0.1, 0, 0))  # shrink margins; looks better
# Fig. 4A
cols <- bart.facets
cols[cols <= 0] <- "darkgrey"
cols[as.numeric(cols) > 0] <- bluegreen
persp(x = 0:200, y = 0:200, z = E.bart,
      theta = 330, phi = 20, col = cols, border = NA, ylab = "No. of branches", 
      xlab = "Leaf width (mm)", zlab = "Evidence (dB)",
      zlim = c(-226, 141)) -> transmtx
xy <- trans3d(x = c(1, 1, 202), 
              y = c(202, 1, 1), 
              z = 0, transmtx)
lines(xy$x, xy$y, lty = 2, col = bluegreen, lwd = 2)
legend("topleft", legend = NA, title = "A", bty = "n", cex = 2)
# Fig. 4B
cols <- mart.facets
cols[cols <= 0] <- "darkgrey"
cols[as.numeric(cols) > 0] <- reddish
persp(x = 0:200, y = 0:200, z = E.mart, 
      theta = 330, phi = 20, col = cols, border = NA, ylab = "No. of branches", 
      xlab = "Leaf width (mm)", zlab = "Evidence (dB)",
      zlim = c(-226, 141))
lines(xy$x, xy$y, lty = 2, col = reddish, lwd = 2)
legend("topleft", legend = "", title = "B", bty = "n", cex = 2)
# Fig. 4C
cols <- pitt.facets
cols[cols <= 0] <- "darkgrey"
cols[as.numeric(cols) > 0] <- skyblue
persp(x = 0:200, y = 0:200, z = E.pitt, 
      theta = 330, phi = 20, col = cols, border = NA, ylab = "No. of branches", 
      xlab = "Leaf width (mm)", zlab = "Evidence (dB)",
      zlim = c(-226, 141))
lines(xy$x, xy$y, lty = 2, col = skyblue, lwd = 2)
legend("topleft", legend = "", title = "C", bty = "n", cex = 2)
# Fig. 4D
cols <- race.facets
cols[cols <= 0] <- "darkgrey"
cols[as.numeric(cols) > 0] <- vermillion
persp(x = 0:200, y = 0:200, z = E.race, 
      theta = 330, phi = 20, col = cols, border = NA, ylab = "No. of branches", 
      xlab = "Leaf width (mm)", zlab = "Evidence (dB)",
      zlim = c(-226, 141))
lines(xy$x, xy$y, lty = 2, col = vermillion, lwd = 2)
legend("topleft", legend = "", title = "D", bty = "n", cex = 2)
dev.off()


## Figure 5 ##
# Get only the positive evidence values for each species
E.bart.pos <- E.bart
E.mart.pos <- E.mart
E.pitt.pos <- E.pitt
E.race.pos <- E.race
# set anything less than 0 to -1 for simplicity of plotting
E.bart.pos[E.bart.pos < 0] <- -1
E.mart.pos[E.mart.pos < 0] <- -1
E.pitt.pos[E.pitt.pos < 0] <- -1
E.race.pos[E.race.pos < 0] <- -1
# combine all four species into one evidence "landscape"
eviscape <- tmp <- pmax(E.bart.pos, E.mart.pos, E.pitt.pos, E.race.pos)
# recode the "landscape" for whichever species has maximum evidence at a particular coordinate
eviscape[eviscape == E.bart.pos] <- 1
eviscape[eviscape == E.mart.pos] <- 2
eviscape[eviscape == E.pitt.pos] <- 3
eviscape[eviscape == E.race.pos] <- 4
eviscape[tmp == -1] <- 0  # "zone of ambiguity" i.e., 'holes' in the positive evidence

fig5cols <- c("white", bluegreen, reddish, skyblue, vermillion)
postscript(file = "Figure_5.eps", width = Fig2col, height = Fig2col)
layout(matrix(c(2,4,1,3), ncol=2, byrow=TRUE), widths=c(4/5,1/5), heights=c(1/5,4/5))
# Level plot with contour lines
par(mar=c(4,4,1,1))
plot(c(min(x1), max(x1)), c(min(x2), max(x2)), type = "n", ylab = "No. of branches", xlab = "Leaf width (mm)")
.filled.contour(x1, x2, eviscape, levels = c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5), col = fig5cols)
contour(x1, x2, E.bart.pos, add=TRUE, lwd = 2, lty = 1)  # solid black
contour(x1, x2, E.mart.pos, add=TRUE, lwd = 2, lty = 4)  # dot-dash black
contour(x1, x2, E.pitt.pos, add=TRUE, lwd = 2, lty = 2)  # dashed black
contour(x1, x2, E.race.pos, add=TRUE, lwd = 2, col = "darkgrey")  # solid grey
# Leaf width evidence margin plot 
par(mar=c(0,4,1,1))
plot(x1, apply(E.pitt, 1, max), xaxt = "n", ylim = c(-25,140), ylab = "Evidence (db)", col = skyblue, pch = 20)
points(x1, apply(E.mart, 1, max), col = reddish, add = TRUE, pch = 20)
points(x1, apply(E.bart, 1, max), col = bluegreen, add = TRUE, pch = 20)
points(x1, apply(E.race, 1, max), col = vermillion, add = TRUE, pch = 20)
abline(h = 0, lty = 3)
# No. Branches evidence margin plot
par(mar=c(4,0,1,1))
plot(apply(E.pitt, 2, max), x2, yaxt = "n", xlim = c(-25,140), xlab = "Evidence (db)", col = skyblue, pch = 20)
points(apply(E.mart, 2, max), x2, col = reddish, add = TRUE, pch = 20)
points(apply(E.bart, 2, max), x2, col = bluegreen, add = TRUE, pch = 20)
points(apply(E.race, 2, max), x2, col = vermillion, add = TRUE, pch = 20)
abline(v = 0, lty = 3)
# legend
par(mar = c(0,0,1,1))
plot(NA, xlim = c(-25,140), ylim = c(-25,140), axes = FALSE)
legend("topright", legend = c("R. bartlettii", "R. martinezii", "R. pittieri", "R. racemiflorum"), text.font = 3, fill = fig5cols[-1], bty = "n")
dev.off()
