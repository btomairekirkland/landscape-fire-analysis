# Example GAM exploring the drivers of fire occurrence

# Clear data
rm(list = ls())

# Load libraries
library(mgcv)
library(corrplot)
library(RColorBrewer)
library(itsadug)

# Read in your gridded fire dataset
setwd("~/BTO projects/Polesia wildfires")
pix <- read.csv("pixel-df.csv", stringsAsFactors = T) 

# Create binary response variable
pix$burn <- 1
pix$burn[grepl("^CP", pix$z)] <- 0 ## Zero for non-fire observations

# Select numeric variables
num.cols <- unlist(lapply(pix, is.numeric))
num.df <- pix[,num.cols]
# Check correlation between variables
cor <- cor(num.df[,-4], method = "pearson")
options(max.print=1000000) 
round(cor,3)

# Calculate p-value
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# Matrix of the p-value of the correlation
p.mat <- cor.mtest(num.df[,-4])

tmp = cor # Copy matrix
tmp[ tmp > -0.6 & tmp < 0.6 ] = 0 # Select correlations greater than .6

# Plot correlogram
corrplot(tmp, type="upper", diag = F, # Add coefficient of correlation
         tl.col="black", tl.cex = .7, p.mat = p.mat, sig.level = .05, insig = "blank",
         col=brewer.pal(n = 8, name = "PuOr")) # Text label color and rotation

# Add your offset term 
ofs <- 5e-0 # Change this based on your sample of non-fire observations
pix$offset <- ifelse(pix$burn == 1, -log(1), -log(ofs)) # 1 for fire observations

# Run binomial GAM with non-collinear variables
## bam used for large datasets
## X and Y are centroid coordinates for voxels to account for spatial dependence
## n_year is an integer column representing the year numerically, starting from 1 
## day is an integer column representing the day of the year to account for temporal dependence
mod <- bam(burn ~ s(X, Y, bs = 'ds', k = 99) + s(n_year, bs = "cr", k = 19) + s(day, k = 100, bs = "cc") + 
             s(Agriculture, bs = 'cr', k = 20) + s(Raised.bog, bs = 'cr', k = 20) + s(Fen.and.transition.mire, bs = 'cr', k = 20) + 
             s(Scrub, bs = 'cr', k = 20) + s(Deciduous.forest, bs = 'cr', k = 20) + s(Coniferous.forest, bs = 'cr', k = 20) + 
             s(Meadows, bs = 'cr', k = 20) + s(Urban, bs = 'cr', k = 20) + s(friction, bs = 'cr', k = 20) + s(dist_road, bs = 'cr', k = 20) + 
             s(pop_dens, bs = 'cr', k = 20) + s(ffmc, bs = 'cr', k = 20) + s(rain, bs = 'cr', k = 20) + s(dc, bs = 'cr', k = 20) + 
             s(ndvi, bs = 'cr', k = 20) + s(fwi, bs = 'cr', k = 20) + s(max_temp, bs = 'cr', k = 20),
           offset = offset, data = pix, family = binomial, discrete = T, select = T)

# Check residuals and values of k
par(mfrow=c(2,2))
gam.check(mod)
dev.off()

# Explore model output 
summary(mod)

# Plot individual effects using temperature as an example
plot_smooth(mod, view = "Agriculture", 
            rm.ranef = F, transform = plogis, # Plot effects on response scale 
            ylim = c(0,1),
            xlab = "Proportion of agriculture", 
            ylab = "Probability of fire", 
            h0 = NULL)