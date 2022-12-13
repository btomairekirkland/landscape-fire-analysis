# Example GAM exploring the drivers of fire occurrence

# Clear data
rm(list = ls())

# Load libraries
library(mgcv)
library(interactions)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Read in your gridded fire data sample
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

# Add weight for each land cover type
# This is a dataframe with the weights to be used in the model, which reflects the proportion of each land cover type sampled 
## One row for each voxel with the weight assigned to the land cover type at the centre point of the voxel 
weights <- read.csv("model-weights.csv")
df <- inner_join(df, weights)

# Run model exploring burning of different key land cover types in relation to moisture 
mod <- bam(burn ~ s(soil_moist, bs = "cr", k = 20) + ti(Meadows, soil_moist, bs = c("cr","cr")) + ti(Meadows, ffmc, bs = c("cr","cr")) + 
             s(Raised.bog, bs = "cr", k = 20) + ti(Raised.bog, ffmc, bs = c("cr","cr")) + ti(Raised.bog, soil_moist, bs = c("cr","cr")) +
             ti(Deciduous.forest, soil_moist, bs = c("cr","cr")) + ti(Deciduous.forest, ffmc, bs = c("cr","cr")) + 
             s(Deciduous.forest, bs = "cr", k = 20) + s(Agriculture, bs = "cr", k = 20) + s(max_temp, bs = "cr") + s(rain, bs = "cr", k = 15) +
             ti(Fen.and.transition.mire, soil_moist, bs = c("cr","cr")) + s(Meadows, bs = "cr", k = 20) + 
             s(Fen.and.transition.mire, bs = "cr") + s(ffmc, bs = "cr", k = 20) + s(ndvi, bs = "cr", k = 20) + 
             protected + fire_season + s(X, Y, bs = "ds", k = 1000) + ti(Fen.and.transition.mire, ffmc, bs = c("cr","cr")),
           data = df[df$fire_season != "Non-fire season",], weights = weight, discrete = T, gamma = 1.4, select = T, family = binomial)

# Check model output
summary(mod)

# Check model residuals
par(mfrow=c(2,2))
gam.check(mod)
dev.off()

# Predict probability of burning and plot effects

# Example based on fen and transition mires and soil moisture 

# First create new dataframe for predictions
new.df <- expand.grid(Fen.and.transition.mire = seq(0,1,by=.001), 
                              soil_moist = c(23,57,82), ## Mean and quartile values 
                              protected = "yes",
                              fire_season = "Early",
                              Agriculture = mean(df$Agriculture),
                              ffmc = mean(df$ffmc),
                              Raised.bog = mean(df$Raised.bog),
                              Deciduous.forest = mean(df$Deciduous.forest),
                              ndvi = mean(df$ndvi),
                              max_temp = mean(df$max_temp), 
                              rain = mean(df$rain),
                              Meadows = mean(df$Meadows),
                              X = mean(df$X),
                              Y = mean(df$Y)) 

# Then predict to new dataframe and calculate confidence intervals
new.df[,15:16]<- predict(mod, newdata=new.df,  type = "response", se = TRUE)
# Convert standard error to confidence intervals
new.df$upr <- new.df$fit + (2 * new.df$se.fit)
new.df$lwr <- new.df$fit - (2 * new.df$se.fit)

# Produce plot
ggplot(new.df) +
  geom_line(aes(x = Fen.and.transition.mire, y = fit, col = as.factor(soil_moist))) +
  geom_ribbon(aes(x = Fen.and.transition.mire, ymin = lwr, ymax = upr, fill = as.factor(soil_moist)), col = NA,
              alpha = .3) +
  theme_classic() +
  scale_fill_manual(values = c('#a6bddb','#3690c0','#023858')) +
  scale_colour_manual(values = c('#a6bddb','#3690c0','#023858')) +
  xlab("Proportion of fen and transition mire") +
  ylab("Probability of burning") +
  theme(axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")))

