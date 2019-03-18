###################################### Packages ######################################
if(!require("kableExtra")) install.packages("kableExtra")
require(kableExtra)
if(!require("qcc")) install.packages("qcc")
require(qcc)
if(!require("gamlss")) install.packages("gamlss")
require(gamlss)
if(!require("devtools")) install.packages("devtools")
require(devtools)
if(!require("DistMom")) install_github("jiperezga/DistMom")
require(DistMom)
if(!require("e1071")) install.packages("e1071")
require(e1071)
if(!require("truncgof")) install.packages("truncgof")
require(truncgof)
if(!require("evmix")) install.packages("evmix")
require(evmix)

###################################### Data analysis ######################################

################### Code 1: Load data set
### Loading data
hospita <- read.table("data/hospitalization.dat", header = T)
hospita$cost <- hospita$cost/1000000
head(hospita)
surgery <- read.table("data/surgery.dat", header = T)
surgery$cost <- surgery$cost/1000000
print(head(surgery))

################### Code 2: Pareto chart and plot hospitalization frequencies
### Pareto chart hospitalization
parhosp <- pareto.chart(table(hospita$year), 
                        main = "Hospitalization", las = 1)
kable(parhosp,
      caption = "Distribution records hospitalization per year
      \\label{tab:frech}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))
plot(parhosp)

################### Code 3: Pareto chart and plot general surgery frequencies
### Pareto chart general surgery
parsurg <- pareto.chart(table(surgery$year),
                        main = "General surgery", las = 1)
kable(parsurg,
      caption = "Distribution records general surgery per year
      \\label{tab:frecs}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

################### Code 4: Histogram, Box-plot, Density and Scatterplot hospitalization individual cost
par(mfrow=c(2,2))
### Histogram
hist(hospita$cost, freq = F, xlab = "Costs", 
     main = "(a) Histogram of \n hospitalization costs")
### Box-plot
boxplot(hospita$cost~factor(hospita$year), col = blues.colors(6), 
        main = "(a) Box-plot of \n hospitalization services",
        xlab = "Costs", ylab = "Years", horizontal = T, las=1)
### Density
plot(density(hospita$cost), lwd = 2,
     main = "(c) Density of \n hospitalization costs")
### Scatterplot
plot(hospita$year,hospita$cost, ylab="Costs", 
     xlab = "Year", main = "(d) Year vs Cost \n hospitalization")

################### Code 5: Histogram, Box-plot, Density and Scatterplot general surgery individual cost
par(mfrow=c(2,2))
### Histogram
hist(surgery$cost, freq = F, xlab = "Costs",
     main = "(a) Histogram of \n general surgery costs")
### Box-plot
boxplot(surgery$cost~factor(surgery$year), col = blues.colors(6), 
        main = "(b) Box-plot of \n hospitalization services", 
        xlab = "Costs", ylab = "Years", horizontal = T, las=1)
### Density
plot(density(surgery$cost), lwd = 2, 
     main = "(c) Density of \n general surgery costs")
### Scatterplot
plot(surgery$year,surgery$cost, ylab="Costs",
     xlab = "Year", main = "(d) Year vs Cost \n general surgery")

################### Code 6: Best fit with gamlss for frequencies of hospitalization services
### The adjustment is made
FitN_hosp1 <- fitDist(y = table(hospita$year), type = "counts")
### Estimation of second and third distribution with best fit
FitN_hosp2 <- gamlssML(formula = table(hospita$year), family = GPO)
FitN_hosp3 <- gamlssML(formula = table(hospita$year), family = NBI)
### The five distributions that present the best fit are
kable(rbind(FitN_hosp1$fits[1:5]),
      caption = "Better fit for frequencies of hospitalization services
      \\label{tab:fitfhosp}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

################### Code 7: Statistical measurements of hospitalization frequencies
### Estimation of mean, variance, skewness and excess kurtosis
#### Mean
MeanNhEmp <- mean(table(hospita$year))
MeanNhPIG <- moments(k = 1, dist = "PIG", domain = "counts",
                     param = c(mu = FitN_hosp1$mu, sigma = FitN_hosp1$sigma))
MeanNhGPO <- moments(k = 1, dist = "GPO", domain = "counts",
                     param = c(mu = FitN_hosp2$mu, sigma = FitN_hosp2$sigma))
MeanNhNBI <- moments(k = 1, dist = "NBI", domain = "counts",
                     param = c(mu = FitN_hosp3$mu, sigma = FitN_hosp3$sigma))
#### Variance
VariNhEmp <- var(table(hospita$year))
VariNhPIG <- moments(k = 2, dist = "PIG", domain = "counts",
                     param = c(mu = FitN_hosp1$mu, sigma = FitN_hosp1$sigma),
                     central = TRUE)
VariNhGPO <- moments(k = 2, dist = "GPO", domain = "counts",
                     param = c(mu = FitN_hosp2$mu, sigma = FitN_hosp2$sigma),
                     central = TRUE)
VariNhNBI <- moments(k = 2, dist = "NBI", domain = "counts",
                     param = c(mu = FitN_hosp3$mu, sigma = FitN_hosp3$sigma),
                     central = TRUE)
### Skewness
SkewNhEmp <- skewness(table(hospita$year), type = 1)
SkewNhPIG <- skew(dist = "PIG", domain = "counts", 
                  param = c(mu = FitN_hosp1$mu, sigma = FitN_hosp1$sigma))
SkewNhGPO <- skew(dist = "GPO", domain = "counts", 
                  param = c(mu = FitN_hosp2$mu, sigma = FitN_hosp2$sigma))
SkewNhNBI <- skew(dist = "NBI", domain = "counts", 
                  param = c(mu = FitN_hosp3$mu, sigma = FitN_hosp3$sigma))
### Excess Kurtosis
KurtNhEmp <- kurtosis(table(hospita$year), type = 1)
KurtNhPIG <- kurt(dist = "PIG", domain = "counts", excess = TRUE,
                  param = c(mu = FitN_hosp1$mu, sigma = FitN_hosp1$sigma))
KurtNhGPO <- kurt(dist = "GPO", domain = "counts", excess = TRUE,
                  param = c(mu = FitN_hosp2$mu, sigma = FitN_hosp2$sigma))
KurtNhNBI <- kurt(dist = "NBI", domain = "counts", excess = TRUE,
                  param = c(mu = FitN_hosp3$mu, sigma = FitN_hosp3$sigma))

kable(cbind(data.frame(Dist = c("Empirical", "PIG", "GPO", "BNI")),
            "Mean" = c(MeanNhEmp, unname(MeanNhPIG), unname(MeanNhGPO),
                       unname(MeanNhNBI)), "Variance" = c(VariNhEmp, 
                                                          unname(VariNhPIG), unname(VariNhGPO), unname(VariNhNBI)),
            "Skewness" = c(SkewNhEmp, unname(SkewNhPIG),
                           unname(SkewNhGPO), unname(SkewNhNBI)), "Excess kurtosis" =
              c(KurtNhEmp, unname(KurtNhPIG), unname(KurtNhGPO),
                unname(KurtNhNBI))), 
      caption = "Statistical measurements of hospitalization
      frequencies \\label{tab:StatisticsNh}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position")) 

################### Code 8:Adjustment cumulative frequencies of hospitalization
### Empirical vs Theorical cumulative distribution function
plot(ecdf(table(hospita$year)), lwd = 3,
     xlab = "Sample quantiles of frequencies of hospitalization", 
     main = "Adjustment cumulative frequencies of hospitalization")
curve(pPIG(x, mu = FitN_hosp1$mu, sigma = FitN_hosp1$sigma),
      add = T, from = 200, to = 600, lwd = 3, col = "blue", lty = 1)
curve(pGPO(x, mu = FitN_hosp2$mu, sigma = FitN_hosp2$sigma),
      add = T, from = 200, to = 600, lwd = 3, col = "green", lty = 2)
curve(pNBI(x, mu = FitN_hosp3$mu, sigma = FitN_hosp3$sigma),
      add = T, from = 200, to = 600, lwd = 3, col = "red", lty = 4)
grid()
legend("bottomright", lty = 1, col = c("black", "blue", "green",
                                       "red"), legend = c("Cumulative empirical distribution",
                                                          "Poisson-Inverse Gaussian", "Generalised Poisson",
                                                          "Negative Binomial type I"), lwd = 2)

################### Code 9: Goodness-of-fit tests for hospitalization frequencies
set.seed(1248) # A seed is established so that the results can be replicated
freqhosp <- table(hospita$year)
### Kolmogorov-Smirnov test
kshPIG <- ks.test(x = freqhosp, distn =  "pPIG", H = min(freqhosp),
                  fit = list(mu = FitN_hosp1$mu, sigma = FitN_hosp1$sigma))
kshGPO <- ks.test(x = freqhosp, distn =  "pGPO", H = min(freqhosp),
                  fit = list(mu = FitN_hosp2$mu, sigma = FitN_hosp2$sigma))
kshNBI <- ks.test(x = freqhosp, distn =  "pNBI", H = min(freqhosp),
                  fit = list(mu = FitN_hosp3$mu, sigma = FitN_hosp3$sigma))
### Cramer-von Mises test
cvmhPIG <- w2.test(x = freqhosp, distn =  "pPIG", H = min(freqhosp),
                   fit = list(mu = FitN_hosp1$mu, sigma = FitN_hosp1$sigma))
cvmhGPO <- w2.test(x = freqhosp, distn =  "pGPO", H = min(freqhosp),
                   fit = list(mu = FitN_hosp2$mu, sigma = FitN_hosp2$sigma))
cvmhNBI <- w2.test(x = freqhosp, distn =  "pNBI", H = min(freqhosp),
                   fit = list(mu = FitN_hosp3$mu, sigma = FitN_hosp3$sigma))
### Kuiper test
kuhPIG <- v.test(x = freqhosp, distn =  "pPIG", H = min(freqhosp),
                 fit = list(mu = FitN_hosp1$mu, sigma = FitN_hosp1$sigma))
kuhGPO <- v.test(x = freqhosp, distn =  "pGPO", H = min(freqhosp),
                 fit = list(mu = FitN_hosp2$mu, sigma = FitN_hosp2$sigma))
kuhNBI <- v.test(x = freqhosp, distn =  "pNBI", H = min(freqhosp),
                 fit = list(mu = FitN_hosp3$mu, sigma = FitN_hosp3$sigma))

kable(cbind(data.frame(Dist = c("PIG", "GPO", "BNI")),
            "ks.test" = c(kshPIG$p.value, kshGPO$p.value,
                          kshNBI$p.value), "w2.test" = c(cvmhPIG$p.value,
                                                         cvmhGPO$p.value, cvmhNBI$p.value), "v.test" = c(
                                                           kuhPIG$p.value, kuhGPO$p.value, kuhNBI$p.value)),
      caption = "Goodness-of-fit tests hospitalization frequencies
      \\label{tab:gftNh}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

###################################### Frequency model estimation for general surgery ######################################

################### Code 10: Best fit with gamlss for frequencies of general surgery services
### The adjustment is made
FitN_surg1 <- fitDist(y = table(surgery$year), type = "counts")
### Estimation of second and third distribution with best fit
FitN_surg2 <- gamlssML(formula = table(surgery$year), family = PIG)
FitN_surg3 <- gamlssML(formula = table(surgery$year), family = GPO)
### The five distributions that present the best fit are
kable(rbind(FitN_surg1$fits[1:5]),
      caption = "Better fit for frequencies of general surgery
      services \\label{tab:fitfhosp}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

################### Code 11: Statistical measurements of general surgery frequencies
### Estimation of mean, variance, skewness and excess kurtosis
#### Mean 
MeanNsEmp <- mean(table(surgery$year))
MeanNsDEL <- moments(k = 1, dist = "DEL", domain = "counts",
                     param = c(mu = FitN_surg1$mu, sigma = FitN_surg1$sigma,
                               nu = FitN_surg1$nu))
MeanNsPIG <- moments(k = 1, dist = "PIG", domain = "counts",
                     param = c(mu = FitN_surg2$mu, sigma = FitN_surg2$sigma))
MeanNsGPO <- moments(k = 1, dist = "GPO", domain = "counts",
                     param = c(mu = FitN_surg3$mu, sigma = FitN_surg3$sigma))
#### Variance
VariNsEmp <- var(table(surgery$year))
VariNsDEL <- moments(k = 2, dist = "DEL", domain = "counts",
                     param = c(mu = FitN_surg1$mu, sigma = FitN_surg1$sigma,
                               nu = FitN_surg1$nu), central = TRUE)
VariNsPIG <- moments(k = 2, dist = "PIG", domain = "counts",
                     param = c(mu = FitN_surg2$mu, sigma = FitN_surg2$sigma),
                     central = TRUE)
VariNsGPO <- moments(k = 2, dist = "GPO", domain = "counts",
                     param = c(mu = FitN_surg3$mu, sigma = FitN_surg3$sigma),
                     central = TRUE)
### Skewness
SkewNsEmp <- skewness(table(surgery$year), type = 1)
SkewNsDEL <- skew(dist = "DEL", domain = "counts", 
                  param = c(mu = FitN_surg1$mu, sigma = FitN_surg1$sigma,
                            nu = FitN_surg1$nu))
SkewNsPIG <- skew(dist = "PIG", domain = "counts", 
                  param = c(mu = FitN_surg2$mu, sigma = FitN_surg2$sigma))
SkewNsGPO <- skew(dist = "GPO", domain = "counts", 
                  param = c(mu = FitN_surg3$mu, sigma = FitN_surg3$sigma))
### Excess Kurtosis
KurtNsEmp <- kurtosis(table(surgery$year), type = 1)
KurtNsDEL <- kurt(dist = "DEL", domain = "counts", excess = TRUE,
                  param = c(mu = FitN_surg1$mu, sigma = FitN_surg1$sigma,
                            nu = FitN_surg1$nu))
KurtNsPIG <- kurt(dist = "PIG", domain = "counts", excess = TRUE,
                  param = c(mu = FitN_surg2$mu, sigma = FitN_surg2$sigma))
KurtNsGPO <- kurt(dist = "GPO", domain = "counts", excess = TRUE,
                  param = c(mu = FitN_surg3$mu, sigma = FitN_surg3$sigma))

kable(cbind(data.frame(Dist = c("Empirical", "DEL", "PIG", "BNI")),
            "Mean" = c(MeanNsEmp, unname(MeanNsDEL), unname(MeanNsPIG),
                       unname(MeanNsGPO)), "Variance" = c(VariNsEmp, 
                                                          unname(VariNsDEL), unname(VariNsPIG), unname(VariNsGPO)),
            "Skewness" = c(SkewNsEmp, unname(SkewNsDEL),
                           unname(SkewNsPIG), unname(SkewNsGPO)), "Excess kurtosis" =
              c(KurtNsEmp, unname(KurtNsDEL), unname(KurtNsPIG),
                unname(KurtNsGPO))), 
      caption = "Statistical measurements of general surgery
                 frequencies \\label{tab:StatisticsNs}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position")) 

################### Code 12: Adjustment cumulative frequencies of general surgery
plot(ecdf(table(surgery$year)), lwd = 3,
     xlab = "Sample quantiles of frequencies of general surgery", 
     main = "Adjustment cumulative frequencies of general surgery")
curve(pDEL(x, mu = FitN_surg1$mu, sigma = FitN_surg1$sigma,
           nu = FitN_surg1$nu), add = T, from = 10, to = 250, lwd = 3,
      col = "blue", lty = 1)
curve(pPIG(x, mu = FitN_surg2$mu, sigma = FitN_surg2$sigma),
      add = T, from = 10, to = 250, lwd = 3, col = "green", lty = 2)
curve(pGPO(x, mu = FitN_surg3$mu, sigma = FitN_surg3$sigma),
      add = T, from = 10, to = 250, lwd = 3, col = "red", lty = 4)
grid()
legend("bottomright", lty = 1, col = c("black", "blue", "green",
                                       "red"), legend = c("Cumulative empirical distribution",
                                                          "Poisson-Inverse Gaussian", "Generalised Poisson",
                                                          "Negative Binomial type I"), lwd = 2)

################### Code 13: Goodness-of-fit for tests general surgery frequencies
set.seed(1248) # A seed is established so that the results can be replicated
freqsurg <- table(surgery$year)
### Kolmogorov-Smirnov test
kshDEL <- ks.test(x = freqsurg, distn =  "pDEL", H = min(freqsurg),
                  fit = list(mu = FitN_surg1$mu, sigma = FitN_surg1$sigma,
                             nu = FitN_surg1$nu))
kshPIG <- ks.test(x = freqsurg, distn =  "pPIG", H = min(freqsurg),
                  fit = list(mu = FitN_surg2$mu, sigma = FitN_surg2$sigma))
kshGPO <- ks.test(x = freqsurg, distn =  "pGPO", H = min(freqsurg),
                  fit = list(mu = FitN_surg3$mu, sigma = FitN_surg3$sigma))
### Cramer-von Mises test
cvmhDEL <- w2.test(x = freqsurg, distn =  "pDEL", H = min(freqsurg),
                   fit = list(mu = FitN_surg1$mu, sigma = FitN_surg1$sigma,
                              nu = FitN_surg1$nu))
cvmhPIG <- w2.test(x = freqsurg, distn =  "pPIG", H = min(freqsurg),
                   fit = list(mu = FitN_surg2$mu, sigma = FitN_surg2$sigma))
cvmhGPO <- w2.test(x = freqsurg, distn =  "pGPO", H = min(freqsurg),
                   fit = list(mu = FitN_surg3$mu, sigma = FitN_surg3$sigma))
### Kuiper test
kuhDEL <- v.test(x = freqsurg, distn =  "pDEL", H = min(freqsurg),
                 fit = list(mu = FitN_surg1$mu, sigma = FitN_surg1$sigma,
                            nu = FitN_surg1$nu))
kuhPIG <- v.test(x = freqsurg, distn =  "pPIG", H = min(freqsurg),
                 fit = list(mu = FitN_surg2$mu, sigma = FitN_surg2$sigma))
kuhGPO <- v.test(x = freqsurg, distn =  "pGPO", H = min(freqsurg),
                 fit = list(mu = FitN_surg3$mu, sigma = FitN_surg3$sigma))

kable(cbind(data.frame(Dist = c("DEL", "PIG", "GPO")),
            "ks.test" = c(kshDEL$p.value, kshPIG$p.value,
                          kshGPO$p.value), "w2.test" = c(cvmhDEL$p.value,
                                                         cvmhPIG$p.value, cvmhGPO$p.value), "v.test" = c(
                                                           kuhDEL$p.value, kuhPIG$p.value, kuhGPO$p.value)),
      caption = "Goodness-of-fit tests general surgery frequencies
      \\label{tab:gftNs}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

###################################### Mean residual life analysis for hospitalization ######################################

################### Code 14: Mean residual life hospitalizacion service
### Graphic adjustment of mean residual life
mrlplot(hospita$cost, legend.loc = NULL, try.thresh = NULL)
abline(v = 35, col = "red", lwd = 2, lty = 2)
abline(v = 60, col = "blue", lwd = 2, lty = 2)

###################################### Mean residual life analysis for general surgery ######################################

################### Code 15: Mean residual life hospitalizacion service
### Graphic adjustment of mean residual life
mrlplot(surgery$cost, legend.loc = NULL, try.thresh = NULL)
abline(v = 16, col = "red", lwd = 2, lty = 2)

###################################### Tail index with Hill plot for hospitalization ######################################

################### Code 16: Hill, AltHill, SmooHill and AltSmooHill for hospitalizacion service
### Hill, AltHill, SmooHill and AltSmooHill
par(mfrow = c(2, 2))
hillplot(data = hospita$cost, alpha = 0.05, legend.loc = "top",
         hill.type = "Hill", ylim = c(0, 2))
title(main = "(a)                   ", line = 3.3)
hillplot(data = hospita$cost, alpha = 0.05, legend.loc = "top",
         hill.type = "Hill", ylim = c(0, 2), x.theta = TRUE)
title(main = "(b)                        ", line = 3.3)
hillplot(data = hospita$cost, alpha = 0.05, legend.loc = "top",
         hill.type = "SmooHill", ylim = c(0, 2))
title(main = "(c)                             ", line = 3.3)
hillplot(data = hospita$cost, alpha = 0.05, legend.loc = "top",
         hill.type = "SmooHill", ylim = c(0, 2), x.theta = TRUE)
title(main = "(d)                                 ", line = 3.3)

###################################### Tail index with Hill plot for general surgery ######################################

################### Code 17:
################### Code 18:
################### Code 19:
################### Code 20:
################### Code 21:
################### Code 22:
################### Code 23:
################### Code 24:
################### Code 25:
################### Code 26:
################### Code 27:
################### Code 28:
################### Code 29:
################### Code 30:
################### Code 31:
################### Code 32:
################### Code 33:
################### Code 34:
################### Code 35:
################### Code 36:
################### Code 37:
################### Code 38:
################### Code 39:
################### Code 40:
################### Code 41:
################### Code 42:
################### Code 43:
################### Code 44:
################### Code 45:
################### Code 46:
################### Code 47:
################### Code 48:
################### Code 49:
################### Code 50:
################### Code 51:
################### Code 52:
################### Code 53:
################### Code 54:
################### Code 55:
################### Code 56:
################### Code 57:
################### Code 58:
################### Code 59:
################### Code 60:
################### Code 61:
################### Code 62:
################### Code 63:
################### Code 64:
################### Code 65:
################### Code 66:
################### Code 67:
################### Code 68: