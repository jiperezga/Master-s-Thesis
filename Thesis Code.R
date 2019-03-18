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
if(!require("car")) install.packages("car")
require(car)

###################################### Data analysis ######################################

################### Code 1: Load data set
### Loading data
hospita <- read.table(file.choose(), header = T) # Search file "hospitalization.dat"
hospita$cost <- hospita$cost/1000000
head(hospita)
surgery <- read.table(file.choose(), header = T) # Search file "surgery.dat"
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
par(mfrow = c(1,1))
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

################### Code 17: Hill, AltHill, SmooHill and AltSmooHill for general surgery service
### Hill, AltHill, SmooHill and AltSmooHill
par(mfrow = c(2, 2))
hillplot(data = surgery$cost, alpha = 0.05, legend.loc = "top",
         hill.type = "Hill", ylim = c(0, 2))
title(main = "(a)                   ", line = 3.3)
hillplot(data = surgery$cost, alpha = 0.05, legend.loc = "top",
         hill.type = "Hill", ylim = c(0, 2), x.theta = TRUE)
title(main = "(b)                        ", line = 3.3)
hillplot(data = surgery$cost, alpha = 0.05, legend.loc = "top",
         hill.type = "SmooHill", ylim = c(0, 2), try.thresh = 
           quantile(surgery$cost, c(0.9, 0.95), na.rm = TRUE))
title(main = "(c)                             ", line = 3.3)
hillplot(data = surgery$cost, alpha = 0.05, legend.loc = "top",
         hill.type = "SmooHill", ylim = c(0, 2), x.theta = TRUE,
         try.thresh = quantile(surgery$cost, c(0.9, 0.95),
                               na.rm = TRUE))
title(main = "(d)                                 ", line = 3.3)

###################################### Adjustment with spliced distributions for hospitalization ######################################

################### Code 18: Adjustment with spliced distributions for hospitalization service
### The adjustment with evmix is made
SpHfit1 <- fgammagpd(x = hospita$cost)
SpHfit2 <- fnormgpd(x = hospita$cost)
SpHfit3 <- fweibullgpd(x = hospita$cost)

################### Code 19: Statistical measurements of spliced distributions for hospitalization services
### Estimation of mean, variance, skewness and excess kurtosis
#### Mean
MeanXhSpE <- mean(hospita$cost)
MeanXhSpG <- moments(k = 1, dist = "gammagpd", domain = "realplus", 
                     param = c(phiu = SpHfit1$phiu, gshape = SpHfit1$gshape,
                               gscale = SpHfit1$gscale, u = SpHfit1$u, xi = SpHfit1$xi,
                               sigmau = SpHfit1$sigmau))
MeanXhSpN <- moments(k = 1, dist = "normgpd", domain = "realline", 
                     param = c(phiu = SpHfit2$phiu, nmean = SpHfit2$nmean,
                               nsd = SpHfit2$nsd, u = SpHfit2$u, xi = SpHfit2$xi,
                               sigmau = SpHfit2$sigmau))
MeanXhSpW <- moments(k = 1, dist = "weibullgpd", domain = "realplus", 
                     param = c(phiu = SpHfit3$phiu, wshape = SpHfit3$wshape,
                               wscale = SpHfit3$wscale, u = SpHfit3$u, xi = SpHfit3$xi,
                               sigmau = SpHfit3$sigmau))
#### Variance
VariXhSpE <- var(hospita$cost)
VariXhSpG <- moments(k = 2, dist = "gammagpd", domain = "realplus", 
                     param = c(phiu = SpHfit1$phiu, gshape = SpHfit1$gshape,
                               gscale = SpHfit1$gscale, u = SpHfit1$u, xi = SpHfit1$xi,
                               sigmau = SpHfit1$sigmau), central = TRUE)
VariXhSpN <- moments(k = 2, dist = "normgpd", domain = "realline", 
                     param = c(phiu = SpHfit2$phiu, nmean = SpHfit2$nmean,
                               nsd = SpHfit2$nsd, u = SpHfit2$u, xi = SpHfit2$xi,
                               sigmau = SpHfit2$sigmau), central = TRUE)
VariXhSpW <- moments(k = 2, dist = "weibullgpd", domain = "realplus", 
                     param = c(phiu = SpHfit3$phiu, wshape = SpHfit3$wshape,
                               wscale = SpHfit3$wscale, u = SpHfit3$u, xi = SpHfit3$xi,
                               sigmau = SpHfit3$sigmau), central = TRUE)
### Skewness
SkewXhSpE <- skewness(hospita$cost, type = 1)
SkewXhSpG <- skew(dist = "gammagpd", domain = "realplus", 
                  param = c(phiu = SpHfit1$phiu, gshape = SpHfit1$gshape,
                            gscale = SpHfit1$gscale, u = SpHfit1$u, xi = SpHfit1$xi,
                            sigmau = SpHfit1$sigmau))
SkewXhSpN <- skew(dist = "normgpd", domain = "realline", 
                  param = c(phiu = SpHfit2$phiu, nmean = SpHfit2$nmean,
                            nsd = SpHfit2$nsd, u = SpHfit2$u, xi = SpHfit2$xi,
                            sigmau = SpHfit2$sigmau))
SkewXhSpW <- skew(dist = "weibullgpd", domain = "realplus", 
                  param = c(phiu = SpHfit3$phiu, wshape = SpHfit3$wshape,
                            wscale = SpHfit3$wscale, u = SpHfit3$u, xi = SpHfit3$xi,
                            sigmau = SpHfit3$sigmau))
### Excess Kurtosis
KurtXhSpE <- kurtosis(hospita$cost, type = 1)
KurtXhSpG <- kurt(dist = "gammagpd", domain = "realplus", 
                  param = c(phiu = SpHfit1$phiu, gshape = SpHfit1$gshape,
                            gscale = SpHfit1$gscale, u = SpHfit1$u, xi = SpHfit1$xi,
                            sigmau = SpHfit1$sigmau), excess = TRUE)
KurtXhSpN <- kurt(dist = "normgpd", domain = "realline", 
                  param = c(phiu = SpHfit2$phiu, nmean = SpHfit2$nmean,
                            nsd = SpHfit2$nsd, u = SpHfit2$u, xi = SpHfit2$xi,
                            sigmau = SpHfit2$sigmau), excess = TRUE)
KurtXhSpW <- kurt(dist = "weibullgpd", domain = "realplus", 
                  param = c(phiu = SpHfit3$phiu, wshape = SpHfit3$wshape,
                            wscale = SpHfit3$wscale, u = SpHfit3$u, xi = SpHfit3$xi,
                            sigmau = SpHfit3$sigmau), excess = TRUE)

kable(cbind(data.frame(Dist = c("Empirical", "gammagdp", "normgdp",
                                "weibullgdp")), "Mean" = c(MeanXhSpE, unname(MeanXhSpG),
                                                           unname(MeanXhSpN), unname(MeanXhSpW)), "Variance" = c(
                                                             VariXhSpE, unname(VariXhSpG), unname(VariXhSpN),
                                                             unname(VariXhSpW)), "Skewness" = c(SkewXhSpE,
                                                                                                unname(SkewXhSpG), unname(SkewXhSpN), unname(SkewXhSpW)),
            "Excess Kurtosis" = c(KurtXhSpE, unname(KurtXhSpG),
                                  unname(KurtXhSpN), unname(KurtXhSpW))), 
      caption = "Statistical measurements of spliced distributions for
      hospitalization services \\label{tab:StatisticsSpXh}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position")) 

################### Code 20: Adjustment of cumulative individual costs of hospitalization services with spliced distributions
FnXh <- ecdf(hospita$cost)
sortXh <- sort(hospita$cost)
### Empirical vs Theorical cumulative distribution function
par(mfrow = c(1,1))
plot(FnXh, lwd = 3,
     xlab = "Sample quantiles of individual costs of hospitalization", 
     main = "Adjustment cumulative individual costs of hospitalization")
fitXhG <- pgammagpd(q = sortXh, phiu = SpHfit1$phiu,
                    gshape = SpHfit1$gshape, gscale = SpHfit1$gscale,
                    u = SpHfit1$u, xi = SpHfit1$xi, sigmau = SpHfit1$sigmau)
lines(sortXh, fitXhG, lwd = 3, lty = 1, col = "blue")
fitXhN <- pnormgpd(q = sortXh, phiu = SpHfit2$phiu,
                   nmean = SpHfit2$nmean, nsd = SpHfit2$nsd, u = SpHfit2$u,
                   xi = SpHfit2$xi, sigmau = SpHfit2$sigmau)
lines(sortXh, fitXhN, lwd = 3, lty = 2, col = "red")
fitXhW <- pweibullgpd(q = sortXh, phiu = SpHfit3$phiu,
                      wshape = SpHfit3$wshape, wscale = SpHfit3$wscale,
                      u = SpHfit3$u, xi = SpHfit3$xi, sigmau = SpHfit3$sigmau)
lines(sortXh, fitXhW, lwd = 3, lty = 4, col = "green")
grid()
legend("bottomright", lty = 1, col = c("black", "blue", "red",
                                       "green"), legend = c("Cumulative empirical distribution",
                                                            "Gamma-Generalized Pareto", "Normal-Generalized Pareto",
                                                            "Weibull-Generalized Pareto"), lwd = 2)


################### Code 21: Adjustment of log-survival distribution of hospitalization with spliced istributions
### Empirical vs Theorical log-survival distribution function
survXh <- 1 - FnXh(sortXh)
plot(x = log(sortXh), y = log(survXh), lwd = 3,
     xlab = "log(Sample quantiles of individual cost of
     hospitalization)", ylab = "log(1 - Fn(x))", 
     main = "Adjustment of log-survival distribution of
     hospitalization", type = "l")
survXhG <- 1 - fitXhG
lines(log(sortXh), log(survXhG), lwd = 3, col = "blue")
survXhN <- 1 - fitXhN
lines(log(sortXh), log(survXhN), lwd = 3, col = "red", lty = 2)
survXhW <- 1 - fitXhW
lines(log(sortXh), log(survXhW), lwd = 3, col = "green", lty = 4)
grid()
legend("bottomleft", lty = 1, col = c("black", "blue", "red", "green"),
       legend = c("log-Survival Distribution",
                  "Gamma-Generalized Pareto", "Normal-Generalized Pareto",
                  "Weibull-Generalized Pareto"), lwd = 2)

################### Code 22: Q-Q plot spliced distribution for hospitalization
par(mai=rep(0.5, 4))
layout(matrix(c(1,1, 2,2, 0, 3,3, 0), ncol = 4, byrow = TRUE))
### QQ-plot 
qqPlot(x = hospita$cost, lwd = 1, distribution = "gammagpd",
       phiu = SpHfit1$phiu, gshape = SpHfit1$gshape,
       gscale = SpHfit1$gscale, u = SpHfit1$u, xi = SpHfit1$xi,
       sigmau = SpHfit1$sigmau, cex = 1, col.lines = "red" ,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(a) Gamma-Generalized Pareto Q-Q Plot", id = FALSE)
qqPlot(x = hospita$cost, lwd = 1, distribution = "normgpd",
       phiu = SpHfit2$phiu, nmean = SpHfit2$nmean,
       nsd = SpHfit2$nsd, u = SpHfit2$u, xi = SpHfit2$xi,
       sigmau = SpHfit2$sigmau, cex = 1, col.lines = "red" ,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(b) Normal-Generalized Pareto Q-Q Plot", id = FALSE)
qqPlot(x = hospita$cost, lwd = 1, distribution = "weibullgpd",
       phiu = SpHfit3$phi, wshape = SpHfit3$wshape,
       wscale = SpHfit3$wscale, u = SpHfit3$u, xi = SpHfit3$xi,
       sigmau = SpHfit3$sigmau, cex = 1, col.lines = "red" ,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(c) Weibull-Generalized Pareto Q-Q Plot", id = FALSE)


################### Code 23: Goodness-of-fit tests for hospitalization services for spliced distributions
### Kolmogorov-Smirnov test
set.seed(1248) # A seed is established so that the results can be replicated
kolmSpXhG <- ks.test(x = hospita$cost, distn =  "pgammagpd", 
                     H = min(hospita$cost), fit = list(phiu = SpHfit1$phiu,
                                                       gshape = SpHfit1$gshape, gscale = SpHfit1$gscale,
                                                       u = SpHfit1$u, xi = SpHfit1$xi, sigmau = SpHfit1$sigmau))
kolmSpXhN <- ks.test(x = hospita$cost, distn =  "pnormgpd", 
                     H = min(hospita$cost), fit = list(phiu = SpHfit2$phiu,
                                                       nmean = SpHfit2$nmean, nsd = SpHfit2$nsd, u = SpHfit2$u,
                                                       xi = SpHfit2$xi, sigmau = SpHfit2$sigmau))
kolmSpXhW <- ks.test(x = hospita$cost, distn =  "pweibullgpd", 
                     H = min(hospita$cost), fit = list(phiu = SpHfit3$phiu,
                                                       wshape = SpHfit3$wshape, wscale = SpHfit3$wscale,
                                                       u = SpHfit3$u, xi = SpHfit3$xi, sigmau = SpHfit3$sigmau))
### Cramer-von Mises test
cramSpXhG <- w2.test(x = hospita$cost, distn =  "pgammagpd", 
                     H = min(hospita$cost), fit = list(phiu = SpHfit1$phiu,
                                                       gshape = SpHfit1$gshape, gscale = SpHfit1$gscale,
                                                       u = SpHfit1$u, xi = SpHfit1$xi, sigmau = SpHfit1$sigmau))
cramSpXhN <- w2.test(x = hospita$cost, distn =  "pnormgpd", 
                     H = min(hospita$cost), fit = list(phiu = SpHfit2$phiu,
                                                       nmean = SpHfit2$nmean, nsd = SpHfit2$nsd, u = SpHfit2$u,
                                                       xi = SpHfit2$xi, sigmau = SpHfit2$sigmau))
cramSpXhW <- w2.test(x = hospita$cost, distn =  "pweibullgpd", 
                     H = min(hospita$cost), fit = list(phiu = SpHfit3$phiu,
                                                       wshape = SpHfit3$wshape, wscale = SpHfit3$wscale,
                                                       u = SpHfit3$u, xi = SpHfit3$xi, sigmau = SpHfit3$sigmau))
### Kuiper test
kuipSpXhG <- v.test(x = hospita$cost, distn =  "pgammagpd", 
                    H = min(hospita$cost), fit = list(phiu = SpHfit1$phiu,
                                                      gshape = SpHfit1$gshape, gscale = SpHfit1$gscale,
                                                      u = SpHfit1$u, xi = SpHfit1$xi, sigmau = SpHfit1$sigmau))
kuipSpXhN <- v.test(x = hospita$cost, distn =  "pnormgpd", 
                    H = min(hospita$cost), fit = list(phiu = SpHfit2$phiu,
                                                      nmean = SpHfit2$nmean, nsd = SpHfit2$nsd, u = SpHfit2$u,
                                                      xi = SpHfit2$xi, sigmau = SpHfit2$sigmau))
kuipSpXhW <- v.test(x = hospita$cost, distn =  "pweibullgpd", 
                    H = min(hospita$cost), fit = list(phiu = SpHfit3$phiu,
                                                      wshape = SpHfit3$wshape, wscale = SpHfit3$wscale,
                                                      u = SpHfit3$u, xi = SpHfit3$xi, sigmau = SpHfit3$sigmau))
### Supremum class Upper Tail Anderson-Darling test
adupSpXhG <- adup.test(x = hospita$cost, distn =  "pgammagpd", 
                       H = min(hospita$cost), fit = list(phiu = SpHfit1$phiu,
                                                         gshape = SpHfit1$gshape, gscale = SpHfit1$gscale,
                                                         u = SpHfit1$u, xi = SpHfit1$xi, sigmau = SpHfit1$sigmau))
adupSpXhN <- adup.test(x = hospita$cost, distn =  "pnormgpd", 
                       H = min(hospita$cost), fit = list(phiu = SpHfit2$phiu,
                                                         nmean = SpHfit2$nmean, nsd = SpHfit2$nsd, u = SpHfit2$u,
                                                         xi = SpHfit2$xi, sigmau = SpHfit2$sigmau))
adupSpXhW <- adup.test(x = hospita$cost, distn =  "pweibullgpd", 
                       H = min(hospita$cost), fit = list(phiu = SpHfit3$phiu,
                                                         wshape = SpHfit3$wshape, wscale = SpHfit3$wscale,
                                                         u = SpHfit3$u, xi = SpHfit3$xi, sigmau = SpHfit3$sigmau))
### Quadratic Class Upper Tail Anderson-Darling test
ad2upSpXhG <- ad2up.test(x = hospita$cost, distn =  "pgammagpd", 
                         H = min(hospita$cost), fit = list(phiu = SpHfit1$phiu,
                                                           gshape = SpHfit1$gshape, gscale = SpHfit1$gscale,
                                                           u = SpHfit1$u, xi = SpHfit1$xi, sigmau = SpHfit1$sigmau))
ad2upSpXhN <- ad2up.test(x = hospita$cost, distn =  "pnormgpd", 
                         H = min(hospita$cost), fit = list(phiu = SpHfit2$phiu,
                                                           nmean = SpHfit2$nmean, nsd = SpHfit2$nsd, u = SpHfit2$u,
                                                           xi = SpHfit2$xi, sigmau = SpHfit2$sigmau))
ad2upSpXhW <- ad2up.test(x = hospita$cost, distn =  "pweibullgpd", 
                         H = min(hospita$cost), fit = list(phiu = SpHfit3$phiu,
                                                           wshape = SpHfit3$wshape, wscale = SpHfit3$wscale,
                                                           u = SpHfit3$u, xi = SpHfit3$xi, sigmau = SpHfit3$sigmau))
### Results table
kable(cbind(data.frame(Dist = c("gammagpd", "normgpd", "weibullgpd")),
            "ks.test" = c(kolmSpXhG$p.value, kolmSpXhN$p.value,
                          kolmSpXhW$p.value), "w2.test" = c(cramSpXhG$p.value,
                                                            cramSpXhN$p.value, cramSpXhW$p.value), "v.test" = c(
                                                              kuipSpXhG$p.value, kuipSpXhN$p.value, kuipSpXhW$p.value),
            "adup.test" = c(adupSpXhG$p.value, adupSpXhN$p.value,
                            adupSpXhW$p.value), "ad2up.test" = c(ad2upSpXhG$p.value,
                                                                 ad2upSpXhN$p.value, ad2upSpXhW$p.value)),
      caption = "Goodness-of-fit tests hospitalization services for
      spliced distributions \\label{tab:gftSpXh}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))


###################################### Adjustment with spliced distributions for general surgery ######################################

################### Code 24: Adjustment with spliced distributions for general surgery service
### The adjustment with evmix is made
SpSfit1 <- fgammagpd(x = surgery$cost)
SpSfit2 <- fnormgpd(x = surgery$cost)
SpSfit3 <- fweibullgpd(x = surgery$cost)

################### Code 25: Statistical measurements of spliced distributions for general surgery services
### Estimation of mean, variance, skewness and excess kurtosis
#### Mean
MeanXsSpE <- mean(surgery$cost)
MeanXsSpG <- moments(k = 1, dist = "gammagpd", domain = "realplus", 
                     param = c(phiu = SpSfit1$phiu, gshape = SpSfit1$gshape,
                               gscale = SpSfit1$gscale, u = SpSfit1$u, xi = SpSfit1$xi,
                               sigmau = SpSfit1$sigmau))
MeanXsSpN <- moments(k = 1, dist = "normgpd", domain = "realline", 
                     param = c(phiu = SpSfit2$phiu, nmean = SpSfit2$nmean,
                               nsd = SpSfit2$nsd, u = SpSfit2$u, xi = SpSfit2$xi,
                               sigmau = SpSfit2$sigmau))
MeanXsSpW <- moments(k = 1, dist = "weibullgpd", domain = "realplus", 
                     param = c(phiu = SpSfit3$phiu, wshape = SpSfit3$wshape,
                               wscale = SpSfit3$wscale, u = SpSfit3$u, xi = SpSfit3$xi,
                               sigmau = SpSfit3$sigmau))
#### Variance
VariXsSpE <- var(surgery$cost)
VariXsSpG <- moments(k = 2, dist = "gammagpd", domain = "realplus", 
                     param = c(phiu = SpSfit1$phiu, gshape = SpSfit1$gshape,
                               gscale = SpSfit1$gscale, u = SpSfit1$u, xi = SpSfit1$xi,
                               sigmau = SpSfit1$sigmau), central = TRUE)
VariXsSpN <- moments(k = 2, dist = "normgpd", domain = "realline", 
                     param = c(phiu = SpSfit2$phiu, nmean = SpSfit2$nmean,
                               nsd = SpSfit2$nsd, u = SpSfit2$u, xi = SpSfit2$xi,
                               sigmau = SpSfit2$sigmau), central = TRUE)
VariXsSpW <- moments(k = 2, dist = "weibullgpd", domain = "realplus", 
                     param = c(phiu = SpSfit3$phiu, wshape = SpSfit3$wshape,
                               wscale = SpSfit3$wscale, u = SpSfit3$u, xi = SpSfit3$xi,
                               sigmau = SpSfit3$sigmau), central = TRUE)
### Skewness
SkewXsSpE <- skewness(surgery$cost, type = 1)
SkewXsSpG <- skew(dist = "gammagpd", domain = "realplus", 
                  param = c(phiu = SpSfit1$phiu, gshape = SpSfit1$gshape,
                            gscale = SpSfit1$gscale, u = SpSfit1$u, xi = SpSfit1$xi,
                            sigmau = SpSfit1$sigmau))
SkewXsSpN <- skew(dist = "normgpd", domain = "realline", 
                  param = c(phiu = SpSfit2$phiu, nmean = SpSfit2$nmean,
                            nsd = SpSfit2$nsd, u = SpSfit2$u, xi = SpSfit2$xi,
                            sigmau = SpSfit2$sigmau))
SkewXsSpW <- skew(dist = "weibullgpd", domain = "realplus", 
                  param = c(phiu = SpSfit3$phiu, wshape = SpSfit3$wshape,
                            wscale = SpSfit3$wscale, u = SpSfit3$u, xi = SpSfit3$xi,
                            sigmau = SpSfit3$sigmau))
### Excess Kurtosis
KurtXsSpE <- kurtosis(surgery$cost, type = 1)
KurtXsSpG <- kurt(dist = "gammagpd", domain = "realplus", 
                  param = c(phiu = SpSfit1$phiu, gshape = SpSfit1$gshape,
                            gscale = SpSfit1$gscale, u = SpSfit1$u, xi = SpSfit1$xi,
                            sigmau = SpSfit1$sigmau), excess = TRUE)
KurtXsSpN <- kurt(dist = "normgpd", domain = "realline", 
                  param = c(phiu = SpSfit2$phiu, nmean = SpSfit2$nmean,
                            nsd = SpSfit2$nsd, u = SpSfit2$u, xi = SpSfit2$xi,
                            sigmau = SpSfit2$sigmau), excess = TRUE)
KurtXsSpW <- kurt(dist = "weibullgpd", domain = "realplus", 
                  param = c(phiu = SpSfit3$phiu, wshape = SpSfit3$wshape,
                            wscale = SpSfit3$wscale, u = SpSfit3$u, xi = SpSfit3$xi,
                            sigmau = SpSfit3$sigmau), excess = TRUE)

kable(cbind(data.frame(Dist = c("Empirical", "gammagdp", "normgdp",
                                "weibullgdp")), "Mean" = c(MeanXsSpE, unname(MeanXsSpG),
                                                           unname(MeanXsSpN), unname(MeanXsSpW)), "Variance" = c(
                                                             VariXsSpE, "does not exist", unname(VariXsSpN),
                                                             "does not exist"), "Skewness" = c(SkewXsSpE,
                                                                                               "does not exist", unname(SkewXsSpN), "does not exist"),
            "Excess Kurtosis" = c(KurtXsSpE, "does not exist",
                                  unname(KurtXsSpN), "does not exist")), 
      caption = "Statistical measurements of spliced distributions for
      general surgery services \\label{tab:StatisticsSpXs}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))  

################### Code 26: Adjustment of cumulative individual costs of general surgery services with spliced distributions
FnXs <- ecdf(surgery$cost)
sortXs <- sort(surgery$cost)
### Empirical vs Theorical cumulative distribution function
par(mfrow = c(1,1))
plot(FnXs, lwd = 3,
     xlab = "Sample quantiles of individual costs of general surgery", 
     main = "Adjustment cumulative individual costs of general surgery")
fitXsG <- pgammagpd(q = sortXs, phiu = SpSfit1$phiu,
                    gshape = SpSfit1$gshape, gscale = SpSfit1$gscale,
                    u = SpSfit1$u, xi = SpSfit1$xi, sigmau = SpSfit1$sigmau)
lines(sortXs, fitXsG, lwd = 3, lty = 1, col = "blue")
fitXsN <- pnormgpd(q = sortXs, phiu = SpSfit2$phiu,
                   nmean = SpSfit2$nmean, nsd = SpSfit2$nsd, u = SpSfit2$u,
                   xi = SpSfit2$xi, sigmau = SpSfit2$sigmau)
lines(sortXs, fitXsN, lwd = 3, lty = 2, col = "red")
fitXsW <- pweibullgpd(q = sortXs, phiu = SpSfit3$phiu,
                      wshape = SpSfit3$wshape, wscale = SpSfit3$wscale,
                      u = SpSfit3$u, xi = SpSfit3$xi, sigmau = SpSfit3$sigmau)
lines(sortXs, fitXsW, lwd = 3, lty = 4, col = "green")
grid()
legend("bottomright", lty = 1, col = c("black", "blue", "red",
                                       "green"), legend = c("Cumulative empirical distribution",
                                                            "Gamma-Generalized Pareto", "Normal-Generalized Pareto",
                                                            "Weibull-Generalized Pareto"), lwd = 2)

################### Code 27: Adjustment of log-survival distribution of general surgery with spliced distributions
### Empirical vs Theorical log-survival distribution function
survXs <- 1 - FnXs(sortXs)
plot(x = log(sortXs), y = log(survXs), lwd = 3,
     xlab = "log(Sample quantiles of individual cost of
     general surgery)", ylab = "log(1 - Fn(x))", 
     main = "Adjustment of log-survival distribution of
     general surgery", type = "l")
survXsG <- 1 - fitXsG
lines(log(sortXs), log(survXsG), lwd = 3, col = "blue")
survXsN <- 1 - fitXsN
lines(log(sortXs), log(survXsN), lwd = 3, col = "red", lty = 2)
survXsW <- 1 - fitXsW
lines(log(sortXs), log(survXsW), lwd = 3, col = "green", lty = 4)
grid()
legend("bottomleft", lty = 1, col = c("black", "blue", "red", "green"),
       legend = c("log-Survival Distribution",
                  "Gamma-Generalized Pareto", "Normal-Generalized Pareto",
                  "Weibull-Generalized Pareto"), lwd = 2)

################### Code 28: AQ-Q plot spliced distribution for general surgery
par(mai=rep(0.5, 4))
layout(matrix(c(1,1, 2,2, 0, 3,3, 0), ncol = 4, byrow = TRUE))
### QQ-plot 
qqPlot(x = surgery$cost, lwd = 1, distribution = "gammagpd",
       phiu = SpSfit1$phiu, gshape = SpSfit1$gshape,
       gscale = SpSfit1$gscale, u = SpSfit1$u, xi = SpSfit1$xi,
       sigmau = SpSfit1$sigmau, cex = 1, col.lines = "red" ,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(a) Gamma-Generalized Pareto Q-Q Plot", id = FALSE)
qqPlot(x = surgery$cost, lwd = 1, distribution = "normgpd",
       phiu = SpSfit2$phiu, nmean = SpSfit2$nmean,
       nsd = SpSfit2$nsd, u = SpSfit2$u, xi = SpSfit2$xi,
       sigmau = SpSfit2$sigmau, cex = 1, col.lines = "red" ,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(b) Normal-Generalized Pareto Q-Q Plot", id = FALSE)
qqPlot(x = surgery$cost, lwd = 1, distribution = "weibullgpd",
       phiu = SpSfit3$phi, wshape = SpSfit3$wshape,
       wscale = SpSfit3$wscale, u = SpSfit3$u, xi = SpSfit3$xi,
       sigmau = SpSfit3$sigmau, cex = 1, col.lines = "red" ,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(c) Weibull-Generalized Pareto Q-Q Plot", id = FALSE)

################### Code 29: Goodness-of-fit tests for general surgery services for spliced distributions
set.seed(1248) # A seed is established so that the results can be replicated
### Kolmogorov-Smirnov test
kolmSpXsG <- ks.test(x = surgery$cost, distn =  "pgammagpd", 
                     H = min(surgery$cost), fit = list(phiu = SpSfit1$phiu,
                                                       gshape = SpSfit1$gshape, gscale = SpSfit1$gscale,
                                                       u = SpSfit1$u, xi = SpSfit1$xi, sigmau = SpSfit1$sigmau))
kolmSpXsN <- ks.test(x = surgery$cost, distn =  "pnormgpd", 
                     H = min(surgery$cost), fit = list(phiu = SpSfit2$phiu,
                                                       nmean = SpSfit2$nmean, nsd = SpSfit2$nsd, u = SpSfit2$u,
                                                       xi = SpSfit2$xi, sigmau = SpSfit2$sigmau))
kolmSpXsW <- ks.test(x = surgery$cost, distn =  "pweibullgpd", 
                     H = min(surgery$cost), fit = list(phiu = SpSfit3$phiu,
                                                       wshape = SpSfit3$wshape, wscale = SpSfit3$wscale,
                                                       u = SpSfit3$u, xi = SpSfit3$xi, sigmau = SpSfit3$sigmau))
### Cramer-von Mises test
cramSpXsG <- w2.test(x = surgery$cost, distn =  "pgammagpd", 
                     H = min(surgery$cost), fit = list(phiu = SpSfit1$phiu,
                                                       gshape = SpSfit1$gshape, gscale = SpSfit1$gscale,
                                                       u = SpSfit1$u, xi = SpSfit1$xi, sigmau = SpSfit1$sigmau))
cramSpXsN <- w2.test(x = surgery$cost, distn =  "pnormgpd", 
                     H = min(surgery$cost), fit = list(phiu = SpSfit2$phiu,
                                                       nmean = SpSfit2$nmean, nsd = SpSfit2$nsd, u = SpSfit2$u,
                                                       xi = SpSfit2$xi, sigmau = SpSfit2$sigmau))
cramSpXsW <- w2.test(x = surgery$cost, distn =  "pweibullgpd", 
                     H = min(surgery$cost), fit = list(phiu = SpSfit3$phiu,
                                                       wshape = SpSfit3$wshape, wscale = SpSfit3$wscale,
                                                       u = SpSfit3$u, xi = SpSfit3$xi, sigmau = SpSfit3$sigmau))
### Kuiper test
kuipSpXsG <- v.test(x = surgery$cost, distn =  "pgammagpd", 
                    H = min(surgery$cost), fit = list(phiu = SpSfit1$phiu,
                                                      gshape = SpSfit1$gshape, gscale = SpSfit1$gscale,
                                                      u = SpSfit1$u, xi = SpSfit1$xi, sigmau = SpSfit1$sigmau))
kuipSpXsN <- v.test(x = surgery$cost, distn =  "pnormgpd", 
                    H = min(surgery$cost), fit = list(phiu = SpSfit2$phiu,
                                                      nmean = SpSfit2$nmean, nsd = SpSfit2$nsd, u = SpSfit2$u,
                                                      xi = SpSfit2$xi, sigmau = SpSfit2$sigmau))
kuipSpXsW <- v.test(x = surgery$cost, distn =  "pweibullgpd", 
                    H = min(surgery$cost), fit = list(phiu = SpSfit3$phiu,
                                                      wshape = SpSfit3$wshape, wscale = SpSfit3$wscale,
                                                      u = SpSfit3$u, xi = SpSfit3$xi, sigmau = SpSfit3$sigmau))
### Supremum class Upper Tail Anderson-Darling test
adupSpXsG <- adup.test(x = surgery$cost, distn =  "pgammagpd", 
                       H = min(surgery$cost), fit = list(phiu = SpSfit1$phiu,
                                                         gshape = SpSfit1$gshape, gscale = SpSfit1$gscale,
                                                         u = SpSfit1$u, xi = SpSfit1$xi, sigmau = SpSfit1$sigmau))
adupSpXsN <- adup.test(x = surgery$cost, distn =  "pnormgpd", 
                       H = min(surgery$cost), fit = list(phiu = SpSfit2$phiu,
                                                         nmean = SpSfit2$nmean, nsd = SpSfit2$nsd, u = SpSfit2$u,
                                                         xi = SpSfit2$xi, sigmau = SpSfit2$sigmau))
adupSpXsW <- adup.test(x = surgery$cost, distn =  "pweibullgpd", 
                       H = min(surgery$cost), fit = list(phiu = SpSfit3$phiu,
                                                         wshape = SpSfit3$wshape, wscale = SpSfit3$wscale,
                                                         u = SpSfit3$u, xi = SpSfit3$xi, sigmau = SpSfit3$sigmau))
### Quadratic Class Upper Tail Anderson-Darling test
ad2upSpXsG <- ad2up.test(x = surgery$cost, distn =  "pgammagpd", 
                         H = min(surgery$cost), fit = list(phiu = SpSfit1$phiu,
                                                           gshape = SpSfit1$gshape, gscale = SpSfit1$gscale,
                                                           u = SpSfit1$u, xi = SpSfit1$xi, sigmau = SpSfit1$sigmau))
ad2upSpXsN <- ad2up.test(x = surgery$cost, distn =  "pnormgpd", 
                         H = min(surgery$cost), fit = list(phiu = SpSfit2$phiu,
                                                           nmean = SpSfit2$nmean, nsd = SpSfit2$nsd, u = SpSfit2$u,
                                                           xi = SpSfit2$xi, sigmau = SpSfit2$sigmau))
ad2upSpXsW <- ad2up.test(x = surgery$cost, distn =  "pweibullgpd", 
                         H = min(surgery$cost), fit = list(phiu = SpSfit3$phiu,
                                                           wshape = SpSfit3$wshape, wscale = SpSfit3$wscale,
                                                           u = SpSfit3$u, xi = SpSfit3$xi, sigmau = SpSfit3$sigmau))
### Results table
kable(cbind(data.frame(Dist = c("gammagpd", "normgpd", "weibullgpd")),
            "ks.test" = c(kolmSpXsG$p.value, kolmSpXsN$p.value,
                          kolmSpXsW$p.value), "w2.test" = c(cramSpXsG$p.value,
                                                            cramSpXsN$p.value, cramSpXsW$p.value), "v.test" = c(
                                                              kuipSpXsG$p.value, kuipSpXsN$p.value, kuipSpXsW$p.value),
            "adup.test" = c(adupSpXsG$p.value, adupSpXsN$p.value,
                            adupSpXsW$p.value), "ad2up.test" = c(ad2upSpXsG$p.value,
                                                                 ad2upSpXsN$p.value, ad2upSpXsW$p.value)),
      caption = "Goodness-of-fit tests general surgery services for
      spliced distributions \\label{tab:gftSpXh}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

###################################### Risk measures estimation for hospitalization services ######################################

################### Code 30: Estimation of the tail index for hospitalization services
### Density function
fhosp <- function(x) dweibullgpd(x, wshape = SpHfit3$wshape, 
                                 wscale = SpHfit3$wscale, sigmau = SpHfit3$sigmau,
                                 u = SpHfit3$u, xi = SpHfit3$xi, phiu = SpHfit3$phiu)
### Cumulative function
Fhosp <- function(x) pweibullgpd(x, wshape = SpHfit3$wshape,
                                 wscale = SpHfit3$wscale, sigmau = SpHfit3$sigmau,
                                 u = SpHfit3$u, xi = SpHfit3$xi, phiu = SpHfit3$phiu)
### Equation to calculate the tail index
par(mfrow = c(1, 1))
tailindexH <- function(x) (1-Fhosp(x))/(x * fhosp(x))
curve(expr = tailindexH, from = 1, to = 350, ylim = c(0,2),
      ylab = "Limit", xlab = "x", mgp=c(2.5,1,0), lwd = 2,
      main = "Estimation of the tail index for hospitalization
services") 

################### Code 31: Value at Risk for hospitalization services
quantileH <- function(kappa) 1 - ((1 - kappa) / MeanNhGPO)
FquantileH <- function(kappa) qweibullgpd(p = quantileH(kappa),
                                          wshape = SpHfit3$wshape, wscale = SpHfit3$wscale,
                                          sigmau = SpHfit3$sigmau, u = SpHfit3$u, xi = SpHfit3$xi,
                                          phiu = SpHfit3$phiu)
correctionH <- MeanXhSpW * (MeanNhGPO + (VariNhGPO / MeanNhGPO) - 1)
VaRH <- function(kappa) FquantileH(kappa) + correctionH

curve(expr = VaRH, from = 0.90, to = 1, 
      ylab = "Value at Risk (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "Value at Risk for hospitalization services")
rect(xleft = 0.95, ybottom = VaRH(0.95) + 15, xright = 0.95,
     ytop = VaRH(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = VaRH(0.99) + 15, xright = 0.99,
     ytop = VaRH(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(VaRH(0.95), VaRH(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.94, y = VaRH(0.95) + 30, bty = "n",
       legend = round(VaRH(0.95), 3))
legend(x = 0.98, y = VaRH(0.99) + 30, bty = "n", 
       legend = round(VaRH(0.99), 3))

################### Code 32: Tail Value at Risk for hospitalization services
TVaRH <- function(kappa) (1/(1-kappa)) * as.numeric(integrate(
         f = VaRH, lower = kappa, upper = 1)$value)
TVaRH <- Vectorize(TVaRH)

curve(expr = TVaRH, from = 0.90, to = 1, 
      ylab = "Tail Value at Risk (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "Tail Value at Risk for hospitalization services")
rect(xleft = 0.95, ybottom = TVaRH(0.95) + 8, xright = 0.95,
     ytop = TVaRH(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = TVaRH(0.99) + 8, xright = 0.99,
     ytop = TVaRH(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(TVaRH(0.95), TVaRH(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.94, y = TVaRH(0.95) + 15, bty = "n",
       legend = round(TVaRH(0.95), 3))
legend(x = 0.98, y = TVaRH(0.99) + 15, bty = "n", 
       legend = round(TVaRH(0.99), 3))

################### Code 33: Expected Shortfall for hospitalization services
ESH <- function(kappa) (1 - kappa)*(TVaRH(kappa) - VaRH(kappa))

curve(expr = ESH, from = 0.90, to = 1, 
      ylab = "Expected Shortfall (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "Expected Shortfall for hospitalization services")
rect(xleft = 0.95, ybottom = ESH(0.95) + 0.2, xright = 0.95,
     ytop = ESH(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = ESH(0.99) + 0.2, xright = 0.99,
     ytop = ESH(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(ESH(0.95), ESH(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.942, y = ESH(0.95) + 0.35, bty = "n",
       legend = round(ESH(0.95), 3))
legend(x = 0.982, y = ESH(0.99) + 0.35, bty = "n", 
       legend = round(ESH(0.99), 3))

###################################### Risk measures estimation for general surgery services ######################################

################### Code 34: Estimation of the tail index for general surgery services
### Density function
fsurg <- function(x) dweibullgpd(x, wshape = SpSfit3$wshape, 
                                 wscale = SpSfit3$wscale, sigmau = SpSfit3$sigmau,
                                 u = SpSfit3$u, xi = SpSfit3$xi, phiu = SpSfit3$phiu)
### Cumulative function
Fsurg <- function(x) pweibullgpd(x, wshape = SpSfit3$wshape,
                                 wscale = SpSfit3$wscale, sigmau = SpSfit3$sigmau,
                                 u = SpSfit3$u, xi = SpSfit3$xi, phiu = SpSfit3$phiu)
### Equation to calculate the tail index
tailindexS <- function(x) (1-Fsurg(x))/(x * fsurg(x))
curve(expr = tailindexS, from = 1, to = 5e11, ylim = c(0.8,1),
      ylab = "Limit", xlab = "x", mgp=c(2.5,1,0), lwd = 2,
      main = "Estimation of the tail index for general surgery
services") 
abline(h = tailindexS(1e11), col = "red")
rect(xleft = 5e10, ybottom = tailindexS(1e11), xright = 5e10,
     ytop = tailindexS(1e11) + 0.02, lwd = 2, border = "red")
points(x = 5e10, y = tailindexS(1e11),
       pch = 19, col = "red", cex = 1.2)
legend(x = 1.2e10, y = tailindexS(1e11) + 0.04, bty = "n",
       legend = round(tailindexS(1e11), 3))

################### Code 35: Value at Risk for general surgery services
quantileS <- function(kappa) 1 - ((1 - kappa) / MeanNsDEL)
FquantileS <- function(kappa) qweibullgpd(p = quantileS(kappa),
                                          wshape = SpSfit3$wshape, wscale = SpSfit3$wscale,
                                          sigmau = SpSfit3$sigmau, u = SpSfit3$u, xi = SpSfit3$xi,
                                          phiu = SpSfit3$phiu)
correctionS <- MeanXhSpW * (MeanNsDEL + (VariNsDEL / MeanNsDEL) - 1)
VaRS <- function(kappa) FquantileS(kappa) + correctionS

curve(expr = VaRS, from = 0.90, to = 0.999, 
      ylab = "Value at Risk (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "Value at Risk for general surgery services")
rect(xleft = 0.95, ybottom = VaRS(0.95) + 1000, xright = 0.95,
     ytop = VaRS(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = VaRS(0.99) + 1000, xright = 0.99,
     ytop = VaRS(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(VaRS(0.95), VaRS(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.94, y = VaRS(0.95) + 1700, bty = "n",
       legend = round(VaRS(0.95), 3))
legend(x = 0.98, y = VaRS(0.99) + 1700, bty = "n", 
       legend = round(VaRS(0.99), 3))

################### Code 36: Tail Value at Risk for general surgery services
TVaRS <- function(kappa) (1/(1-kappa)) * as.numeric(integrate(
         f = VaRS, lower = kappa, upper = 1)$value)
TVaRS <- Vectorize(TVaRS)

curve(expr = TVaRS, from = 0.90, to = 0.999, 
      ylab = "Tail Value at Risk (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "Tail Value at Risk for general surgery services")
rect(xleft = 0.95, ybottom = TVaRS(0.95) + 10000, xright = 0.95,
     ytop = TVaRS(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = TVaRS(0.99) + 10000, xright = 0.99,
     ytop = TVaRS(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(TVaRS(0.95), TVaRS(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.94, y = TVaRS(0.95) + 15000, bty = "n",
       legend = round(TVaRS(0.95), 3))
legend(x = 0.98, y = TVaRS(0.99) + 15000, bty = "n", 
       legend = round(TVaRS(0.99), 3))

################### Code 37: Expected Shortfall for general surgery services
ESS <- function(kappa) (1 - kappa)*(TVaRS(kappa) - VaRS(kappa))

curve(expr = ESS, from = 0.90, to = 0.999, 
      ylab = "Expected Shortfall (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "Expected Shortfall for general surgery services")
rect(xleft = 0.95, ybottom = ESS(0.95) - 8, xright = 0.95,
     ytop = ESS(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = ESS(0.99) - 8, xright = 0.99,
     ytop = ESS(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(ESS(0.95), ESS(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.942, y = ESS(0.95) - 8, bty = "n",
       legend = round(ESS(0.95), 3))
legend(x = 0.982, y = ESS(0.99) - 8, bty = "n", 
       legend = round(ESS(0.99), 3))

###################################### Optimum retention point estimation for hospitalization service ######################################

################### Code 38: Optimum retention point estimation for hospitalization service
rho <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 3, 4, 5, 7, 10, 20, 50)
VaRTH <- function(rho, kappa) VaRH(kappa) + (1 + rho)*ESH(kappa)
ResultH <- function(rho){
  kapparho <- 1 - 1/(1 + rho)
  VaRHrho <- VaRH(kapparho)
  DeltaHrho <- (1+rho)*ESH(kapparho)
  VaRTHrho <- VaRTH(rho, kapparho)
  return(c(rho, kapparho, VaRHrho, DeltaHrho, VaRTHrho))
}

tableH <- round(t(sapply(X = rho, FUN = ResultH)), 6)

kable(cbind(data.frame(tableH)), escape = FALSE,
      col.names = c("$\\rho$", "$\\kappa_{\\rho^*}$",
                    "$M_{hosp}^*$", "$\\delta(M_{hosp}^*)$",
                    "$VaR_{T_{hosp}}(\\kappa_{\\rho^*})$"),
      caption = "Optimum retention point estimation for
      hospitalization services \\label{tab:retentionH}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"),
                position = "center") 

###################################### Optimum retention point estimation for general surgery services ######################################

################### Code 39: Optimum retention point estimation for hospitalization service
VaRTS <- function(rho, kappa) VaRS(kappa) + (1 + rho)*ESS(kappa)
ResultS <- function(rho){
  kapparho <- 1 - 1/(1 + rho)
  VaRSrho <- VaRS(kapparho)
  DeltaSrho <- (1+rho)*ESS(kapparho)
  VaRTSrho <- VaRTS(rho, kapparho)
  return(c(rho, kapparho, VaRSrho, DeltaSrho, VaRTSrho))
}

tableS <- round(t(sapply(X = rho, FUN = ResultS)), 6)

kable(cbind(data.frame(tableS)), escape = FALSE,
      col.names = c("$\\rho$", "$\\kappa_{\\rho^*}$",
                    "$M_{surg}^*$", "$\\delta(M_{surg}^*)$",
                    "$VaR_{T_{surg}}(\\kappa_{\\rho^*})$"),
      caption = "Optimum retention point estimation for
      general surgery services \\label{tab:retentionS}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"),
                position = "center") 

###################################### Premium approximation through Black-Scholes model for hospitalization services ######################################

################### Code 40: Adjustment for hospitalization services with Log-Normal distribution
### The adjustment with gamlss is made
LNHfit <- gamlssML(formula = hospita$cost, family = LOGNO)
### Mean
MeanXhLN <- moments(k = 1, dist = "LOGNO", domain = "realplus", 
                    param = c(mu = LNHfit$mu, sigma = LNHfit$sigma))

par(mai=rep(0.5, 4))
layout(matrix(c(1,1, 2,2, 0, 3,3, 0), ncol = 4, byrow = TRUE))
### Empirical vs Theorical cumulative distribution function
FnXh <- ecdf(hospita$cost)
sortXh <- sort(hospita$cost)
plot(FnXh, lwd = 3,
     xlab = "Sample quantiles of individual costs of hospitalization", 
     main = "(a) Adjustment cumulative individual costs 
     of hospitalization with Log-Normal",
     mgp=c(2, 1, 0))

fitXhL <- pLOGNO(q = sortXh, mu = LNHfit$mu, sigma = LNHfit$sigma)
lines(sortXh, fitXhL, lwd = 3, lty = 4, col = "red")
grid()
legend("bottomright", lty = 1, col = c("black", "red"),
       legend = c("Cumulative empirical distribution",
                  "Log-Normal"), lwd = 2)
### Empirical vs Theorical log-survival distribution function
survXh <- 1 - FnXh(sortXh)
plot(x = log(sortXh), y = log(survXh), lwd = 3,
     xlab = "log(Sample quantiles of individual cost of
     hospitalization)", ylab = "log(1 - Fn(x))", 
     main = "(b) Adjustment of log-survival distribution 
     of hospitalization with Log-Normal", type = "l")
survXhL <- 1 - fitXhL
lines(log(sortXh), log(survXhL), lwd = 3, col = "red", lty = 2)
grid()
legend("bottomleft", lty = 1, col = c("black", "red"),
       legend = c("Cumulative empirical distribution",
                  "Log-Normal"), lwd = 2)
### QQ-plot 
par(mgp=c(2, 1, 0))
qqPlot(x = hospita$cost, lwd = 1, distribution = "LOGNO",
       mu = LNHfit$mu, sigma = LNHfit$sigma, cex = 1, col.lines = "red",
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(c) Log-Normal", id = FALSE, ylim = c(0, 500))

###################################### Premium approximation through Black-Scholes model for hospitalization services ######################################

################### Code 41: Risk measures for hospitalization services with Log-Normal distribution
### Density function
fhospL <- function(x) dLOGNO(x, mu = LNHfit$mu, sigma = LNHfit$sigma)
### Cumulative function
FhospL <- function(x) pLOGNO(x, mu = LNHfit$mu, sigma = LNHfit$sigma)
### Equation to calculate the tail index
tailindexHL <- function(x) (1-FhospL(x))/(x * fhospL(x))
### Value at Risk
quantileHL <- function(kappa) 1 - ((1 - kappa) / MeanNhGPO)
FquantileHL <- function(kappa) qLOGNO(p = quantileHL(kappa), mu = LNHfit$mu,
                                      sigma = LNHfit$sigma)
correctionHL <- MeanXhLN * (MeanNhGPO + (VariNhGPO / MeanNhGPO) - 1)
VaRHL <- function(kappa) FquantileHL(kappa) + correctionHL
### Tail Value at Risk
TVaRHL <- function(kappa) (1/(1-kappa)) * as.numeric(integrate(f = VaRHL,
                                                               lower = kappa, upper = 1)$value)
TVaRHL <- Vectorize(TVaRHL)
### Expected Shortfall
ESHL <- function(kappa) (1 - kappa)*(TVaRHL(kappa) - VaRHL(kappa))

### Plots
par(mai=rep(0.5, 4), mfrow = c(2,2))
curve(expr = tailindexHL, from = 1, to = 4e5, ylim = c(0,2),
      ylab = "Limit", xlab = "x", mgp=c(2,1,0), lwd = 2,
      main = "(a) Estimation of the tail index for hospitalization
      services with Log-Normal distribution", cex.main = 0.9) 

curve(expr = VaRHL, from = 0.90, to = 1, 
      ylab = "Value at Risk (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "(b) Value at Risk for hospitalization services
      with Log-Normal distribution", 
      mgp=c(2,1,0), cex.main = 0.9)
rect(xleft = 0.95, ybottom = VaRHL(0.95) + 300, xright = 0.95,
     ytop = VaRHL(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = VaRHL(0.99) + 300, xright = 0.99,
     ytop = VaRHL(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(VaRHL(0.95), VaRHL(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.93, y = VaRHL(0.95) + 650, bty = "n",
       legend = round(VaRHL(0.95), 3))
legend(x = 0.97, y = VaRHL(0.99) + 650, bty = "n", 
       legend = round(VaRHL(0.99), 3))

curve(expr = TVaRHL, from = 0.90, to = 0.9999, 
      ylab = "Tail Value at Risk (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "(c) Tail Value at Risk for hospitalization services
with Log-Normal distribution", cex.main = 0.9)
rect(xleft = 0.95, ybottom = TVaRHL(0.95) + 800, xright = 0.95,
     ytop = TVaRHL(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = TVaRHL(0.99) + 800, xright = 0.99,
     ytop = TVaRHL(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(TVaRHL(0.95), TVaRHL(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.93, y = TVaRHL(0.95) + 1850, bty = "n",
       legend = round(TVaRHL(0.95), 3))
legend(x = 0.97, y = TVaRHL(0.99) + 1850, bty = "n", 
       legend = round(TVaRHL(0.99), 3))

curve(expr = ESHL, from = 0.90, to = 0.9999, 
      ylab = "Expected Shortfall (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "(d) Expected Shortfall for hospitalization services
with Log-Normal distribution", cex.main = 0.9)
rect(xleft = 0.95, ybottom = ESHL(0.95) + 5, xright = 0.95,
     ytop = ESHL(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = ESHL(0.99) + 5, xright = 0.99,
     ytop = ESHL(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(ESHL(0.95), ESHL(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.935, y = ESHL(0.95) + 10.3, bty = "n",
       legend = round(ESHL(0.95), 3))
legend(x = 0.975, y = ESHL(0.99) + 10.3, bty = "n", 
       legend = round(ESHL(0.99), 3))

################### Code 42: Optimum retention point estimation for hospitalization service with Log-Normal distribution
VaRTHL <- function(rho, kappa) VaRHL(kappa) + (1 + rho)*ESHL(kappa)
ResultHL <- function(rho){
  kapparho <- 1 - 1/(1 + rho)
  VaRHLrho <- VaRHL(kapparho)
  DeltaHLrho <- (1+rho)*ESHL(kapparho)
  VaRTHLrho <- VaRTHL(rho, kapparho)
  return(c(rho, kapparho, VaRHLrho, DeltaHLrho, VaRTHLrho))
}

tableHL <- round(t(sapply(X = rho, FUN = ResultHL)), 6)

kable(cbind(data.frame(tableHL)), escape = FALSE,
      col.names = c("$\\rho$", "$\\kappa_{\\rho^*}$",
                    "$M_{hosp}^*$", "$\\delta(M_{hosp}^*)$",
                    "$VaR_{T_{hosp}}(\\kappa_{\\rho^*})$"),
      caption = "Optimum retention point estimation for
      hospitalization services with Log-Normal distribution
      \\label{tab:retentionHL}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"),
                position = "center") 

###################################### Premium approximation through Black-Scholes model for general surgery services ######################################

################### Code 43: Adjustment for general surgery services with Log-Normal distribution
### The adjustment with evmix is made
LNSfit <- gamlssML(formula = surgery$cost, family = LOGNO)
### Mean
MeanXsLN <- moments(k = 1, dist = "LOGNO", domain = "realplus", 
                    param = c(mu = LNSfit$mu, sigma = LNSfit$sigma))

par(mai=rep(0.5, 4))
layout(matrix(c(1,1, 2,2, 0, 3,3, 0), ncol = 4, byrow = TRUE))
### Empirical vs Theorical cumulative distribution function
FnXh <- ecdf(surgery$cost)
sortXh <- sort(surgery$cost)
plot(FnXh, lwd = 3,
     xlab = "Sample quantiles of individual costs of general surgery", 
     main = "(a) Adjustment cumulative individual costs of
general surgery with Log-Normal", mgp=c(2,1,0))

fitXhL <- pLOGNO(q = sortXh, mu = LNSfit$mu, sigma = LNSfit$sigma)
lines(sortXh, fitXhL, lwd = 3, lty = 4, col = "red")
grid()
legend("bottomright", lty = 1, col = c("black", "red"),
       legend = c("Cumulative empirical distribution",
                  "Log-Normal"), lwd = 2)
### Empirical vs Theorical log-survival distribution function
survXh <- 1 - FnXh(sortXh)
plot(x = log(sortXh), y = log(survXh), lwd = 3,
     xlab = "log(Sample quantiles of individual cost of
     general surgery)", ylab = "log(1 - Fn(x))", 
     main = "(b) Adjustment of log-survival distribution
of general sugery with Log-Normal", type = "l")
survXhL <- 1 - fitXhL
lines(log(sortXh), log(survXhL), lwd = 3, col = "red", lty = 2)

grid()
legend("bottomleft", lty = 1, col = c("black", "red"),
       legend = c("Cumulative empirical distribution",
                  "Log-Normal"), lwd = 2)
### QQ-plot 
par(mgp=c(2, 1, 0))
qqPlot(x = surgery$cost, lwd = 1, distribution = "LOGNO",
       mu = LNSfit$mu, sigma = LNSfit$sigma, cex = 1, col.lines = "red" ,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(c) Log-Normal", id = FALSE, ylim = c(0, 500))

################### Code 44: Risk measures for general surgery services with Log-Normal distribution
### Density function
fsurgL <- function(x) dLOGNO(x, mu = LNSfit$mu, sigma = LNSfit$sigma)
### Cumulative function
FsurgL <- function(x) pLOGNO(x, mu = LNSfit$mu, sigma = LNSfit$sigma)
### Equation to calculate the tail index
tailindexSL <- function(x) (1-FsurgL(x))/(x * fsurgL(x))
### Value at Risk
quantileSL <- function(kappa) 1 - ((1 - kappa) / MeanNsDEL)
FquantileSL <- function(kappa) qLOGNO(p = quantileSL(kappa), mu = LNSfit$mu,
                                      sigma = LNSfit$sigma)
correctionSL <- MeanXsLN * (MeanNsDEL + (VariNsDEL / MeanNsDEL) - 1)
VaRSL <- function(kappa) FquantileSL(kappa) + correctionSL
### Tail Value at Risk
TVaRSL <- function(kappa) (1/(1-kappa)) * as.numeric(integrate(f = VaRSL,
                                                               lower = kappa, upper = 1)$value)
TVaRSL <- Vectorize(TVaRSL)
### Expected Shortfall
ESSL <- function(kappa) (1 - kappa)*(TVaRSL(kappa) - VaRSL(kappa))

### Plots
par(mai=rep(0.5, 4), mfrow = c(2,2))
curve(expr = tailindexSL, from = 1, to = 4e5, ylim = c(0,2),
      ylab = "Limit", xlab = "x", mgp=c(2,1,0), lwd = 2,
      main = "(a) Estimation of the tail index for general surgery
      services with Log-Normal distribution", cex.main = 0.9) 

curve(expr = VaRSL, from = 0.90, to = 1, 
      ylab = "Value at Risk (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "(b) Value at Risk for general surgery services
      with Log-Normal distribution", cex.main = 0.9)
rect(xleft = 0.95, ybottom = VaRSL(0.95) + 300, xright = 0.95,
     ytop = VaRSL(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = VaRSL(0.99) + 300, xright = 0.99,
     ytop = VaRSL(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(VaRSL(0.95), VaRSL(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.93, y = VaRSL(0.95) + 450, bty = "n",
       legend = round(VaRSL(0.95), 3))
legend(x = 0.97, y = VaRSL(0.99) + 450, bty = "n", 
       legend = round(VaRSL(0.99), 3))

curve(expr = TVaRSL, from = 0.90, to = 0.9999, 
      ylab = "Tail Value at Risk (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "(c) Tail Value at Risk for general surgery services
with Log-Normal distribution", cex.main = 0.9)
rect(xleft = 0.95, ybottom = TVaRSL(0.95) + 600, xright = 0.95,
     ytop = TVaRSL(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = TVaRSL(0.99) + 600, xright = 0.99,
     ytop = TVaRSL(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(TVaRSL(0.95), TVaRSL(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.93, y = TVaRSL(0.95) + 1200, bty = "n",
       legend = round(TVaRSL(0.95), 3))
legend(x = 0.97, y = TVaRSL(0.99) + 1200, bty = "n", 
       legend = round(TVaRSL(0.99), 3))

curve(expr = ESSL, from = 0.90, to = 0.9999, 
      ylab = "Expected Shortfall (in millions of pesos)",
      xlab = expression(kappa), mgp=c(2,1,0), lwd = 2,
      main = "(d) Expected Shortfall for general surgery services
with Log-Normal distribution", cex.main = 0.9)
rect(xleft = 0.95, ybottom = ESSL(0.95) + 2, xright = 0.95,
     ytop = ESSL(0.95), lwd = 2, border = "red")
rect(xleft = 0.99, ybottom = ESSL(0.99) + 2, xright = 0.99,
     ytop = ESSL(0.99), lwd = 2, border = "red")
points(x = c(0.95, 0.99), y = c(ESSL(0.95), ESSL(0.99)),
       pch = 19, col = "red", cex = 1.2)
legend(x = 0.933, y = ESSL(0.95) + 4.5, bty = "n",
       legend = round(ESSL(0.95), 3))
legend(x = 0.975, y = ESSL(0.99) + 4.5, bty = "n", 
       legend = round(ESSL(0.99), 3))

################### Code 45: Optimum retention point estimation for general surgery service with Log-Normal distribution
VaRTSL <- function(rho, kappa) VaRSL(kappa) + (1 + rho)*ESSL(kappa)
ResultSL <- function(rho){
  kapparho <- 1 - 1/(1 + rho)
  VaRSLrho <- VaRSL(kapparho)
  DeltaSLrho <- (1+rho)*ESSL(kapparho)
  VaRTSLrho <- VaRTSL(rho, kapparho)
  return(c(rho, kapparho, VaRSLrho, DeltaSLrho, VaRTSLrho))
}

tableSL <- round(t(sapply(X = rho, FUN = ResultSL)), 6)

kable(cbind(data.frame(tableSL)), escape = FALSE,
      col.names = c("$\\rho$", "$\\kappa_{\\rho^*}$",
                    "$M_{surg}^*$", "$\\delta(M_{surg}^*)$",
                    "$VaR_{T_{surg}}(\\kappa_{\\rho^*})$"),
      caption = "Optimum retention point estimation for
      general surgery services with Log-Normal distribution
      \\label{tab:retentionSL}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"),
                position = "center") 

###################################### GAMLSS distributions to adjust severity distributions ######################################

################### Code 46: Figure showing the shape of the tail for diferent types of GAMLSS distributions for $k_1; k_3; k_5$ = 1; 2, and $k_2; k_4; k_6$. Smaller values $i$ the $k$'s result heavier tails. Rigby et al. (2014)
### Types of distributions
TypeI <- function(x, k1, k2) - k2 * log(abs(x))^k1
TypeII <- function(x, k3, k4) - k4 * abs(x)^k3
TypeIII <- function(x, k5, k6) - k6 * exp(k5 * abs(x))
x <- seq(0, 10, 0.001)
x2 <- seq(1, 10, 0.001)
par(mfrow = c(2, 2), cex.main = 0.9, lwd = 2)

plot(x = x, y = TypeIII(x = x, k5 = 1, k6 = 1), type = "l",
     main = expression(paste(k[1], ", ", k[3], ", ", k[5], " = ", 1,
                             "  and  ", k[2], ", ", k[4], ", ", k[6], " = ", 1)),  lty = 4,
     ylab = "Types of heavy  tails", col = "blue", ylim = c(-10, 1))
legend("topright", legend = c("Type I", "Type II", "Type III"),
       col = c("black", "red", "blue"), lty = c(1,2,4))
lines(x = x, y = TypeII(x = x, k3 = 1, k4 = 1), type = "l",
      col = "red", lwd = 2, lty = 2)
lines(x = x, y = TypeI(x = x, k1 = 1, k2 = 1), type = "l", 
      col = "black", lwd = 2, lty = 1)

plot(x = x, y = TypeIII(x = x, k5 = 1, k6 = 2), type = "l",
     main = expression(paste(k[1], ", ", k[3], ", ", k[5], " = ", 1,
                             "  and  ", k[2], ", ", k[4], ", ", k[6], " = ", 2)), lty = 4,
     ylab = "Types of heavy  tails", col = "blue", ylim = c(-10, 1))
legend("topright", legend = c("Type I", "Type II", "Type III"),
       col = c("black", "red", "blue"), lty = c(1,2,4))
lines(x = x, y = TypeII(x = x, k3 = 1, k4 = 2), type = "l",
      col = "red", lwd = 2, lty = 2)
lines(x = x, y = TypeI(x = x, k1 = 1, k2 = 2), type = "l",
      col = "black", lwd = 2, lty = 1)

plot(x = x, y = TypeIII(x = x, k5 = 2, k6 = 1), type = "l",
     main = expression(paste(k[1], ", ", k[3], ", ", k[5], " = ", 2,
                             "  and  ", k[2], ", ", k[4], ", ", k[6], " = ", 1)), lty = 4,
     ylab = "Types of heavy  tails", col = "blue", ylim = c(-10, 1))
legend("topright", legend = c("Type I", "Type II", "Type III"),
       col = c("black", "red", "blue"), lty = c(1,2,4))
lines(x = x, y = TypeII(x = x, k3 = 2, k4 = 1), type = "l",
      col = "red", lwd = 2, lty = 2)
lines(x = x2, y = TypeI(x = x2, k1 = 2, k2 = 1), type = "l", 
      col = "black", lwd = 2, lty = 1)

plot(x = x, y = TypeIII(x = x, k5 = 2, k6 = 2), type = "l",
     main = expression(paste(k[1], ", ", k[3], ", ", k[5], " = ", 2,
                             "  and  ", k[2], ", ", k[4], ", ", k[6], " = ", 2)), lty = 4, 
     ylab = "Types of heavy  tails", col = "blue", ylim = c(-10, 1))
legend("topright", legend = c("Type I", "Type II", "Type III"),
       col = c("black", "red", "blue"), lty = c(1,2,4))
lines(x = x, y = TypeII(x = x, k3 = 2, k4 = 2), type = "l",
      col = "red", lwd = 2, lty = 2)
lines(x = x2, y = TypeI(x = x2, k1 = 2, k2 = 2), type = "l",
      col = "black", lwd = 2, lty = 1)

###################################### GAMLSS distribution for hospitalization ######################################

################### Code 47: Best fit with GAMLSS distributions for individual cost of hospitalization services
### The adjustment is made
GdHfit1 <- fitDist(y = hospita$cost, type = "realplus")
### Estimation of second and third distribution with best fit
GdHfit2 <- gamlssML(formula = hospita$cost, family = GB2)
GdHfit3 <- gamlssML(formula = hospita$cost, family = BCPE)
### The five distributions that present the best fit are
kable(rbind(GdHfit1$fits[1:5]),
      caption = "Better fit for individual cost of hospitalization
      services with GAMLSS distributions \\label{tab:fitfhosp}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))


################### Code 48: Statistical measurements of GAMLSS distributions for hospitalization services
### Estimation of mean, variance, skewness and excess kurtosis
#### Mean
MeanXhGdE <- mean(hospita$cost)
MeanXhGdGG <- moments(k = 1, dist = "GG", domain = "realplus", 
                      param = c(mu = GdHfit1$mu, sigma = GdHfit1$sigma,
                                nu = GdHfit1$nu))
MeanXhGdGB2 <- moments(k = 1, dist = "GB2", domain = "realplus", 
                       param = c(mu = GdHfit2$mu, sigma = GdHfit2$sigma,
                                 nu = GdHfit2$nu, tau = GdHfit2$tau))
MeanXhGdBCPE <- moments(k = 1, dist = "BCPE", domain = "realplus", 
                        param = c(mu = GdHfit3$mu, sigma = GdHfit3$sigma,
                                  nu = GdHfit3$nu, tau = GdHfit3$tau))
#### Variance
VariXhGdE <- var(hospita$cost)
VariXhGdGG <- moments(k = 2, dist = "GG", domain = "realplus", 
                      param = c(mu = GdHfit1$mu, sigma = GdHfit1$sigma,
                                nu = GdHfit1$nu), central = TRUE)
VariXhGdGB2 <- moments(k = 2, dist = "GB2", domain = "realplus", 
                       param = c(mu = GdHfit2$mu, sigma = GdHfit2$sigma,
                                 nu = GdHfit2$nu, tau = GdHfit2$tau), central = TRUE)
VariXhGdBCPE <- moments(k = 2, dist = "BCPE", domain = "realplus", 
                        param = c(mu = GdHfit3$mu, sigma = GdHfit3$sigma,
                                  nu = GdHfit3$nu, tau = GdHfit3$tau), central = TRUE)
### Skewness
SkewXhGdE <- skewness(hospita$cost, type = 1)
SkewXhGdGG <- skew(dist = "GG", domain = "realplus", 
                   param = c(mu = GdHfit1$mu, sigma = GdHfit1$sigma,
                             nu = GdHfit1$nu))
SkewXhGdGB2 <- skew(dist = "GB2", domain = "realplus", 
                    param = c(mu = GdHfit2$mu, sigma = GdHfit2$sigma,
                              nu = GdHfit2$nu, tau = GdHfit2$tau))
SkewXhGdBCPE <- skew(dist = "BCPE", domain = "realplus", 
                     param = c(mu = GdHfit3$mu, sigma = GdHfit3$sigma,
                               nu = GdHfit3$nu, tau = GdHfit3$tau))
### Excess Kurtosis
KurtXhGdE <- kurtosis(hospita$cost, type = 1)
KurtXhGdGG <- kurt(dist = "GG", domain = "realplus", 
                   param = c(mu = GdHfit1$mu, sigma = GdHfit1$sigma,
                             nu = GdHfit1$nu), excess = TRUE)
KurtXhGdGB2 <- kurt(dist = "GB2", domain = "realplus", 
                    param = c(mu = GdHfit2$mu, sigma = GdHfit2$sigma,
                              nu = GdHfit2$nu, tau = GdHfit2$tau), excess = TRUE)
KurtXhGdBCPE <- kurt(dist = "BCPE", domain = "realplus", 
                     param = c(mu = GdHfit3$mu, sigma = GdHfit3$sigma,
                               nu = GdHfit3$nu, tau = GdHfit3$tau), excess = TRUE)

kable(cbind(data.frame(Dist = c("Empirical", "GG", "GB2",
                                "BCPE")), "Mean" = c(MeanXhGdE, unname(MeanXhGdGG),
                                                     unname(MeanXhGdGB2), unname(MeanXhGdBCPE)),
            "Variance" = c(VariXhGdE, unname(VariXhGdGG),
                           unname(VariXhGdGB2), unname(VariXhGdBCPE)),
            "Skewness" = c(SkewXhGdE, unname(SkewXhGdGG),
                           unname(SkewXhGdGB2), unname(SkewXhGdBCPE)),
            "Excess Kurtosis" = c(KurtXhGdE, unname(KurtXhGdGG),
                                  unname(KurtXhGdGB2), unname(KurtXhGdBCPE))), 
      caption = "Statistical measurements of GAMLSS distributions for
      hospitalization services \\label{tab:StatisticsGdXh}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position")) 

################### Code 49: Adjustment of cumulative individual costs of hospitalization services with GAMLSS distributions
FnXh <- ecdf(hospita$cost)
sortXh <- sort(hospita$cost)
### Empirical vs Theorical cumulative distribution function
plot(FnXh, lwd = 3,
     xlab = "Sample quantiles of individual costs of hospitalization", 
     main = "Adjustment cumulative individual costs of hospitalization")
fitXhGG <- pGG(q = sortXh, mu = GdHfit1$mu, sigma = GdHfit1$sigma,
               nu = GdHfit1$nu)
lines(sortXh, fitXhGG, lwd = 3, lty = 1, col = "blue")
fitXhGB2 <- pGB2(q = sortXh, mu = GdHfit2$mu, sigma = GdHfit2$sigma,
                 nu = GdHfit2$nu, tau = GdHfit2$tau)
lines(sortXh, fitXhGB2, lwd = 3, lty = 2, col = "red")
fitXhBCPE <- pBCPE(q = sortXh, mu = GdHfit3$mu,
                   sigma = GdHfit3$sigma, nu = GdHfit3$nu, tau = GdHfit3$tau)
lines(sortXh, fitXhBCPE, lwd = 3, lty = 4, col = "green")
grid()
legend("bottomright", lty = 1, col = c("black", "blue", "red",
                                       "green"), legend = c("Cumulative empirical distribution",
                                                            "Generalized Gamma", "Generalized Beta type 2",
                                                            "Box-Cox Power Exponential"), lwd = 2)

################### Code 50: Adjustment of log-survival distribution of hospitalization with GAMLSS distributions
### Empirical vs Theorical log-survival distribution function
survXh <- 1 - FnXh(sortXh)
plot(x = log(sortXh), y = log(survXh), lwd = 3,
     xlab = "log(Sample quantiles of individual cost of
     hospitalization)", ylab = "log(1 - Fn(x))", 
     main = "Adjustment of log-survival distribution of
     hospitalization", type = "l")
survXhGG <- 1 - fitXhGG
lines(log(sortXh), log(survXhGG), lwd = 3, col = "blue")
survXhGB2 <- 1 - fitXhGB2
lines(log(sortXh), log(survXhGB2), lwd = 3, col = "red", lty = 2)
survXhBCPE <- 1 - fitXhBCPE
lines(log(sortXh), log(survXhBCPE), lwd = 3, col = "green", lty = 4)
grid()
legend("bottomleft", lty = 1, col = c("black", "blue", "red", "green"),
       legend = c("log-Survival Distribution", "Generalized Gamma",
                  "Generalized Beta type 2", "Box-Cox Power Exponential"),
       lwd = 2)

################### Code 51: Q-Q plot GAMLSS distribution for hospitalization
par(mai=rep(0.5, 4))
layout(matrix(c(1,1, 2,2, 0, 3,3, 0), ncol = 4, byrow = TRUE))
### QQ-plot 
qqPlot(x = hospita$cost, lwd = 1, distribution = "GG",
       mu = GdHfit1$mu, sigma = GdHfit1$sigma, nu = GdHfit1$nu,
       cex = 1, col.lines = "red" , xlab = "Theorical Quantiles",
       ylab = "Sample Quantiles", id = FALSE, 
       main = "(a) Generalized Gamma Q-Q Plot")
qqPlot(x = hospita$cost, lwd = 1, distribution = "GB2",
       mu = GdHfit2$mu, sigma = GdHfit2$sigma, nu = GdHfit2$nu,
       tau = GdHfit2$tau, cex = 1, col.lines = "red",
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(b) Generalized Beta type 2 Q-Q Plot", id = FALSE)
qqPlot(x = hospita$cost, lwd = 1, distribution = "BCPE",
       mu = GdHfit3$mu, sigma = GdHfit3$sigma, nu = GdHfit3$nu,
       tau = GdHfit3$tau, cex = 1, col.lines = "red" ,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(c) Box-Cox Power Exponential Q-Q Plot", id = FALSE)

################### Code 52: Goodness-of-fit tests for hospitalization services for GAMLSS
### Kolmogorov-Smirnov test
set.seed(1248) # A seed is established so that the results can be replicated
kolmGdXhGG <- ks.test(x = hospita$cost, distn =  "pGG", 
                      H = min(hospita$cost), fit = list(mu = GdHfit1$mu,
                                                        sigma = GdHfit1$sigma, nu = GdHfit1$nu))
kolmGdXhGB2 <- ks.test(x = hospita$cost, distn =  "pGB2", 
                       H = min(hospita$cost), fit = list(mu = GdHfit2$mu,
                                                         sigma = GdHfit2$sigma, nu = GdHfit2$nu,
                                                         tau = GdHfit2$tau))
kolmGdXhBCPE <- ks.test(x = hospita$cost, distn =  "pBCPE", 
                        H = min(hospita$cost), fit = list(mu = GdHfit3$mu,
                                                          sigma = GdHfit3$sigma, nu = GdHfit3$nu,
                                                          tau = GdHfit3$tau))
### Cramer-von Mises test
cramGdXhGG <- w2.test(x = hospita$cost, distn =  "pGG", 
                      H = min(hospita$cost), fit = list(mu = GdHfit1$mu,
                                                        sigma = GdHfit1$sigma, nu = GdHfit1$nu))
cramGdXhGB2 <- w2.test(x = hospita$cost, distn =  "pGB2", 
                       H = min(hospita$cost), fit = list(mu = GdHfit2$mu,
                                                         sigma = GdHfit2$sigma, nu = GdHfit2$nu,
                                                         tau = GdHfit2$tau))
cramGdXhBCPE <- w2.test(x = hospita$cost, distn =  "pBCPE", 
                        H = min(hospita$cost), fit = list(mu = GdHfit3$mu,
                                                          sigma = GdHfit3$sigma, nu = GdHfit3$nu,
                                                          tau = GdHfit3$tau))
### Kuiper test
kuipGdXhGG <- v.test(x = hospita$cost, distn =  "pGG", 
                     H = min(hospita$cost), fit = list(mu = GdHfit1$mu,
                                                       sigma = GdHfit1$sigma, nu = GdHfit1$nu))
kuipGdXhGB2 <- v.test(x = hospita$cost, distn =  "pGB2", 
                      H = min(hospita$cost), fit = list(mu = GdHfit2$mu,
                                                        sigma = GdHfit2$sigma, nu = GdHfit2$nu,
                                                        tau = GdHfit2$tau))
kuipGdXhBCPE <- v.test(x = hospita$cost, distn =  "pBCPE", 
                       H = min(hospita$cost), fit = list(mu = GdHfit3$mu,
                                                         sigma = GdHfit3$sigma, nu = GdHfit3$nu,
                                                         tau = GdHfit3$tau))
### Supremum class Upper Tail Anderson-Darling test
adupGdXhGG <- adup.test(x = hospita$cost, distn =  "pGG", 
                        H = min(hospita$cost), fit = list(mu = GdHfit1$mu,
                                                          sigma = GdHfit1$sigma, nu = GdHfit1$nu))
adupGdXhGB2 <- adup.test(x = hospita$cost, distn =  "pGB2", 
                         H = min(hospita$cost), fit = list(mu = GdHfit2$mu,
                                                           sigma = GdHfit2$sigma, nu = GdHfit2$nu,
                                                           tau = GdHfit2$tau))
adupGdXhBCPE <- adup.test(x = hospita$cost, distn =  "pBCPE", 
                          H = min(hospita$cost), fit = list(mu = GdHfit3$mu,
                                                            sigma = GdHfit3$sigma, nu = GdHfit3$nu,
                                                            tau = GdHfit3$tau))
### Quadratic Class Upper Tail Anderson-Darling test
ad2upGdXhGG <- ad2up.test(x = hospita$cost, distn =  "pGG", 
                          H = min(hospita$cost), fit = list(mu = GdHfit1$mu,
                                                            sigma = GdHfit1$sigma, nu = GdHfit1$nu))
ad2upGdXhGB2 <- ad2up.test(x = hospita$cost, distn =  "pGB2", 
                           H = min(hospita$cost), fit = list(mu = GdHfit2$mu,
                                                             sigma = GdHfit2$sigma, nu = GdHfit2$nu,
                                                             tau = GdHfit2$tau))
ad2upGdXhBCPE <- ad2up.test(x = hospita$cost, distn =  "pBCPE", 
                            H = min(hospita$cost), fit = list(mu = GdHfit3$mu,
                                                              sigma = GdHfit3$sigma, nu = GdHfit3$nu,
                                                              tau = GdHfit3$tau))
### Results table
kable(cbind(data.frame(Dist = c("GG", "GB2", "BCPE")),
            "ks.test" = c(kolmGdXhGG$p.value, kolmGdXhGB2$p.value,
                          kolmGdXhBCPE$p.value), "w2.test" = c(cramGdXhGG$p.value,
                                                               cramGdXhGB2$p.value, cramGdXhBCPE$p.value), "v.test" = c(
                                                                 kuipGdXhGG$p.value, kuipGdXhGB2$p.value,
                                                                 kuipGdXhBCPE$p.value), "adup.test" = c(adupGdXhGG$p.value,
                                                                                                        adupGdXhGB2$p.value, adupGdXhBCPE$p.value), 
            "ad2up.test" = c(ad2upGdXhGG$p.value, 
                             ad2upGdXhGB2$p.value, ad2upGdXhBCPE$p.value)),
      caption = "Goodness-of-fit tests hospitalization services
      for GAMLSS distributions \\label{tab:gftGdXh}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

###################################### GAMLSS distribution for general surgery ######################################

################### Code 53: Best fit with GAMLSS distributions for individual cost of general surgery services
### The adjustment is made
GdSfit1 <- fitDist(y = surgery$cost, type = "realplus")
### Estimation of second and third distribution with best fit
GdSfit2 <- gamlssML(formula = surgery$cost, family = BCPE)
GdSfit3 <- gamlssML(formula = surgery$cost, family = GG)
### The five distributions that present the best fit are
kable(rbind(GdSfit1$fits[1:5]),
      caption = "Better fit for individual cost of general surgery
      services with GAMLSS distributions \\label{tab:GdSfit}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

################### Code 54: Statistical measurements of GAMLSS distributions for general surgery services
### Estimation of mean, variance, skewness and excess kurtosis
#### Mean
MeanXsGdE <- mean(surgery$cost)
MeanXsGdBCPEo <- moments(k = 1, dist = "BCPEo", domain = "realplus", 
                         param = c(mu = GdSfit1$mu, sigma = GdSfit1$sigma,
                                   nu = GdSfit1$nu, tau = GdSfit1$tau))
MeanXsGdBCPE <- moments(k = 1, dist = "BCPE", domain = "realplus", 
                        param = c(mu = GdSfit2$mu, sigma = GdSfit2$sigma,
                                  nu = GdSfit2$nu, tau = GdSfit2$tau))
MeanXsGdGG <- moments(k = 1, dist = "GG", domain = "realplus", 
                      param = c(mu = GdSfit3$mu, sigma = GdSfit3$sigma,
                                nu = GdSfit3$nu))
#### Variance
VariXsGdE <- var(surgery$cost)
VariXsGdBCPEo <- moments(k = 2, dist = "BCPEo", domain = "realplus", 
                         param = c(mu = GdSfit1$mu, sigma = GdSfit1$sigma,
                                   nu = GdSfit1$nu, tau = GdSfit1$tau), central = TRUE)
VariXsGdBCPE <- moments(k = 2, dist = "BCPE", domain = "realplus", 
                        param = c(mu = GdSfit2$mu, sigma = GdSfit2$sigma,
                                  nu = GdSfit2$nu, tau = GdSfit2$tau), central = TRUE)
VariXsGdGG <- moments(k = 2, dist = "GG", domain = "realplus", 
                      param = c(mu = GdSfit3$mu, sigma = GdSfit3$sigma,
                                nu = GdSfit3$nu), central = TRUE)
### Skewness
SkewXsGdE <- skewness(surgery$cost, type = 1)
SkewXsGdBCPEo <- skew(dist = "BCPEo", domain = "realplus", 
                      param = c(mu = GdSfit1$mu, sigma = GdSfit1$sigma,
                                nu = GdSfit1$nu, tau = GdSfit1$tau))
SkewXsGdBCPE <- skew(dist = "BCPE", domain = "realplus", 
                     param = c(mu = GdSfit2$mu, sigma = GdSfit2$sigma,
                               nu = GdSfit2$nu, tau = GdSfit2$tau))
SkewXsGdGG <- skew(dist = "GG", domain = "realplus", 
                   param = c(mu = GdSfit3$mu, sigma = GdSfit3$sigma,
                             nu = GdSfit3$nu))
### Excess Kurtosis
KurtXsGdE <- kurtosis(surgery$cost, type = 1)
KurtXsGdBCPEo <- kurt(dist = "BCPEo", domain = "realplus", 
                      param = c(mu = GdSfit1$mu, sigma = GdSfit1$sigma,
                                nu = GdSfit1$nu, tau = GdSfit1$tau), excess = TRUE)
KurtXsGdBCPE <- kurt(dist = "BCPE", domain = "realplus", 
                     param = c(mu = GdSfit2$mu, sigma = GdSfit2$sigma,
                               nu = GdSfit2$nu, tau = GdSfit2$tau), excess = TRUE)
KurtXsGdGG <- kurt(dist = "GG", domain = "realplus", 
                   param = c(mu = GdSfit3$mu, sigma = GdSfit3$sigma,
                             nu = GdSfit3$nu), excess = TRUE)

kable(cbind(data.frame(Dist = c("Empirical", "GG", "BCPEo",
                                "BCPE")), "Mean" = c(MeanXsGdE, unname(MeanXsGdBCPEo),
                                                     unname(MeanXsGdBCPE), unname(MeanXsGdGG)),
            "Variance" = c(VariXsGdE, unname(VariXsGdBCPEo),
                           unname(VariXsGdBCPE), unname(VariXsGdGG)),
            "Skewness" = c(SkewXsGdE, unname(SkewXsGdBCPEo),
                           unname(SkewXsGdBCPE), unname(SkewXsGdGG)),
            "Excess Kurtosis" = c(KurtXsGdE, unname(KurtXsGdBCPEo),
                                  unname(KurtXsGdBCPE), unname(KurtXsGdGG))), 
      caption = "Statistical measurements of GAMLSS distributions for
      general surgery services \\label{tab:StatisticsGdXs}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position")) 

################### Code 55: Adjustment of cumulative individual costs of general surgery services with GAMLSS distributions
FnXs <- ecdf(surgery$cost)
sortXs <- sort(surgery$cost)
### Empirical vs Theorical cumulative distribution function
par(mfrow = c(1, 1))
plot(FnXs, lwd = 3,
     xlab = "Sample quantiles of individual costs of general surgery", 
     main = "Adjustment cumulative individual costs of general surgery")
fitXsBCPEo <- pBCPEo(q = sortXs, mu = GdSfit1$mu,
                     sigma = GdSfit1$sigma, nu = GdSfit1$nu,
                     tau = GdSfit1$tau)
lines(sortXs, fitXsBCPEo, lwd = 3, lty = 1, col = "blue")
fitXsBCPE <- pBCPE(q = sortXs, mu = GdSfit2$mu, sigma = GdSfit2$sigma,
                   nu = GdSfit2$nu, tau = GdSfit2$tau)
lines(sortXs, fitXsBCPE, lwd = 3, lty = 2, col = "red")
fitXsGG <- pGG(q = sortXs, mu = GdSfit3$mu, sigma = GdSfit3$sigma,
               nu = GdSfit3$nu)
lines(sortXs, fitXsGG, lwd = 3, lty = 4, col = "green")
grid()
legend("bottomright", lty = 1, col = c("black", "blue", "red",
                                       "green"), legend = c("Cumulative empirical distribution",
                                                            "Box-Cox Power Exponential-Original",
                                                            "Box-Cox Power Exponential", "Generalized Gamma"), lwd = 2)


################### Code 56: Adjustment of log-survival distribution of general surgery with GAMLSS distributions
### Empirical vs Theorical log-survival distribution function
survXs <- 1 - FnXs(sortXs)
plot(x = log(sortXs), y = log(survXs), lwd = 3,
     xlab = "log(Sample quantiles of individual cost of
     general surgery)", ylab = "log(1 - Fn(x))", 
     main = "Adjustment of log-survival distribution of
     general surgery", type = "l")
survXsBCPEo <- 1 - fitXsBCPEo
lines(log(sortXs), log(survXsBCPEo), lwd = 3, col = "blue")
survXsBCPE <- 1 - fitXsBCPE
lines(log(sortXs), log(survXsBCPE), lwd = 3, col = "red", lty = 2)
survXsGG <- 1 - fitXsGG
lines(log(sortXs), log(survXsGG), lwd = 3, col = "green", lty = 4)
grid()
legend("bottomleft", lty = 1, col = c("black", "blue", "red", "green"),
       legend = c("log-Survival Distribution",
                  "Box-Cox Power Exponential-Original",
                  "Box-Cox Power Exponential", "Generalized Gamma"),
       lwd = 2)

################### Code 57: Q-Q plot GAMLSS distribution for general surgery
par(mai=rep(0.5, 4))
layout(matrix(c(1,1, 2,2, 0, 3,3, 0), ncol = 4, byrow = TRUE))
### QQ-plot 
qqPlot(x = surgery$cost, lwd = 1, distribution = "BCPEo",
       mu = GdSfit1$mu, sigma = GdSfit1$sigma, nu = GdSfit1$nu,
       tau = GdSfit1$tau, cex = 1, col.lines = "red", id = FALSE,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(a) Box-Cox Power Exponential-Orig. Q-Q Plot")
qqPlot(x = surgery$cost, lwd = 1, distribution = "BCPE",
       mu = GdSfit2$mu, sigma = GdSfit2$sigma, nu = GdSfit2$nu,
       tau = GdSfit2$tau, cex = 1, col.lines = "red", id = FALSE,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(b) Box-Cox Power Exponential Q-Q Plot")
qqPlot(x = surgery$cost, lwd = 1, distribution = "GG",
       mu = GdSfit3$mu, sigma = GdSfit3$sigma, nu = GdSfit3$nu,
       cex = 1, col.lines = "red", id = FALSE,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles", 
       main = "(c) Generalized Gamma Q-Q Plot")

################### Code 58: Goodness-of-fit tests for general surgery services for GAMLSS distributions
set.seed(1248) # A seed is established so that the results can be replicated
### Kolmogorov-Smirnov test
kolmGdXsBCPEo <- ks.test(x = surgery$cost, distn =  "pBCPEo", 
                         H = min(surgery$cost), fit = list(mu = GdSfit1$mu,
                                                           sigma = GdSfit1$sigma,nu = GdSfit1$nu,
                                                           tau = GdSfit1$tau))
kolmGdXsBCPE <- ks.test(x = surgery$cost, distn =  "pBCPE", 
                        H = min(surgery$cost), fit = list(mu = GdSfit2$mu,
                                                          sigma = GdSfit2$sigma, nu = GdSfit2$nu,
                                                          tau = GdSfit2$tau))
kolmGdXsGG <- ks.test(x = surgery$cost, distn =  "pGG", 
                      H = min(surgery$cost), fit = list(mu = GdSfit3$mu,
                                                        sigma = GdSfit3$sigma, nu = GdSfit3$nu))
### Cramer-von Mises test
cramGdXsBCPEo <- w2.test(x = surgery$cost, distn =  "pBCPEo", 
                         H = min(surgery$cost), fit = list(mu = GdSfit1$mu,
                                                           sigma = GdSfit1$sigma,nu = GdSfit1$nu,
                                                           tau = GdSfit1$tau))
cramGdXsBCPE <- w2.test(x = surgery$cost, distn =  "pBCPE", 
                        H = min(surgery$cost), fit = list(mu = GdSfit2$mu,
                                                          sigma = GdSfit2$sigma, nu = GdSfit2$nu,
                                                          tau = GdSfit2$tau))
cramGdXsGG <- w2.test(x = surgery$cost, distn =  "pGG", 
                      H = min(surgery$cost), fit = list(mu = GdSfit3$mu,
                                                        sigma = GdSfit3$sigma, nu = GdSfit3$nu))
### Kuiper test
kuipGdXsBCPEo <- v.test(x = surgery$cost, distn =  "pBCPEo", 
                        H = min(surgery$cost), fit = list(mu = GdSfit1$mu,
                                                          sigma = GdSfit1$sigma,nu = GdSfit1$nu,
                                                          tau = GdSfit1$tau))
kuipGdXsBCPE <- v.test(x = surgery$cost, distn =  "pBCPE", 
                       H = min(surgery$cost), fit = list(mu = GdSfit2$mu,
                                                         sigma = GdSfit2$sigma, nu = GdSfit2$nu,
                                                         tau = GdSfit2$tau))
kuipGdXsGG <- v.test(x = surgery$cost, distn =  "pGG", 
                     H = min(surgery$cost), fit = list(mu = GdSfit3$mu,
                                                       sigma = GdSfit3$sigma, nu = GdSfit3$nu))
### Supremum class Upper Tail Anderson-Darling test
adupGdXsBCPEo <- adup.test(x = surgery$cost, distn =  "pBCPEo", 
                           H = min(surgery$cost), fit = list(mu = GdSfit1$mu,
                                                             sigma = GdSfit1$sigma,nu = GdSfit1$nu,
                                                             tau = GdSfit1$tau))
adupGdXsBCPE <- adup.test(x = surgery$cost, distn =  "pBCPE", 
                          H = min(surgery$cost), fit = list(mu = GdSfit2$mu,
                                                            sigma = GdSfit2$sigma, nu = GdSfit2$nu,
                                                            tau = GdSfit2$tau))
adupGdXsGG <- adup.test(x = surgery$cost, distn =  "pGG", 
                        H = min(surgery$cost), fit = list(mu = GdSfit3$mu,
                                                          sigma = GdSfit3$sigma, nu = GdSfit3$nu))
### Quadratic Class Upper Tail Anderson-Darling test
ad2upGdXsBCPEo <- ad2up.test(x = surgery$cost, distn =  "pBCPEo", 
                             H = min(surgery$cost), fit = list(mu = GdSfit1$mu,
                                                               sigma = GdSfit1$sigma,nu = GdSfit1$nu,
                                                               tau = GdSfit1$tau))
ad2upGdXsBCPE <- ad2up.test(x = surgery$cost, distn =  "pBCPE", 
                            H = min(surgery$cost), fit = list(mu = GdSfit2$mu,
                                                              sigma = GdSfit2$sigma, nu = GdSfit2$nu,
                                                              tau = GdSfit2$tau))
ad2upGdXsGG <- ad2up.test(x = surgery$cost, distn =  "pGG", 
                          H = min(surgery$cost), fit = list(mu = GdSfit3$mu,
                                                            sigma = GdSfit3$sigma, nu = GdSfit3$nu))
### Results table
kable(cbind(data.frame(Dist = c("BCPEo", "BCPE", "GG")),
            "ks.test" = c(kolmGdXsBCPEo$p.value, kolmGdXsBCPE$p.value,
                          kolmGdXsGG$p.value), "w2.test" = c(cramGdXsBCPEo$p.value,
                                                             cramGdXsBCPE$p.value, cramGdXsGG$p.value), "v.test" = c(
                                                               kuipGdXsBCPEo$p.value, kuipGdXsBCPE$p.value,
                                                               kuipGdXsGG$p.value), "adup.test" = c(adupGdXsBCPEo$p.value,
                                                                                                    adupGdXsBCPE$p.value, adupGdXsGG$p.value), "ad2up.test" =
              c(ad2upGdXsBCPEo$p.value, ad2upGdXsBCPE$p.value,
                ad2upGdXsGG$p.value)),
      caption = "Goodness-of-fit tests general surgery services for
      GAMLSS distributions \\label{tab:gftGdXs}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

###################################### Comparison of adjustment between the spliced distribution and GAMLSS distribution for hospitalization ######################################

################### Code 59: Comparison of statistical measurements for hospitalization services
kable(cbind(data.frame(Dist = c("Empirical", "weibullgpd", "GB2")),
            "Mean" = c(MeanXhGdE, unname(MeanXhSpW),
                       unname(MeanXhGdGB2)), "Variance" = c(VariXhGdE,
                                                            unname(VariXhSpW), unname(VariXhGdGB2)),
            "Skewness" = c(SkewXhGdE, unname(SkewXhSpW),
                           unname(SkewXhGdGB2)), "Excess Kurtosis" = c(KurtXhGdE,
                                                                       unname(KurtXhSpW), unname(KurtXhGdGB2))), 
      caption = "Comparison of statistical measurements for
      hospitalization services \\label{tab:StatisticsCpXs}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position")) 

################### Code 60: Comparison of adjustment of cumulative individual costs of hospitalization services
### Empirical vs Theorical cumulative distribution function
par(mfrow = c(1,1))
plot(FnXh, lwd = 3,
     xlab = "Sample quantiles of individual costs of hospitalization", 
     main = "Adjustment cumulative individual costs of hospitalization")
lines(sortXh, fitXhW, lwd = 3, lty = 1, col = "blue")
lines(sortXh, fitXhGB2, lwd = 3, lty = 2, col = "red")
grid()
legend("bottomright", lty = 1, col = c("black", "blue", "red"),
       lwd = 2, legend = c("Cumulative empirical distribution",
                           "Weibull-Generalized Pareto", "Generalized Beta type 2"))

################### Code 61: Comparison of adjustment of log-survival costs of hospitalization services
### Empirical vs Theorical log-survival distribution function
plot(x = log(sortXh), y = log(survXh), lwd = 3,
     xlab = "log(Sample quantiles of individual cost of
     hospitalization)", ylab = "log(1 - Fn(x))", 
     main = "Adjustment of log-survival distribution of
     hospitalization", type = "l")
lines(log(sortXh), log(survXhW), lwd = 3, col = "blue", lty = 1)
lines(log(sortXh), log(survXhGB2), lwd = 3, col = "red", lty = 2)
grid()
legend("bottomleft", lty = 1, col = c("black", "blue", "red"),
       lwd = 2, legend = c("log-Survival Distribution",
                           "Weibull-Generalized Pareto", "Generalized Beta type 2"))

################### Code 62: Comparison Q-Q plot for hospitalization services
par(mfrow = c(2,1))
### QQ-plot 
qqPlot(x = hospita$cost, lwd = 1, distribution = "weibullgpd",
       phiu = SpHfit3$phi, wshape = SpHfit3$wshape,
       wscale = SpHfit3$wscale, u = SpHfit3$u, xi = SpHfit3$xi,
       sigmau = SpHfit3$sigmau, cex = 1, col.lines = "red" ,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(a) Weibull-Generalized Pareto Q-Q Plot", id = FALSE)
qqPlot(x = hospita$cost, lwd = 1, distribution = "GB2",
       mu = GdHfit2$mu, sigma = GdHfit2$sigma, nu = GdHfit2$nu,
       tau = GdHfit2$tau, cex = 1, col.lines = "red",
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(b) Generalized Beta type 2 Q-Q Plot", id = FALSE)

################### Code 63: Comparision goodness-of-fit tests for hospitalization services
### Results table
kable(cbind(data.frame(Dist = c("weibullgpd", "GB2")),
            "ks.test" = c(kolmSpXhW$p.value, kolmGdXhGB2$p.value),
            "w2.test" = c(cramSpXhW$p.value, cramGdXhGB2$p.value),
            "v.test" = c(kuipSpXhW$p.value, kuipGdXhGB2$p.value), 
            "adup.test" = c(adupSpXhW$p.value, adupGdXhGB2$p.value), 
            "ad2up.test" = c(ad2upSpXhW$p.value, ad2upGdXhGB2$p.value)),
      caption = "Comparision goodness-of-fit tests hospitalization
      services
      \\label{tab:gftCpXh}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

###################################### Comparison of adjustment between the spliced distribution and GAMLSS distribution for general surgery ######################################

################### Code 64: Comparison of statistical measurements for general surgery services
kable(cbind(data.frame(Dist = c("Empirical", "weibullgpd",
                                "BCPEo")), "Mean" = c(MeanXsGdE, unname(MeanXsSpW),
                                                      unname(MeanXsGdBCPEo)), "Variance" = c(round(VariXsGdE,
                                                                                                   5), "does not exist", round(unname(VariXsGdBCPEo), 5)),
            "Skewness" = c(round(SkewXsGdE, 6), "does not exist",
                           round(unname(SkewXsGdBCPEo), 6)), "Excess Kurtosis" =
              c(round(KurtXsGdE, 5), "does not exist", round(unname(
                KurtXsGdBCPEo, 5)))), 
      caption = "Comparison of statistical measurements for
      general surgery services \\label{tab:StatisticsCpXs}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position")) 

################### Code 65: Comparison of adjustment of cumulative individual costs of general surgery services
### Empirical vs Theorical cumulative distribution function
par(mfrow = c(1,1))
plot(FnXs, lwd = 3,
     xlab = "Sample quantiles of individual costs of general surgery", 
     main = "Adjustment cumulative individual costs of general surgery")
lines(sortXs, fitXsW, lwd = 3, lty = 1, col = "blue")
lines(sortXs, fitXsBCPEo, lwd = 3, lty = 2, col = "red")
grid()
legend("bottomright", lty = 1, col = c("black", "blue", "red"),
       lwd = 2, legend = c("Cumulative empirical distribution",
                           "Weibull-Generalized Pareto",
                           "Box-Cox Power Exponential-Original"))

################### Code 66: Comparison of adjustment of log-survival costs of general surgery services
### Empirical vs Theorical log-survival distribution function
plot(x = log(sortXs), y = log(survXs), lwd = 3,
     xlab = "log(Sample quantiles of individual cost of
     general surgery)", ylab = "log(1 - Fn(x))", 
     main = "Adjustment of log-survival distribution of
     general surgery", type = "l")
lines(log(sortXs), log(survXsW), lwd = 3, col = "blue", lty = 1)
lines(log(sortXs), log(survXsBCPEo), lwd = 3, col = "red", lty = 2)
grid()
legend("bottomleft", lty = 1, col = c("black", "blue", "red"),
       lwd = 2, legend = c("log-Survival Distribution",
                           "Weibull-Generalized Pareto",
                           "Box-Cox Power Exponential-Original"))

################### Code 67: Comparison Q-Q plot for general surgery services
par(mfrow = c(2,1))
### QQ-plot 
qqPlot(x = surgery$cost, lwd = 1, distribution = "weibullgpd",
       phiu = SpSfit3$phi, wshape = SpSfit3$wshape,
       wscale = SpSfit3$wscale, u = SpSfit3$u, xi = SpSfit3$xi,
       sigmau = SpSfit3$sigmau, cex = 1, col.lines = "red" ,
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(a) Weibull-Generalized Pareto Q-Q Plot", id = FALSE)
qqPlot(x = surgery$cost, lwd = 1, distribution = "BCPEo",
       mu = GdSfit2$mu, sigma = GdSfit2$sigma, nu = GdSfit2$nu,
       tau = GdSfit2$tau, cex = 1, col.lines = "red",
       xlab = "Theorical Quantiles", ylab = "Sample Quantiles",
       main = "(b) Box-Cox Power Exponential-Orig. Q-Q Plot",
       id = FALSE)

################### Code 68: Comparision goodness-of-fit tests for general surgery services
### Results table
kable(cbind(data.frame(Dist = c("weibullgpd", "BCPEo")),
            "ks.test" = c(kolmSpXsW$p.value, kolmGdXsBCPEo$p.value),
            "w2.test" = c(cramSpXsW$p.value, cramGdXsBCPEo$p.value),
            "v.test" = c(kuipSpXsW$p.value, kuipGdXsBCPEo$p.value),  
            "adup.test" = c(adupSpXsW$p.value, adupGdXsBCPEo$p.value),   
            "ad2up.test" = c(ad2upSpXsW$p.value,
                             ad2upGdXsBCPEo$p.value)),
      caption = "Comparision goodness-of-fit tests general surgery
      services 
      \\label{tab:gftCpXs}",
      "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

