###################################################
## Experimental Design and Statistical Analsysis ##
###################################################

####################################
## PART 1 The experimental design ##
####################################
# Computation of the Z value for the 'Polio Experiment'
n1 <- 200745      # sample size of those who received Salk vaccine as treatment
n2 <- 201229      # sample size of those who received placebo as treatment
n  <- n1 + n2     # Total sample size
num_pol1 <- 56    # Number of children who contracted polio and received Salk vaccine
num_pol2 <- 142   # Number of children who contracted polio and received Placebo
num_pol  <- num_pol1+num_pol2 # Total number of children who contracted polio
P1 <- num_pol1/n1 # Proportion of children who were diagnosed with polio and received Salk vaccine 1-((n1-num_pol1)/n1)
P2 <- num_pol2/n2 # Proportion of children who were diagnosed with polio and received Placebo      1-((n2-num_pol2)/n2)
P  <- num_pol/n   
z  <- (P1 - P2)/sqrt((P*(1-P))*((1/n1)+(1/n2)))

# Power Requirement Data Frame
library(EMSaov) # Required library
power_cut <- structure(list(cut = structure(c(rep(2L,16),c(rep(1L,16))),
                                            .Label = c("Continuous","Interrupted"),
                                            class = "factor"),
                            tool   = structure(c(rep(c(rep(2L,8),rep(1L,8)),2)),
                                               .Label = c("Tool1","Tool2"),
                                               class = "factor"),
                            bevel = structure(c(rep(c(rep(2L,4),rep(1L,4)),4)),
                                              .Label = c("15","30"),
                                              class = "factor"),
                            mil_def = as.numeric(c(27.5, 28, 27, 26, 24.5, 25, 28,
                                                   26, 27, 29, 27.5, 27.5, 28, 25,
                                                   26.5, 26.5, 29.5, 32, 29, 28, 28,
                                                   28.5, 28, 25, 28.5, 28.5, 30,
                                                   32.5, 29, 26.5, 30.5, 27))),
                       .Names = c("cut","tool","bevel","mil_del"),
                       class = "data.frame",
                       row.names = c(NA, -32L))

# ANOVA for the 'Power Requirement Experiment' #change!
Anova_T <- EMSanova(mil_del~tool*bevel*cut,data=power_cut,
                    type=c("F","F","F"),nested=c(NA,NA,NA))

# Figure for the interaction of two factors of the 'Power Requirement Experiment'
library(ggplot2,plyr) # Required libraries
desc_power <- ddply(power_cut,.(cut,bevel),summarise,pw_mean = mean(mil_del),
                    pw_ci = 1.96*sd(mil_del)/sqrt(length(mil_del)))
pd <- position_dodge(width = 0.2)
(ansPlot <-  ggplot(data = desc_power,aes(x=bevel,y=pw_mean,group = cut))+
    geom_line(aes(linetype = cut), position = pd) +
    geom_errorbar(aes(ymin= pw_mean - pw_ci, ymax= pw_mean + pw_ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Type of Cut")) +
    labs(title = c(),x = c("Bevel (degrees)"),
         y = c("Sample mean")) +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)
##################################
## PART 2 Statistical Inference ##
##################################
# Figure on the Beta distribution (Type II Error)
par(mfrow=c(3,1)) ; 
n = 300
x1 <- rnorm(n,mean = 19.5, sd=0.4)
h1 <- hist(x1,breaks = 1, xlim = c(17,21),ylim=c(0,700),
           border = "white",main="",xlab = expression(mu == 19.5),ylab = "") 
x1fit <- seq(min(x1),max(x1),length=length(x1)) ; abline(v=19.5,lwd=3) ; abline(v=18.842)
y1fit <- dnorm(x1fit,mean=mean(x1),sd=sd(x1)) ; y1fit <- y1fit*diff(h1$mids[1:2])*n
lines(x1fit,y1fit,col="blue", lwd=2)
text(18,150,expression(alpha == 0.05)) ; text(21,100,expression(bar(y)))

x2 <- rnorm(n,mean = 19, sd=0.4)
h2 <- hist(x2,breaks = 1, xlim = c(17,21),ylim=c(0,700),
           border = "white",main="",xlab = expression(mu == 19.0),ylab = "") 
x2fit <- seq(min(x2),max(x2),length=length(x2)) ; abline(v=19.0,lwd=3) ; abline(v=18.842)
y2fit <- dnorm(x2fit,mean=mean(x2),sd=sd(x2)) ; y2fit <- y2fit*diff(h2$mids[1:2])*n
lines(x2fit,y2fit,col="blue", lwd=2)
text(20,600,expression(beta)) ; text(21,100,expression(bar(y)))

x3 <- rnorm(n,mean = 18.5, sd=0.4)
h3 <- hist(x3,breaks = 1, xlim = c(17,21),ylim=c(0,700),
           border = "white",main="",xlab = expression(mu == 18.5),ylab = "") 
x3fit <- seq(min(x3),max(x3),length=length(x3)) ; abline(v=18.5,lwd=3) ; abline(v=18.842)
y3fit <- dnorm(x3fit,mean=mean(x3),sd=sd(x3)) ; y3fit <- y3fit*diff(h3$mids[1:2])*n
lines(x3fit,y3fit,col="blue", lwd=2)
text(19.5,400,expression(beta)) ; text(21,100,expression(bar(y)))


############################
# # R CODE for Exercises # #
############################
rm(list = ls())
options(scipen = 99) # turn off scientific notation 
x <- c('ggplot2','plyr','EMSaov') # Listado de paquetes minimos para correr el Algoritmo
for (i in 1:length(x)) {
  library(x[i],character.only = TRUE)
}
cat("\014")
#####################
# # # Chapter 1 # # #
#####################
#' Example 1.5 # #

n1 <- 200745      # sample size of those who received Salk vaccine as treatment
n2 <- 201229      # sample size of those who received placebo as treatment
n <- n1 + n2      # Total sample size
num_pol1 <- 56    # Number of children who contracted polio and received Salk vaccine
num_pol2 <- 142   # Number of children who contracted polio and received Placebo
num_pol <- num_pol1+num_pol2 # Total number of children who contracted polio
P1 <- num_pol1/n1 # Proportion of children who were diagnosed with polio and received Salk vaccine 1-((n1-num_pol1)/n1)
P2 <- num_pol2/n2 # Proportion of children who were diagnosed with polio and received Placebo      1-((n2-num_pol2)/n2)
P  <- num_pol/n   
z  <- (P1 - P2)/sqrt((P*(1-P))*((1/n1)+(1/n2)))
pnorm(z)

#' When the null hypothesis is true, Z has an approximate normal distribution with mean 0 and variance 1
#' The probability of observing a value of Z as extreme as -6.09 when the null hypothesis is true is P(Z ≤ -6.09)
#' is equivalent to 4.126 x 10^-10. Thus the likelihood of observing a difference of this magnitude when the Salk
#' vaccine has no effect is less than one in a billion. This gives a very strong support to reject H0 and conclude
#' that the vaccine is effective.
#' 
# # Example 1.6 # #
#' The results in millimiter deflection are given below
library(plyr,EMSaov)
power_cut <- structure(list(cut = structure(c(rep(2L,16),c(rep(1L,16))),
                                            .Label = c("Continuous","Interrupted"),
                                            class = "factor"),
                            tool   = structure(c(rep(c(rep(2L,8),rep(1L,8)),2)),
                                               .Label = c("Tool1","Tool2"),
                                               class = "factor"),
                            bevel = structure(c(rep(c(rep(2L,4),rep(1L,4)),4)),
                                              .Label = c("15","30"),
                                              class = "factor"),
                            mil_def = as.numeric(c(27.5, 28, 27, 26, 24.5, 25, 28,
                                                   26, 27, 29, 27.5, 27.5, 28, 25,
                                                   26.5, 26.5, 29.5, 32, 29, 28, 28,
                                                   28.5, 28, 25, 28.5, 28.5, 30,
                                                   32.5, 29, 26.5, 30.5, 27))),
                       .Names = c("cut","tool","bevel","mil_del"),
                       class = "data.frame",
                       row.names = c(NA, -32L))
Anova_T <- EMSanova(mil_del~tool*bevel*cut,data=power_cut,type=c("F","F","F")) # nested=c(NA,NA,NA)
library(ggplot2,plyr)
desc_power <- ddply(power_cut,.(cut,bevel),summarise,pw_mean = mean(mil_del),
                    pw_ci = 1.96*sd(mil_del)/sqrt(length(mil_del)))
pd <- position_dodge(width = 0.2)
(ansPlot <-  ggplot(data = desc_power,aes(x=bevel,y=pw_mean,group = cut))+
    geom_line(aes(linetype = cut), position = pd) +
    geom_errorbar(aes(ymin= pw_mean - pw_ci, ymax= pw_mean + pw_ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Type of Cut")) +
    labs(title = c(),x = c("Bevel (degrees)"),
         y = c("Sample mean")) +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

#####################
# # Chapter 2 # #
#################### 
#' Figure about
#' the \beta distribution
par(mfrow=c(3,1)) ; 
n = 300
x1 <- rnorm(n,mean = 19.5, sd=0.4)
h1 <- hist(x1,breaks = 1, xlim = c(17,21),ylim=c(0,700),
           border = "white",main="",xlab = expression(mu == 19.5),ylab = "") 
x1fit <- seq(min(x1),max(x1),length=length(x1)) ; abline(v=19.5,lwd=3) ; abline(v=18.842)
y1fit <- dnorm(x1fit,mean=mean(x1),sd=sd(x1)) ; y1fit <- y1fit*diff(h1$mids[1:2])*n
lines(x1fit,y1fit,col="blue", lwd=2)
text(18,150,expression(alpha == 0.05)) ; text(21,100,expression(bar(y)))

x2 <- rnorm(n,mean = 19, sd=0.4)
h2 <- hist(x2,breaks = 1, xlim = c(17,21),ylim=c(0,700),
           border = "white",main="",xlab = expression(mu == 19.0),ylab = "") 
x2fit <- seq(min(x2),max(x2),length=length(x2)) ; abline(v=19.0,lwd=3) ; abline(v=18.842)
y2fit <- dnorm(x2fit,mean=mean(x2),sd=sd(x2)) ; y2fit <- y2fit*diff(h2$mids[1:2])*n
lines(x2fit,y2fit,col="blue", lwd=2)
text(20,600,expression(beta)) ; text(21,100,expression(bar(y)))

x3 <- rnorm(n,mean = 18.5, sd=0.4)
h3 <- hist(x3,breaks = 1, xlim = c(17,21),ylim=c(0,700),
           border = "white",main="",xlab = expression(mu == 18.5),ylab = "") 
x3fit <- seq(min(x3),max(x3),length=length(x3)) ; abline(v=18.5,lwd=3) ; abline(v=18.842)
y3fit <- dnorm(x3fit,mean=mean(x3),sd=sd(x3)) ; y3fit <- y3fit*diff(h3$mids[1:2])*n
lines(x3fit,y3fit,col="blue", lwd=2)
text(19.5,400,expression(beta)) ; text(21,100,expression(bar(y)))

x4 <- rnorm(n,mean = 18, sd=0.4)
h4 <- hist(x4,breaks = 1, xlim = c(17,21),ylim=c(0,700),
           border = "white",main="",xlab = expression(mu == 18.0),ylab = "") 
x4fit <- seq(min(x4),max(x4),length=length(x4)) ; abline(v=18.0,lwd=3) ; abline(v=18.842)
y4fit <- dnorm(x4fit,mean=mean(x4),sd=sd(x4)) ; y4fit <- y4fit*diff(h4$mids[1:2])*n
lines(x4fit,y4fit,col="blue", lwd=2)
text(19.0,200,expression(beta)) ; text(21,100,expression(bar(y)))

# # Example 2.1 # #

# # Problems # #
y  <- 18472.90
mu <- 18470
var<- 40 
n  <- 50
z  <- (y-mu)/(var/sqrt(n)) 
z2 <- (sqrt(n)*(y-mu))/var

tmp1 <- var/sqrt(n)
tmp2 <- mu-1.645*(tmp1) # the desirable value for stating that is significantly different

dta <- data.frame(matrix(data = c(seq(18450,18470,1)),nrow = 21,ncol = 1));colnames(dta) <- c('u')
for(i in 1:nrow(dta)){
  dta$z1[i]       <- (tmp2-dta$u[i])/tmp1
  dta$one_beta[i] <- pnorm(dta$z1[i])
  dta$beta[i]     <- 1-dta$one_beta[i]
}

(ansPlot <- ggplot(data = dta,aes(x=u, y=beta)) +
    labs(title = 'Operating characteristic curve', x = 'mu', y = 'beta') +
    geom_line() +
    theme_bw(base_size = 12)+
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2)) 
)

y1 <- 42.2
y2 <- 44.5
v1 <- 10
v2 <- 18
n1 <- 8
n2 <- 15
yd <- y1-y2
vd <- v1-v2
t <- yd/(vd/sqrt(n1+n2))
f <- v1/v2

# # 2.11 # #
stu <- structure(list(before  = as.numeric(c(14,12,20,8,11,15,17,18,9,7)),
                      after   = as.numeric(c(17,16,21,10,10,14,20,22,14,12))),
                 .Names = c("before","after"),class="data.frame",
                 row.names= c(NA,-10L))
stu$dif <- stu$after-stu$before
tt <- t.test(stu$before,stu$after,paired=TRUE)

# # 2.13 # #
y  <- as.numeric(c(12,8,14,20,26,26,20,21,18,24,30,21,18,16,10,20))
y1 <- mean(y)
mu1 <- 12
mu2 <- 16
mu3 <- 18
var <- 7 
n   <- length(y) 
z  <- (y1-mu2)/(var/sqrt(n)) 
z2 <- (sqrt(n)*(y1-mu2))/var

# # 2.14 # #
y  <- as.numeric(c(12,8,14,20,26,26,20,21,18,24,30,21,18,16,10,20))
y  <- mean(y)
mu <- seq(12,26,1)
var<- 7
n  <- length(y)
dt <- data.frame(matrix(data = mu,nrow = length(mu),ncol = 1));colnames(dt) <- c('mu')
for(i in 1:nrow(dt)){
  dt$z[i] <- (y-mu[i])/(var/sqrt(n))
  dt$one_beta[i] <- pnorm(dt$z[i])
  dt$beta[i]     <- 1-dt$one_beta[i]
}
j   <- c(13,14,15,18,20,22,23)
dt2 <- subset(dt,dt$mu %in% j)
(ansPlot <- ggplot(data = dt,aes(x=mu, y=one_beta)) +
    labs(title = 'Power', x = 'mu', y = '1-beta(mu)') +
    geom_line() +
    theme_bw(base_size = 12)+
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
    scale_x_continuous(limits = c(dt$mu[1],dt$mu[15]),breaks = seq(12,26,2))
)

# # 2.16 # #
dt <- c(8.2,8.3,8.2,8.6,8.0,8.4,8.5,7.3)
sigma_fx <- 0.35
s2 <- sd(dt)^2
w <- (length(dt)-1)*(s2/sigma_fx)
val1 <- (length(dt)-1)*s2/14.067
val2 <- (length(dt)-1)*s2/2.167
round(sqrt(val1),digits = 1)
round(sqrt(val2),digits = 1)

# # 2.19 # #
exp <- c(9.333,4.945,12)
ctr <- c(8.375,1.187,8)
f <- sqrt(exp[2])/sqrt(ctr[2])
z <- exp[1]-ctr[1]/sqrt(((exp[2]^2)/exp[3])+((ctr[2]^2)/ctr[3]))

exp2 <- rnorm(exp[3],mean=exp[1],sd=exp[2])
ctr2 <- rnorm(ctr[3],mean=ctr[1],sd=ctr[2])
t.test(x=exp2,y=ctr2)

#'compute t estimate
par1 <- exp[1]-ctr[1]
par2 <- ((((exp[3]-1))*exp[2]^2)+(((ctr[3]-1))*ctr[2]^2))/(exp[3]+ctr[3]-2)
par3 <- (1/exp[3])+(1/ctr[3])
t <- par1/sqrt(par2*par3)

# # 2.20 # #
boys <- c(101,86,72,129,99,118,104,125,90,107) # values of WISC for boys
girls <- c(97,107,94,101,90,108,108,86)        # values of WISC for girls
t.test(x=boys,y=girls)
#'optionally
wisc <- structure(list(sex    = structure(c(rep(1L,10),rep(2L,8)),.Label = c("Boy","Girl"),class = "factor"),
                       scores = as.numeric(c(101,86,72,129,99,118,104,125,90,107,97,107,94,101,90,108,108,86))),
                  .Names = c("sex","wisc_scores"),
                  class = "data.frame",row.names = c(NA, -18L))

# # 2.24 # #
pop1 <- c(16,5,9) # values of population 1: [1] bar{y_1}, [2] s_1^2, [3] n_1
pop2 <- c(12,2,4) # values of population 2: [1] bar{y_2}, [2] s_2^2, [3] n_2
#' First run an F test to see if variances are unequal, if so no test on means may be necessary
#' as the two samples are so heterogeneous with respect to their variances
#' Apply the case B samples are independent and variances are known. First conduct an F test for
#' equal variances...
f <- pop1[2]/pop2[2] # F=S_1^2/S_2^2 = Following table D with v_1 = n_1 - 1 & v_2 = n_2 - 1
#' equal of variances can be assumed

tmp1 <- pop1[1]-pop2[1]# \bar{Y_1}-\bar{Y_2}
numT <- ((pop1[3]-1)*pop1[2])+((pop2[3]-1)*pop2[2]) # (S_1^2/n_1)+(S_2^2/n_2)
denT <- pop1[3]+pop2[3]-2
tmp2 <- sqrt((numT/denT)*((1/pop1[3])+(1/pop2[3])))
t <- tmp1/tmp2
df <- pop1[3]+pop2[3]-2

# # 2.28 # #
dta <- data.frame(matrix(data = c(0.1776,0.1830,0.2635,0.1995,0.1954,
                                  0.2142,0.2039,0.2343,0.2247,0.1656,
                                  0.2428,0.1710,0.2022,0.1623,0.1735,
                                  0.1607,0.1969,0.1465,0.1443,0.1582,
                                  0.1616,0.1628,0.1509,0.2202,0.2303,
                                  0.1684,0.1719,0.1623,0.1755,0.2015,
                                  0.1860,0.1842,0.1990,0.1770,0.1814,
                                  0.2921,0.2261,0.1992,0.2080,0.2718,
                                  0.2042,0.2762,0.1683,0.2414,0.2386,
                                  0.2072,0.2618,0.1948,0.2351,0.1771,
                                  0.1678,0.2260,0.2212,0.2241,0.2282,
                                  0.1661,0.2114,0.2036,0.2244,0.1855,
                                  0.2311,0.1952,0.1795,0.2380,0.2171,
                                  0.2171,0.1812,0.2072,0.1731,0.2556,
                                  0.1644,0.1923,0.2200,0.1672,0.2892,
                                  0.1840,0.1762,0.1713,0.1698,0.1786,
                                  0.1891,0.2238,0.1852,0.1709,0.1977,
                                  0.2041,0.1644,0.1880,0.1801,0.2245,
                                  0.1993,0.1726,0.1666,0.1621,0.1882,
                                  0.2109,0.1827,0.2369,0.1959,0.2059,
                                  0.1879,0.2097,0.2680,0.1905,0.2232,
                                  0.2912,0.1881,0.1659,0.1682,0.2086,
                                  0.1736,0.1616,0.1940,0.1442,0.2297,
                                  0.2021,0.2186,0.1968,0.1765,0.1579,
                                  0.1883,0.1986,0.1585,0.1719,0.1688
),nrow = 25,ncol = 5, byrow = TRUE))
colnames(dta)[1:ncol(dta)] <- c(paste0('x',seq(1:ncol(dta))))
for(i in 1:nrow(dta)){
  dta$x_mean[i] <- mean(dta$x1[i]:dta$x5[i])
}
# Add a Normal Curve (Thanks to Peter Dalgaard)
x <- dta$x_mean
h<-hist(x, breaks=10, col="gray", xlab="residual binder", 
        main="Histogram with Normal Curve") 
xfit<-seq(min(x),max(x),length=40) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)

t <- ((mean(dta$x_mean)-0.1985)/(sd(dta$x_mean)/sqrt(nrow(dta)-1)))
y <- mean(dta$x_mean)
tmp1 <- 2.80*(sd(dta$x_mean)/sqrt(nrow(dta)-1))

# # 2.30 # #
y <- dta$x_mean
sd_y <- sd(y)
mn_y <- mean(y)
p1 <- subset(dta,dta$x_mean >= mn_y-sd_y & dta$x_mean <= mn_y+sd_y,select = c(6))
p1 <- nrow(p1)/nrow(dta)
p2 <- subset(dta,dta$x_mean >= mn_y-(sd_y)*2 & dta$x_mean <= mn_y+(sd_y)*2,select = c(6))
p2 <- nrow(p2)/nrow(dta)
p3 <- subset(dta,dta$x_mean >= mn_y-(sd_y)*3 & dta$x_mean <= mn_y+(sd_y)*3,select = c(6))
p3 <- nrow(p3)/nrow(dta)

#####################
# # # Chapter 3 # # #
#####################
dta <- iris
j   <- ddply(dta, .(Species), summarise, SUM=sum(Sepal.Length),MEAN=mean(Sepal.Length),n=length(Sepal.Length))
t_av<- mean(dta$Sepal.Length)
for(i in 1:nrow(j)){
  j$residual[i] <- (j$MEAN[i]-t_av)^2
}

res <- summary(aov(Sepal.Length ~ Species,data = dta))
#' The F value is then calculated from the MeanSq/Residual MeanSq

(ansPlot <-  ggplot(data = j,aes(x=Species, y=MEAN,group=3)) +
    labs(title = 'Sepal length', x = element_blank(), y = 'length') +
    geom_point(aes(colour=Species),na.rm = TRUE)+
    scale_y_continuous(limits = c(3,8),breaks = seq(3,8,1)) 
)
dev.off()

# Matrix for the orthogonal contrast
dta_2 <- data.frame(matrix(data=c(1,rep(-1/3,3),
                                  0,1,rep(-1/2,2),
                                  0,0,1,-1),nrow = 3,ncol = 4,byrow = TRUE));

# # Example 3.2 # #
dta_aov <- structure(list(fact = structure(c(rep(1L,4),rep(2L,4),rep(3L,4),rep(4L,4)), .Label = c("Factor1","Factor2","Factor3","Factor4"), class = "factor"),
                          mil_def = as.numeric(c(1.93,2.38,2.20,2.25,2.55,2.72,2.75,2.70,2.40,2.68,2.31,2.28,2.33,2.40,2.28,2.25))),
                     .Names = c("fact","vals"), class = "data.frame",
                     row.names = c(NA, -20L)) 

j   <- ddply(dta_aov, .(fact), summarise, t_total=sum(vals),mean=mean(vals),n=length(vals))
t_av<- mean(dta_aov$vals)
for(i in 1:nrow(dta_aov)){ # computing the sum of squares
  #tmp <- subset(j,j$fact == dta_aov$fact[i])
  dta_aov$squares[i] <- (dta_aov$vals[i])^2
}
j2 <- ddply(dta_aov,.(fact),summarise,ssq=sum(squares))

j3 <- cbind(j,j2[,2])

# # 3.1
dta_aov <- structure(list(fact = structure(c(rep(1L,5),rep(2L,5),rep(3L,5),rep(4L,5),rep(5L,5)), .Label = c("Factor1","Factor2","Factor3","Factor4","Factor5"), class = "factor"),
                          mil_def = as.numeric(c(8,6,7,5,8,4,-2,0,-2,3,1,2,0,-1,-3,4,6,5,5,4,10,8,7,4,9))),
                     .Names = c("fact","vals"), class = "data.frame",
                     row.names = c(NA, -25L)) 

j   <- ddply(dta_aov, .(fact), summarise, t_total=sum(vals),fact_mean=mean(vals),n=length(vals))
t_av<- mean(dta_aov$vals)
for(i in 1:nrow(dta_aov)){ # computing the sum of squares
  #tmp <- subset(j,j$fact == dta_aov$fact[i])
  dta_aov$squares[i] <- (dta_aov$vals[i])^2
}
j2 <- ddply(dta_aov,.(fact),summarise,ssq=sum(squares))

j3 <- cbind(j,j2[,2])
colnames(j3)[5] <- c('ssquares')

ss_total <- sum(j3$ssquares)-((sum(j3$t_total)^2)/sum(j3$n))
ss_treat <- j3$t_total[1]^2/5+j3$t_total[2]^2/5+j3$t_total[3]^2/5+j3$t_total[4]^2/5+j3$t_total[5]^2/5-sum(j3$t_total)^2/25
ss_error <- ss_total-ss_treat
res <- summary(aov(vals ~ fact,data = dta_aov))

# # 3.6 # #
dta_aov <- structure(list(fact = structure(c(rep(1L,2),rep(2L,6),rep(3L,11),rep(4L,4),rep(5L,2)), .Label = c("ang_67","ang_71","ang_75","ang_79","ang_83"), class = "factor"),
                          vals = as.numeric(c(83,85,84,85,85,86,86,87,86,87,87,87,88,88,88,88,88,89,90,89,90,90,91,90,92))),
                     .Names = c("fact","vals"), class = "data.frame",
                     row.names = c(NA, -25L)) 

j   <- ddply(dta_aov, .(fact), summarise, t_total=sum(vals),fact_mean=mean(vals),n=length(vals))
t_av<- mean(dta_aov$vals)
for(i in 1:nrow(dta_aov)){ # computing the sum of squares
  #tmp <- subset(j,j$fact == dta_aov$fact[i])
  dta_aov$squares[i] <- (dta_aov$vals[i])^2
}
j2 <- ddply(dta_aov,.(fact),summarise,ssq=sum(squares))

j3 <- cbind(j,j2[,2])
colnames(j3)[5] <- c('ssquares')

ss_total <- sum(j3$ssquares)-((sum(j3$t_total)^2)/sum(j3$n))
ss_treat <- j3$t_total[1]^2/2+j3$t_total[2]^2/6+j3$t_total[3]^2/11+j3$t_total[4]^2/4+j3$t_total[5]^2/2-sum(j3$t_total)^2/25
ss_error <- ss_total-ss_treat
res <- summary(aov(vals ~ fact,data = dta_aov))

# # 3.7 # #
ver1 <- j3$n[4]*j3$t_total[1]-j3$n[1]*j3$t_total[4]
ver2 <- j3$n[5]*j3$t_total[2]-j3$n[2]*j3$t_total[5]
ver3 <- j3$n[3]*j3$t_total[1]+j3$n[3]*j3$t_total[2]-sum(j3$n[1],j3$n[2],j3$n[4],j3$n[5])*j3$t_total[3]+j3$n[3]*j3$t_total[4]+j3$n[3]*j3$t_total[5]
ver4 <- sum(j3$n[2],j3$n[5])*j3$t_total[1]-sum(j3$n[1],j3$n[4])*j3$t_total[2]+sum(j3$n[2],j3$n[5])*j3$t_total[4]-sum(j3$n[1],j3$n[4])*j3$t_total[5]

ss1 <- (ver1)^2/(2*4^2+4*(-2)^2)
ss2 <- (ver2)^2/(6*2^2+2*(-6)^2)
ss3 <- (ver3)^2/(2*11^2+6*11^2+11*(-14)^2+4*11^2+2*11^2)
ss4 <- (ver4)^2/(2*8^2+6*(-6)^2+4*8^2+2*(-6)^2)

# # 3.15 # #
dta_aov <- structure(list(fact = structure(c(rep(1L,5),rep(2L,5),rep(3L,5),rep(4L,5)), .Label = c("f_1","f_2","f_3","f_4","f_5"), class = "factor"),
                          vals = as.numeric(c(56,55,62,59,60,64,61,50,55,56,45,46,45,39,43,42,39,45,43,41))),
                     .Names = c("fact","vals"), class = "data.frame",
                     row.names = c(NA, -20L)) 

j   <- ddply(dta_aov, .(fact), summarise, t_total=sum(vals),fact_mean=mean(vals),n=length(vals))
t_av<- mean(dta_aov$vals)
for(i in 1:nrow(dta_aov)){ # computing the sum of squares
  #tmp <- subset(j,j$fact == dta_aov$fact[i])
  dta_aov$squares[i] <- (dta_aov$vals[i])^2
}
j2 <- ddply(dta_aov,.(fact),summarise,ssq=sum(squares))

j3 <- cbind(j,j2[,2])
colnames(j3)[5] <- c('ssquares')

ss_total <- sum(j3$ssquares)-((sum(j3$t_total)^2)/sum(j3$n))
ss_treat <- j3$t_total[1]^2/5+j3$t_total[2]^2/5+j3$t_total[3]^2/5+j3$t_total[4]^2/5-sum(j3$t_total)^2/20
ss_error <- ss_total-ss_treat
res <- summary(aov(vals ~ fact,data = dta_aov))

#####################
# # # Chapter 5 # # #
#####################
# # 5.1
phosphor_table <- structure(list(phos_type = structure(c(rep(1L,6),rep(2L,6),rep(3L,6)), .Label = c("A","B","C"), class = "factor"),
                                 glass_type = structure(c(rep(c(rep(1L,3),rep(2L,3)),3)), .Label = c("I","II"), class = "factor"),
                                 light = as.numeric(c(280,290,285,230,235,240,
                                                      300,310,295,260,240,235,
                                                      270,285,290,220,225,230))),
                            .Names = c("phos_type","glass_type","light"), class = "data.frame",
                            row.names = c(NA, -18L)) # last number in row.names is just nrow()

# # 5.8
thrust_table <- structure(list(speed = structure(c(rep(1L,12),rep(2L,12),rep(3L,12),rep(4L,12),rep(5L,12)), .Label = c("100","220","475","715","870"), class = "factor"),
                               mat   = structure(c(rep(c(rep(1L,6),rep(2L,6)),5)), .Label = c("B10","V10"), class = "factor"),
                               feed  = structure(c(rep(c(rep(1L,2),rep(2L,2),rep(3L,2)),10)), .Label = c("0.004","0.008","0.014"), class = "factor"),
                               force = as.numeric(c(122,110,332,330,640,500,192,170,386,365,810,725,
                                                    108,85,276,310,612,500,136,130,333,330,779,670,
                                                    108,60,248,295,543,450,122,85,318,330,810,750,
                                                    66,50,248,275,612,610,108,75,472,350,893,890,
                                                    80,60,276,310,696,610,136,75,499,390,1820,890))),
                          .Names = c("speed","mat","feed","force"), class = "data.frame",
                          row.names = c(NA, -60L)) # last number in row.names is just nrow()

# # 5.18
job_table <- structure(list(analyst = structure(c(rep(1L,8),rep(2L,8),rep(3L,8),rep(4L,8),rep(5L,8)), .Label = c("Ana1","Ana2","Ana3","Ana4","Ana5"), class = "factor"),
                            dimension  = structure(c(rep(c(rep(1L,4),rep(2L,4)),5)), .Label = c("2D","3D"), class = "factor"),
                            job  = structure(c(rep(c(rep(1L,2),rep(2L,2)),10)), .Label = c("Job1","Job2"), class = "factor"),
                            rating = as.numeric(c(1.42,1.25,1.40,1.44,1.00,0.90,0.92,0.93,
                                                  1.59,2.09,1.83,2.02,1.46,1.68,1.43,1.02,
                                                  1.26,1.48,1.80,1.75,0.85,1.43,1.33,1.48,
                                                  1.76,1.47,0.98,1.23,1.21,1.25,0.73,1.22,
                                                  1.43,1.72,1.17,1.18,1.60,2.03,1.31,1.12))),
                       .Names = c("analyst","dimension","job","rating"), class = "data.frame",
                       row.names = c(NA, -40L)) # last number in row.names is just nrow()

# # 5.26
new_table <- structure(list(bonder = structure(c(rep(1L,15),rep(2L,15),rep(3L,15),rep(4L,15)), .Label = c("A","B","C","D"), class = "factor"),
                            position  = structure(c(rep(c(rep(1L,5),rep(2L,5),rep(3L,5)),4)), .Label = c("Pos1","Pos2","Pos3"), class = "factor"),
                            pres = as.numeric(c(204,181,201,203,214,262,246,230,256,288,220,232,235,231,220,
                                                197,223,206,232,213,207,259,223,195,197,214,248,191,197,222,
                                                264,226,228,249,246,255,186,237,236,240,215,176,171,208,180,
                                                248,138,273,220,186,304,330,268,295,276,208,248,247,220,241))),
                       .Names = c("bonder","position","pres"), class = "data.frame",
                       row.names = c(NA, -60L)) # last number in row.names is just nrow()
# # 5.29
heat_sag <- structure(list(time_delay = structure(c(rep(1L,36),rep(2L,36),rep(3L,36)), .Label = c("0","12","24"), class = "factor"),
                           postcure_time  = structure(c(rep(c(rep(1L,12),rep(2L,12),rep(3L,12)),3)), .Label = c("60","90","120"), class = "factor"),
                           heat = as.numeric(c(0.55,0.64,0.77,0.73,0.82,0.63,0.93,0.91,1.08,0.92,1.07,0.73,
                                               0.67,0.67,0.76,0.90,0.87,0.84,0.54,0.75,0.73,0.80,0.56,0.66,
                                               0.47,0.68,0.71,0.54,0.67,0.55,0.69,0.59,0.65,0.72,0.62,0.53,
                                               0.71,0.90,0.78,0.81,1.22,0.92,1.04,0.79,1.02,1.02,0.76,1.04,
                                               0.86,0.74,0.83,0.82,0.72,0.61,0.79,0.92,0.82,0.85,0.63,0.67,
                                               0.86,0.67,0.76,0.53,0.62,0.88,0.65,0.72,0.88,0.77,0.47,0.70,
                                               0.90,0.67,0.65,0.92,0.81,0.72,0.68,0.75,0.66,0.67,0.82,0.81,
                                               0.95,0.77,1.35,1.44,0.75,0.78,0.56,1.03,0.99,0.47,0.68,0.98,
                                               0.53,0.47,0.19,0.63,0.46,rep(c(NA),7)))),
                      .Names = c("time_delay","postcure_time","heat"), class = "data.frame",
                      row.names = c(NA, -108L)) # last number in row.names is just nrow()

# # 5.30
electronic_db <- structure(list(machine = structure(c(rep(1L,6),rep(2L,6),rep(3L,6)), .Label = c("I","II","III"), class = "factor"),
                                manufacturer  = structure(c(rep(c(rep(1L,3),rep(2L,3)),3)), .Label = c("Man1","Man2"), class = "factor"),
                                bond = as.numeric(c(13.5,14.0,14.5,13.0,12.0,11.0,
                                                    17.0,16.0,15.0,15.5,14.5,15.0,
                                                    18.5,18.0,17.5,10.0,12.0,13.0))),
                           .Names = c("machine","manufacturer","bond"), class = "data.frame",
                           row.names = c(NA, -18L)) # last number in row.names is just nrow()

#####################
# # # Chapter 6 # # #
#####################
test_db <- structure(list(tester = structure(c(rep(1L,4),rep(2L,4),rep(3L,4),rep(4L,4)), .Label = c("I","II","III","IV"), class = "factor"),
                          cards  = structure(c(rep(c(rep(1L,2),rep(2L,2)),4)), .Label = c("Card1","Card2"), class = "factor"),
                          current = as.numeric(c(40.498,40.513,40.382,40.375,
                                                 40.235,40.253,40.158,40.131,
                                                 40.252,40.277,40.185,40.215,
                                                 40.458,40.503,40.460,40.507))),
                     .Names = c("tester","cards","current"), class = "data.frame",
                     row.names = c(NA, -16L)) # last number in row.names is just nrow()

# I like this package
fit <- anovaVCA(current~tester*cards, test_db)

# # 6.18
test2_db <- structure(list(operator = structure(c(rep(1L,20),rep(2L,20),rep(3L,20)), .Label = c("I","II","III"), class = "factor"),
                           part  = structure(c(rep(c(rep(1L,2),rep(2L,2),rep(3L,2),rep(4L,2),rep(5L,2),rep(6L,2),rep(7L,2),rep(8L,2),rep(9L,2),rep(10L,2)),3)), .Label = c(LETTERS[seq(from=1,to=10)]), class = "factor"),
                           current = as.numeric(c(0.417,0.419,0.417,0.417,0.423,0.418,0.412,0.410,0.407,0.409,
                                                  0.408,0.413,0.409,0.408,0.408,0.411,0.412,0.408,0.410,0.404,
                                                  0.394,0.398,0.387,0.399,0.389,0.407,0.389,0.405,0.386,0.405,
                                                  0.382,0.400,0.385,0.400,0.384,0.401,0.387,0.401,0.386,0.407,
                                                  0.404,0.401,0.398,0.402,0.407,0.402,0.407,0.411,0.400,0.410,
                                                  0.405,0.410,0.407,0.400,0.402,0.405,0.412,0.405,0.418,0.404))),
                      .Names = c("operator","part","current"), class = "data.frame",
                      row.names = c(NA, -60L)) # last number in row.names is just nrow()

#####################
# # # Chapter 7 # # #
#####################
library('EMSaov')
library('VCA')
## Example for data in Table 7.2
strain_db <- structure(list(machine = structure(c(rep(1L,16),rep(2L,16),rep(3L,16),rep(4L,16),rep(5L,16)), .Label = c(LETTERS[seq(from=1,to=5)]), class = "factor"),
                            head  = structure(c(rep(c(rep(1L,4),rep(2L,4),rep(3L,4),rep(4L,4)),5)), .Label = c(letters[seq(from=1,to=4)]), class = "factor"),
                            strain = as.numeric(c(6,2,0,8,13,3,9,8,1,10,0,6,7,4,7,9,
                                                  10,9,7,12,2,1,1,10,4,1,7,9,0,3,4,1,
                                                  0,0,5,5,10,11,6,7,8,5,0,7,7,2,5,4,
                                                  11,0,6,4,5,10,8,3,1,8,9,4,0,8,6,5,
                                                  1,4,7,9,6,7,0,3,3,0,2,2,3,7,4,0))),
                       .Names = c("machine","head","strain"), class = "data.frame",
                       row.names = c(NA, -80L)) # last number in row.names is just nrow()

fit1 <- anovaVCA(strain~machine*head, strain_db) #VCA
fit2 <- EMSanova(strain~machine*head,data=strain_db,type=c("F","R"),nested=c(NA,"machine"))
desc1 <- ddply(strain_db,.(machine,head),summarise,m=sum(strain))
desc2 <- ddply(strain_db,.(machine),summarise,m=sum(strain))

# # Example 7.1
guns_db <- structure(list(group = structure(c(rep(1L,12),rep(2L,12),rep(3L,12)), .Label = c("I","II","III"), class = "factor"),
                          method  = structure(c(rep(c(rep(1L,2),rep(2L,2)),9)), .Label = c("Method1","Method2"), class = "factor"),
                          team  = as.numeric(rep(c(c(rep(1L,4),rep(2L,4),rep(3L,4))),3)),
                          gun = as.numeric(c(20.2,24.1,14.2,16.2,26.2,26.9,18.0,19.1,23.8,24.9,12.5,15.4,
                                             22.0,23.5,14.1,16.1,22.6,24.6,14.0,18.1,22.9,25.0,13.7,16.0,
                                             23.1,22.9,14.1,16.1,22.9,23.7,12.2,13.8,21.8,23.5,12.7,15.1))),
                     .Names = c("group","method","team","gun"), class = "data.frame",
                     row.names = c(NA, -36L)) # last number in row.names is just nrow()

guns_db$team <- as.factor(guns_db$team)

library('EMSaov')
fit <- EMSanova(gun~group*method*team,data=guns_db,type=c("F","F","R"),nested=c(NA,NA,"group"))

## Example 7.3
training2_db <- structure(list(group = structure(c(rep(1L,14),rep(2L,14),rep(3L,14)), .Label = c("I","II","III"), class = "factor"),
                               subject  = structure(c(rep(seq(1,7,1),6)), .Label = c(letters[seq(from=1,to=7)]), class = "factor"),
                               test  = structure(c(rep(c(rep(1L,7),rep(2L,7)),3)), .Label = c('Pretest','Postest'), class = "factor"),
                               velocity = as.numeric(c(seq(1,42,1)))),
                          .Names = c("group","subject","test","velocity"), class = "data.frame",
                          row.names = c(NA, -42L)) # last number in row.names is just nrow()

fit <- EMSanova(velocity~group*subject*test,data=training2_db,type=c("F","R","F"),nested=c(NA,"group",NA))



# # PROBLEM 7.1
porosity_db <- structure(list(lot = structure(c(rep(1L,12),rep(2L,12),rep(3L,12)), .Label = c("I","II","III"), class = "factor"),
                              roll = as.numeric(rep(c(c(rep(1L,3),rep(2L,3),rep(3L,3),rep(4L,3))),3)),
                              porosity = as.numeric(c(1.5,1.7,1.6,1.5,1.6,1.7,2.7,1.9,2.0,3.0,2.4,2.6,
                                                      1.9,1.5,2.1,2.3,2.4,2.4,1.8,2.9,4.7,1.9,3.5,2.8,
                                                      2.5,2.9,3.3,3.2,5.5,7.1,1.4,1.5,3.4,7.8,5.2,5.0))),
                         .Names = c("lot","roll","porosity"), class = "data.frame",
                         row.names = c(NA, -36L)) # last number in row.names is just nrow()
porosity_db$roll <- as.factor(porosity_db$roll)

fit <- EMSanova(porosity~lot*roll,data=porosity_db,type=c("R","R"),nested=c(NA,"lot"))
desc <- ddply(porosity_db,.(lot,roll),summarise, suma=sum(porosity))
desc1 <- ddply(porosity_db,.(lot),summarise, suma=sum(porosity))
pf(q=, df1=1, df2=10, lower.tail=FALSE)

# # PROBLEM 7.3
pseudo_db <- structure(list(A = structure(c(rep(c(seq(1,5,1)),24)), .Label = c(LETTERS[seq(from=1,to=5)]), class = "factor"),
                            B = structure(c(rep(c(seq(1,4,1)),30)), .Label = c(letters[seq(from=1,to=4)]), class = "factor"),
                            C = structure(c(rep(c(seq(1,3,1)),40)), .Label = c('I','II','III'), class = "factor"),
                            n = as.numeric(c(seq(1,120,1)))),
                       .Names = c("A","B","C","n"), class = "data.frame",
                       row.names = c(NA, -120L)) # last number in row.names is just nrow()
pseudo_db <- arrange(pseudo_db,A,B,C)
pseudofit <- EMSanova(n~A*B*C,data=pseudo_db,type=c("F","R","R"),nested=c(NA,"A","B"))

# # PROBLEM 7.10
methods_db <- structure(list(method = structure(c(rep(1L,30),rep(2L,30)), .Label = c("New","Traditional"), class = "factor"),
                             subject  = structure(c(rep(seq(1,10,1),6)), .Label = c(letters[seq(from=1,to=10)]), class = "factor"),
                             test  = structure(c(rep(c(rep(1L,10),rep(2L,10),rep(3L,10)),2)), .Label = c('Pretest','Postest','Aftertest'), class = "factor"),
                             grade = as.numeric(c(rnorm(10,mean=60,sd=9.3),rnorm(10,mean=60,sd=9.5),rnorm(10,mean=60,sd=10),
                                                  rnorm(10,mean=50,sd=10),rnorm(10,mean=67,sd=9),rnorm(10,mean=68,sd=11)))),
                        .Names = c("method","subject","test","grade"), class = "data.frame",
                        row.names = c(NA, -60L)) # last number in row.names is just nrow()

fit <- EMSanova(grade~method*subject*test,data=methods_db,type=c("F","R","F"),nested=c(NA,"method",NA))
pd <- position_dodge(width = 0.2)
methods_desc <- ddply(methods_db,.(test,method),summarise,score_mean=mean(grade),score_ci=1.96*sd(grade)/sqrt(length(grade)))
methods_desc <- arrange(methods_desc,method)
(ansPlot <-  ggplot(data = methods_desc,aes(x=test,y=score_mean,group = method))+
    geom_line(aes(linetype = method), position = pd) +
    geom_errorbar(aes(ymin= score_mean - score_ci, ymax= score_mean + score_ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Methods")) +
    labs(title = paste("Mean score depending on",
                       "teaching method.",
                       "Error bars represent 95% Confidence Intervals",
                       sep = "\n"),
         x = "Test",
         y = "Achievement score") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

# #7.14
creative_db <- structure(list(method = structure(c(rep(1L,4),rep(2L,4),rep(3L,4)), .Label = c("Restricted","Permisive","Control"), class = "factor"),
                              class  = structure(c(rep(seq(1,4,1),3)), .Label = c(letters[seq(from=1,to=4)]), class = "factor"),
                              cog_style  = structure(c(rep(c(rep(1L,3),rep(2L,3)),2)), .Label = c('field-dependent','field-independent'), class = "factor"),
                              grade = as.numeric(c(rnorm(12,mean=60,sd=9)))),
                         .Names = c("method","class","cog_style","grade"), class = "data.frame",
                         row.names = c(NA, -12L)) # last number in row.names is just nrow()

fit <- EMSanova(grade~method*class*cog_style,data=creative_db,type=c("F","R","F"),nested=c(NA,"method",NA))
pf(q=15.97, df1=2, df2=9, lower.tail=FALSE)

# #7.19
xmeasure_db <- structure(list(day = structure(c(rep(c(seq(1,3,1)),22)), .Label = c('I','II','III'), class = "factor"),
                              grader = structure(c(rep(c(seq(1,2,1)),33)), .Label = c('A','B'), class = "factor"),
                              sample = structure(c(rep(c(seq(1,11,1)),6)), .Label = c(letters[seq(from=1,to=11)]), class = "factor")),
                         .Names = c("day","grader","sample"), class = "data.frame",
                         row.names = c(NA, -66L)) # last number in row.names is just nrow()
xmeasure_db$measure <- as.numeric(c(4,6,6,13,7,7,14,12,9,6,8,
                                    11,7,10,11,10,11,16,10,12,9,13,
                                    5,17,8,3,14,11,6,11,16,-1,3,
                                    11,13,15,14,20,19,11,17,4,9,14,
                                    0,-1,2,8,8,4,5,10,16,8,7,
                                    6,-2,5,2,6,10,18,13,17,15,11))
desc1 <- ddply(xmeasure_db,.(day),summarise,suma=sum(measure))
desc2 <- ddply(xmeasure_db,.(day,grader),summarise,suma=sum(measure))
fit2 <- EMSanova(measure~day*grader*sample,data=xmeasure_db,type=c("R","F","R"),nested=c(NA,NA,"day"))

# #7.24
resistance_db <- structure(list(filler = structure(c(rep(1L,54),rep(2L,54),rep(3L,54),rep(4L,54),rep(5L,54)), .Label = c("Alumina","Graphite","Iron Filings","Iron Oxide","Copper"), class = "factor"),
                                ratio  = structure(c(rep(c(rep(1L,27),rep(2L,27)),5)), .Label = c("0.5:1","1:1"), class = "factor"),
                                sample = structure(c(rep(c(rep(1L,9),rep(2L,9),rep(3L,9)),10)), .Label = c(letters[seq(from=1,to=3)]), class = "factor"),
                                plate  = structure(c(rep(c(rep(1L,3),rep(2L,3),rep(3L,3)),30)), .Label = c("I","II","III"), class = "factor"),
                                resistance = as.numeric(c(7,7.3,6.6,5.4,5.1,5,6.4,5.7,7.7,
                                                          6.8,6.9,6.2,5.1,5.5,4.8,4.6,5.8,5.8,
                                                          7.3,7.6,5.6,6.3,5.3,4.8,5.1,6.5,7.1,
                                                          7.1,6.9,6.5,5,5.7,4,6.4,6.6,6.2,
                                                          6.8,7.2,4.8,4,5.2,4.4,7.5,7.8,6,
                                                          6.2,6.6,4.4,4.4,4.8,3.8,5.7,6.6,5.6,
                                                          8.8,8.3,6.9,5.5,5.9,4.7,9,11,11.3,
                                                          7.7,6.5,5.9,6.8,7,7.2,10.1,11.2,10.8,
                                                          7.2,6.8,6.2,7,6.9,5.2,7.7,9.3,9.8,
                                                          18,17.6,18.4,14.4,12.4,12.6,18.5,19.2,19.4,
                                                          15.6,14.6,15.1,11.6,12.8,13.4,17.6,18.7,19.3,
                                                          13.3,13.8,13.6,11.7,13.4,12.1,15.5,18,19.6,
                                                          2.1,2.1,1,1.1,1.1,0.9,1.3,1.7,1.7,
                                                          1.7,1.8,1.3,1.1,1,0.8,1.7,2,2,
                                                          3.2,2.8,2,0.8,0.8,0.5,1.4,2,1.9,
                                                          3.6,3.8,2.9,1.6,1,1.6,2.3,2.5,2.4,
                                                          2.1,2.6,1.6,1.1,1.1,0.9,1.7,1,1.4,
                                                          2.3,2.6,2.1,1.1,1.2,0.7,1.8,2,2.4,
                                                          3,3.8,3,1.1,1.8,1.1,3.2,4.1,3.3,
                                                          3.1,3,2.7,1.6,1.7,1.2,2.2,3,2.7,
                                                          2.1,2.9,1.6,1.5,1.9,2.55,2,2.5,2.4,
                                                          1.8,2.1,1.9,1.3,1,1.6,2.4,1.9,1.8,
                                                          2.1,2,2.3,1,0.8,1,1.5,2,2,
                                                          1.6,1.6,1.9,1.2,1.1,1.2,1.6,1.4,2.1,
                                                          2.7,2.9,3,1.4,2.2,1.8,2.8,3.8,3.4,
                                                          2.5,3,2,1.8,2.6,2.1,2.4,3.3,2.4,
                                                          3.3,3.2,3.2,2.1,2,2,3.4,4.1,3.8,
                                                          2.8,2.6,1.9,1.5,1.1,1.4,2.4,2.1,2,
                                                          2.1,2.5,1.5,1,1,1.2,1.5,1.6,1.7,
                                                          2.6,2.9,2.8,1.3,1.5,1.2,2.6,3.3,2.7))),
                           .Names = c("filler","ratio","sample","plate","resistance"), class = "data.frame",
                           row.names = c(NA, -270L)) # last number in row.names is just nrow()
fit3 <- EMSanova(resistance~filler*ratio*sample*plate,data=resistance_db,type=c("F","F","R","F"),nested=c(NA,NA,'filler',NA))
fit4 <- EMSanova(resistance~filler*ratio*sample*plate,data=resistance_db,type=c("F","F","R","F"),nested=c(NA,NA,c(1,2),NA))
pf(q=5.4252615, df1=20, df2=180, lower.tail=FALSE)
desc_3 <- ddply(resistance_db,.(filler,ratio),summarise, m=mean(resistance))
res_desc <- ddply(resistance_db,.(filler,ratio),summarise,wear_mean=mean(resistance),score_ci=1.96*sd(resistance)/sqrt(length(resistance)))

pd <- position_dodge(width = 0.2)
(ansPlot <-  ggplot(data = res_desc,aes(x=filler,y=wear_mean,group=ratio))+
    geom_line(aes(linetype = ratio), position = pd) +
    geom_errorbar(aes(ymin= wear_mean - score_ci, ymax= wear_mean + score_ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    geom_hline(yintercept = 3,linetype = "dashed",colour='blue') +
    guides(linetype = guide_legend("Ratio")) +
    labs(title = paste("Mean wear depending on",
                       "fillers and Ratio Concentration.",
                       "Error bars represent 95% Confidence Intervals",
                       sep = "\n"),
         x = "Fillers",
         y = "Wear resistance") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

# # 7.32
y_db <- structure(list(repet = structure(c(rep(1L,27),rep(2L,27),rep(3L,27)), .Label = c("I","II","III"), class = "factor"),
                       temp  = structure(c(rep(c(rep(1L,9),rep(2L,9),rep(3L,9)),3)), .Label = c("temp1","temp2","temp3"), class = "factor"),
                       acid  = structure(c(rep(c(rep(1L,3),rep(2L,3),rep(3L,3)),9)), .Label = c("acid1","acid2","acid3"), class = "factor"),
                       loc = structure(c(rep(c(seq(1,3)),27)), .Label = c("loc1","loc2","loc3"), class = "factor"),
                       y = as.numeric(c(1.18,1.48,2.20,-0.45,0.94,1.58,-1.02,0.36,1.06,
                                        1.67,1.91,1.99,-0.21,0.36,0.81,0.36,0.67,0.94,
                                        1.06,1.58,1.83,0.52,0.52,1.48,0.81,1.28,1.48,
                                        0.67,0.81,1.58,0.81,1.28,1.39,-0.21,0.94,1.28,
                                        1.28,1.83,2.06,0.52,0.52,0.52,0.19,0.52,1.06,
                                        1.48,1.75,1.99,0.19,0.19,0.36,1.06,1.28,1.75,
                                        0.81,1.83,2.13,0.81,1.28,1.67,0.00,0.36,0.81,
                                        1.75,1.83,2.13,0.36,0.52,0.81,0.52,0.52,0.81,
                                        1.18,1.39,1.83,-0.21,0.52,0.67,1.06,1.18,1.28))),
                  .Names = c("repet","temp","acid","location","y"), class = "data.frame",
                  row.names = c(NA, -81L)) # last number in row.names is just nrow()
fity <- EMSanova(y~temp*acid*location*repet,data=y_db,type=c("F","F","F","R"),nested=c(NA,NA,"acid",NA))
desc_3 <- ddply(y_db,.(temp,acid),summarise, m=mean(y))

#test
fity <- EMSanova(y~temp*acid*location,data=y_db,type=c("F","F","F"),nested=c(NA,NA,NA))
#####################
# # # Chapter 8 # # #
#####################
library(EMSaov)
library(VCA)
library(lme4)
## Example 8.1
nozzle_db <- structure(list(repe = structure(c(rep(c(1,2,3),15)), .Label = c("I","II","III"), class = "factor"),
                            operator = structure(c(rep(1L,9),rep(2L,9),rep(3L,9),rep(4L,9),rep(5L,9)), .Label = c(letters[seq(from=1,to=5)]), class = "factor"),
                            nozzle = structure(c(rep(c(rep(1L,3),rep(2L,3),rep(3L,3)),5)), .Label = c(LETTERS[seq(from=1, to=3)]), class = "factor"),
                            flow = as.numeric(c(6,6,-15,13,6,13,10,10,-11,26,12,5,4,4,11,-35,0,-14,11,4,4,17,10,17,11,-10,-17,
                                                21,14,7,-5,2,-5,12,-2,-16,25,18,25,15,8,1,-4,10,24))),
                       .Names = c("repe","operator","nozzle","flow"), class = "data.frame",
                       row.names = c(NA, -45L)) # last number in row.names is just nrow()

fit_noz <- EMSanova(flow~repe*operator*nozzle,data=nozzle_db,type=c("R","F","F")) # replicated the textbook result but you have to sum up the sources that contain the factor repe
#fit_noze <- anovaVCA(flow~repe*operator*nozzle,nozzle_db) #VCA failed to replicate the textbook
fitest_noz <- EMSanova(flow~operator*nozzle*repe,data=nozzle_db,type=c("R","F","F"))
fitest_noz <- aov(flow~(operator*nozzle),data=nozzle_db)


## Problem8.1
glass_db <- structure(list(station = structure(c(rep(1L,27),rep(2L,27)), .Label = c("I","II"), class = "factor"),
                           week = structure(c(rep(c(rep(1L,9),rep(2L,9),rep(3L,9)),2)), .Label = c("week1","week2","week3"), class = "factor"),
                           shift = structure(c(rep(c(rep(1L,3),rep(2L,3),rep(3L,3)),6)), .Label = c("1","2","3"), class = "factor"),
                           rating = as.numeric(c(3,6,6,3,4,7,3,6,7,
                                                 14,16,19,8,8,9,11,12,17,
                                                 2,3,6,2,3,4,2,4,6,
                                                 6,8,11,3,9,11,6,8,13,
                                                 4,6,7,15,15,17,4,7,10,
                                                 2,5,7,2,4,6,10,12,13))),
                      .Names = c("station","week","shift","rating"), class = "data.frame",
                      row.names = c(NA, -54L)) # last number in row.names is just nrow()

fit_glass <- EMSanova(rating~station*shift*week,data=glass_db,type=c("F","F","R"))
pf(q=1.2384, df1=1, df2=46, lower.tail=FALSE)


## Problem 8.5
alum_db <- structure(list(machine = structure(c(rep(1L,6),rep(2L,6),rep(3L,6),rep(4L,6),rep(5L,6)), .Label = c("I","II","III","IV","V"), class = "factor"),
                          rep = structure(c(rep(c(rep(1L,2),rep(2L,2),rep(3L,2)),5)), .Label = c("rep1","rep2","rep3"), class = "factor"),
                          slug = structure(c(rep(c(rep(1L,1),rep(2L,1)),15)), .Label = c("Rivet","Staple"), class = "factor"),
                          thickness = as.numeric(c(175,165,190,170,185,190,
                                                   95,165,185,160,165,160,
                                                   180,175,180,175,175,200,
                                                   170,185,200,165,195,185,
                                                   155,130,190,190,200,200))),
                     .Names = c("machine","rep","slug","thickness"), class = "data.frame",
                     row.names = c(NA, -30L)) # last number in row.names is just nrow()

fit_alum <- EMSanova(thickness~rep*machine*slug,data=alum_db,type=c("R","F","F"))
pf(q=5.299, df1=2, df2=8, lower.tail=FALSE)

## Problem 8.10
cereal_db <- structure(list(run = structure(c(rep(1L,36),rep(2L,36)), .Label = c("I","II"), class = "factor"),
                            position = structure(c(rep(c(rep(1L,12),rep(2L,12),rep(3L,12)),2)), .Label = c("B","M","T"), class = "factor"),
                            cycle = structure(c(rep(c(rep(1L,6),rep(2L,6)),6)), .Label = c("C1","C2"), class = "factor"),
                            head = structure(c(rep(c(seq(1,6,1)),12)), .Label = c(letters[seq(from=1,to=6)]),class = "factor"),
                            filling_capacity = as.numeric(c(19,13,13,13,19,6,
                                                            rep(0,4),19,0,
                                                            31,19,31,19,19,0,
                                                            25,38,19,31,38,0,
                                                            25,31,13,25,19,0,
                                                            13,25,0,13,25,0,
                                                            6,0,28,28,16,13,
                                                            6,13,6,0,13,6,
                                                            0,0,25,6,13,6,
                                                            6,31,19,0,19,19,
                                                            0,6,13,28,0,0,
                                                            0,0,28,25,0,19))),
                       .Names = c("run","position","cycle","head","filling_capacity"), class = "data.frame",
                       row.names = c(NA, -72L)) # last number in row.names is just nrow()
cereal_db <- cereal_db[,c(2,4,3,1,5)]
lol <- count(cereal_db,.(position,head));lol$seq <- c(seq(1,nrow(lol),1));lol$seq <- as.factor(lol$seq)
cereal_db$combi <- c(NA)
for(i in 1:nrow(cereal_db)){
  x <- subset(lol,lol$position == cereal_db$position[i] & lol$head == cereal_db$head[i],select = c('seq'))
  cereal_db$combi[i] <- letters[as.numeric(x[1])]
  rm(x)
}
cereal_db$combi <- as.factor(cereal_db$combi)

fit_cereal_1 <- EMSanova(filling_capacity~run*combi,data=cereal_db,type=c("R","F"),nested = c(NA,NA))

fit_cereal_2 <- EMSanova(filling_capacity~run*position*head,data=cereal_db,type=c("R","F","F"),nested = c(NA,NA,NA))

## Problem 8.17
ash_db <- structure(list(repli = structure(c(rep(1L,20),rep(2L,20)), .Label = c("I","II"), class = "factor"),
                         test = structure(c(rep(c(rep(1L,10),rep(2L,10)),2)), .Label = c("supplier","Plant"), class = "factor"),
                         sample = structure(c(rep(c(seq(1,10,1)),4)), .Label = c(letters[seq(from=1,to=10)]),class = "factor"),
                         ash = as.numeric(c(5.43,5.35,5.53,5.55,5.93,5.95,6.34,6.05,5.87,5.60,
                                            5.21,5.48,5.48,5.61,5.96,6.09,6.49,6.12,6.04,5.70,
                                            5.47,5.31,5.46,5.55,5.93,5.97,6.32,6.02,5.87,5.58,
                                            5.13,5.46,5.54,5.54,6.00,5.99,6.43,6.13,6.00,5.67))),
                    .Names = c("repli","test","sample","ash"), class = "data.frame",
                    row.names = c(NA, -40L)) # last number in row.names is just nrow()

fit_ash <- EMSanova(ash~repli+sample+test,data=ash_db,type=c("R","R","F"),nested = c(NA,"repli",NA))

pf(q=3.02014304481098, df1=1, df2=19, lower.tail=FALSE)
pd <- position_dodge(width = 0.2)
ash_desc <- ddply(ash_db,.(repli,test),summarise,ash_reading=mean(ash),score_ci=2.576*sd(ash)/sqrt(length(ash)))
ash_desc2 <- ddply(ash_db,.(test),summarise,ash_reading=mean(ash),score_ci=2.576*sd(ash)/sqrt(length(ash)))

(ansPlot <-  ggplot(data = ash_desc,aes(x=repli,y=ash_reading,group = test))+
    geom_line(aes(linetype = test), position = pd) +
    geom_errorbar(aes(ymin= ash_reading - score_ci, ymax= ash_reading + score_ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Laboratories")) +
    labs(title = paste("Mean score depending on",
                       "Laboratory readings.",
                       "Error bars represent 99% Confidence Intervals",
                       sep = "\n"),
         x = "Replicates",
         y = "Mean Ash Readings") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

(ansPlot <-  ggplot(data = ash_desc2,aes(x=test,y=ash_reading))+
    geom_errorbar(aes(ymin= ash_reading - score_ci, ymax= ash_reading + score_ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    labs(title = paste("Mean score depending on",
                       "Laboratory readings.",
                       "Error bars represent 99% Confidence Intervals",
                       sep = "\n"),
         x = "Laboratories",
         y = "Mean Ash Readings") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

#####################
# # # Chapter 9 # # #
#####################
library(EMSaov)

## Example 9.1
c_b <- (-1)*(210)+(-1)*(240)+(+1)*(180)+(+1)*(200)
ss <- (c_b)^2/4# total sum of square

## Example 9.2
nfactor_db <- structure(list(factor1 = structure(c(rep(c(rep(1L,5),rep(2L,5)),2)), .Label = c('La','Ha'), class = "factor"),
                             factor2 = structure(c(rep(1L,10),rep(2L,10)), .Label = c('Lb','Hb'), class = "factor"),
                             measure  = as.numeric(c(125,50,40,0,50,
                                                     175,150,50,100,175,
                                                     150,250,200,240,325,
                                                     250,275,200,350,400))),
                        .Names = c("factor1","factor2","measure"), class = "data.frame",
                        row.names = c(NA, -20L)) # last number in row.names is just nrow()
desc_nfactor <- ddply(nfactor_db,.(factor1,factor2),summarise,suma=sum(measure))
ss_ap <- (nrow(nfactor_db)-1)*((sd(nfactor_db$measure))^2) # approximated sum of squares
fit <- EMSanova(measure~factor1*factor2,data=nfactor_db,type=c("F","F"),nested=c(NA,NA))

##Example 9.3

##Problem 9.1
xfactor_db <- structure(list(factor1 = structure(c(rep(1L,6),rep(2L,6)), .Label = c('A1','A2'), class = "factor"),
                             factor2 = structure(c(rep(c(rep(1L,3),rep(2L,3)),2)), .Label = c('B1','B2'), class = "factor"),
                             measure  = as.numeric(c(0,2,1,-1,-3,1,
                                                     4,6,2,-1,-3,-7))),
                        .Names = c("factor1","factor2","measure"), class = "data.frame",
                        row.names = c(NA, -12L)) # last number in row.names is just nrow()

fit_x <- EMSanova(measure~factor1*factor2,data=xfactor_db,type=c("F","F"),nested=c(NA,NA))

##Problem 9.4
chem_db <- structure(list(factor1 = structure(c(rep(1L,16),rep(2L,16)), .Label = c('A1','A2'), class = "factor"),
                          factor2 = structure(c(rep(c(rep(1L,8),rep(2L,8)),2)), .Label = c('B1','B2'), class = "factor"),
                          factor3 = structure(c(rep(c(rep(1L,4),rep(2L,4)),4)), .Label = c('C1','C2'), class = "factor"),
                          factor4 = structure(c(rep(c(rep(1L,2),rep(2L,2)),8)), .Label = c('D1','D2'), class = "factor"),
                          measure  = as.numeric(c(1985,1592,2156,2032,1694,1712,2184,1921,
                                                  1765,1700,1923,2007,1806,1758,1957,1717,
                                                  1595,2067,1578,1733,2243,1745,1745,1818,
                                                  1835,1823,1863,1910,1614,1838,1917,1922))),
                     .Names = c("factor1","factor2","factor3","factor4","measure"), class = "data.frame",
                     row.names = c(NA, -32L)) # last number in row.names is just nrow()
chem_x <- EMSanova(measure~factor1*factor2*factor3*factor4,data=chem_db,type=c("F","F","F","F"),nested=c(NA,NA,NA,NA))

#graph
eff <- c(-41.4375,-27.8125,52.5625,1.6875,53.0625,-38.8125,-50.9375,100.6875,-134.9375,33.9375,125.8125,-4.3125,29.0625,-6.0625,49.3125)
qqnorm(eff);qqline(eff,col=2)
q <- quantile(eff)
y  <- qqnorm(eff)
db <- structure(list(eff=c(eff),
                     x=c(as.numeric(y$x[]))),
                .Names =c("effect","theoretical_quantile"),class = "data.frame",
                row.names=c(NA,-15L))
db$source <- c('A','B','AB','C','AC','BC','ABC','D','AD','BD','ABD','CD','ACD','BCD','ABCD')
ggplot(db,aes(theoretical_quantile,effect,label=source))+
  labs(title = '               Normal Q-Q Plot', x = 'Theoretical Quantiles', y = 'Effect') +
  geom_text()+
  theme_bw()

##Problem 9.17
steel_db <- structure(list(factor1 = structure(c(rep(1L,4),rep(2L,4)), .Label = c('A1','A2'), class = "factor"),
                           factor2 = structure(c(rep(c(rep(1L,2),rep(2L,2)),2)), .Label = c('B1','B2'), class = "factor"),
                           factor3 = structure(c(rep(c(rep(1L,1),rep(2L,1)),4)), .Label = c('C1','C2'), class = "factor"),
                           measure  = as.numeric(c(7.8,6.8,4.3,5.9,4.5,11.0,8.2,6.0))),
                      .Names = c("factor1","factor2","factor3","measure"), class = "data.frame",
                      row.names = c(NA, -8L)) # last number in row.names is just nrow()

steel_x <- EMSanova(measure~factor1*factor2*factor3,data=steel_db,type=c("F","F","F"),nested=c(NA,NA,NA))
pf(q=10.348241301059, df1=1, df2=3, lower.tail=FALSE)


desc_steel <- ddply(steel_db,.(factor1,factor2,factor3),summarise,suma=sum(measure))
eff <- c(1.225,-1.425,1.225,0.775,0.925,-1.525,-2.825)
q1 <- quantile(eff)
qqnorm(eff);qqline(eff,col=2)
y <- qqnorm(eff)
db <- structure(list(eff=c(eff),
                     x=c(as.numeric(y$x[]))),
                .Names =c("effect","theoretical_quantile"),class = "data.frame",
                row.names=c(NA,-7L))
db$source <- c('A','B','C','AB','AC','BC','ABC')
ggplot(db,aes(theoretical_quantile,effect,label=source))+
  labs(title = '                Normal Q-Q Plot', x = 'Theoretical Quantiles', y = 'Effect') +
  geom_text()+
  theme_bw()

steel_db$combi <- c(NA)
director <- count(steel_db,.(factor2,factor3))
for(i in 1:nrow(steel_db)){
  if(steel_db$factor2[i] == director$factor2[1] & steel_db$factor3[i] == director$factor3[1]){
    steel_db$combi[i] <- c('(-1,-1)')
  }else if(steel_db$factor2[i] == director$factor2[2] & steel_db$factor3[i] == director$factor3[2]){
    steel_db$combi[i] <- c('(-1,1)')
  }else if(steel_db$factor2[i] == director$factor2[3] & steel_db$factor3[i] == director$factor3[3]){
    steel_db$combi[i] <- c('(1,-1)')
  }else if(steel_db$factor2[i] == director$factor2[4] & steel_db$factor3[i] == director$factor3[4]){
    steel_db$combi[i] <- c('(1,1)')
  }
}
steel_db$combi <- as.factor(steel_db$combi)

pd <- position_dodge(width = 0.2)
(ansPlot <-  ggplot(data = steel_db,aes(x=combi,y=measure,group = factor1))+
    geom_line(aes(linetype = factor1), position = pd) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Factor A")) +
    labs(title = paste("Interaction Plot",
                       "between three factors.",
                       sep = "\n"),
         x = "(B/C)",
         y = "Heat") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

##Problem 9.19
steel2_db <- structure(list(factor1 = structure(c(rep(1L,4),rep(2L,4)), .Label = c('A1','A2'), class = "factor"),
                            factor2 = structure(c(rep(c(rep(1L,2),rep(2L,2)),2)), .Label = c('B1','B2'), class = "factor"),
                            factor3 = structure(c(rep(c(rep(1L,1),rep(2L,1)),4)), .Label = c('C1','C2'), class = "factor"),
                            measure  = as.numeric(c(0,0,4.3,4.4,9.8,10.1,4,3.9))),
                       .Names = c("factor1","factor2","factor3","measure"), class = "data.frame",
                       row.names = c(NA, -8L)) # last number in row.names is just nrow()
steel2_x <- EMSanova(measure~factor1*factor2*factor3,data=steel2_db,type=c("F","F","F"),nested=c(NA,NA,NA))


steel2_db$combi <- c(NA)
director <- count(steel_db,.(factor1,factor2))
for(i in 1:nrow(steel2_db)){
  if(steel2_db$factor1[i] == director$factor1[1] & steel2_db$factor2[i] == director$factor2[1]){
    steel2_db$combi[i] <- c('(-1,-1)')
  }else if(steel2_db$factor1[i] == director$factor1[2] & steel2_db$factor2[i] == director$factor2[2]){
    steel2_db$combi[i] <- c('(-1,1)')
  }else if(steel2_db$factor1[i] == director$factor1[3] & steel2_db$factor2[i] == director$factor2[3]){
    steel2_db$combi[i] <- c('(1,-1)')
  }else if(steel2_db$factor1[i] == director$factor1[4] & steel2_db$factor2[i] == director$factor2[4]){
    steel2_db$combi[i] <- c('(1,1)')
  }
}
steel2_db$combi <- as.factor(steel2_db$combi)
model <- aov(measure~combi,data = steel2_db)
out <- SNK.test(model, 'combi', alpha = 0.05, group=TRUE, main = NULL,console=FALSE)
out$means$combi <- c('(-1,-1)','(-1,1)','(1,-1)','(1,1)')
out$means$combi <- as.factor(out$means$combi)
out$means$ci <- 1.96*out$means$std/sqrt(2)

pd <- position_dodge(width = 0.2)
(ansPlot <-  ggplot(data = out$means,aes(x=combi,y=measure))+
    geom_errorbar(aes(ymin= measure - ci, ymax= measure + ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    labs(title = paste("            Mean heat of two factors interaction",
                       sep = "\n"),
         x = "(A/B)",
         y = "Heat") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)
##quantile plot
desc_steel <- ddply(steel2_db,.(factor1,factor2,factor3),summarise,suma=sum(measure))
eff <- c(4.775,-0.825,0.075,-5.175,0.025,-0.075,-0.125)
q1 <- quantile(eff)
qqnorm(eff);qqline(eff,col=2)
y <- qqnorm(eff)
db <- structure(list(eff=c(eff),
                     x=c(as.numeric(y$x[]))),
                .Names =c("effect","theoretical_quantile"),class = "data.frame",
                row.names=c(NA,-7L))
db$source <- c('A','B','C','AB','AC','BC','ABC')
ggplot(db,aes(theoretical_quantile,effect,label=source))+
  labs(title = '           Normal Q-Q Plot', x = 'Theoretical Quantiles', y = 'Effect') +
  geom_text()+
  theme_bw()

#Problem 9.22
chip_db <- structure(list(lot = structure(c(rep(1L,4),rep(2L,4)), .Label = c('A','B'), class = "factor"),
                          temperature = structure(c(rep(c(rep(1L,2),rep(2L,2)),2)), .Label = c('Cold','Room'), class = "factor"),
                          rate  = as.numeric(c(3548,3932,10107,9562,6704,6656,34653,31097))),
                     .Names = c("lot","temperature","rate"), class = "data.frame",
                     row.names = c(NA, -8L)) # last number in row.names is just nrow()

chip_x <- EMSanova(rate~lot*temperature,data=chip_db,type=c("F","F"),nested=c(NA,NA))
chip_x <- EMSanova(rate~temperature*lot,data=chip_db,type=c("F","R"),nested=c(NA,'temperature'))

##Problem 9.30
library('olsrr')
bonding_db <- structure(list(country = structure(c(rep(1L,10),rep(2L,10)), .Label = c('Domestic','Foreign'), class = "factor"),
                             design = structure(c(rep(c(rep(1L,5),rep(2L,5)),2)), .Label = c('Standard','New'), class = "factor"),
                             grams  = as.numeric(c(204,181,201,203,214,
                                                   264,226,228,249,246,
                                                   244,270,208,235,216,
                                                   197,223,206,232,213))),
                        .Names = c("country","design","grams"), class = "data.frame",
                        row.names = c(NA, -20L)) # last number in row.names is just nrow()
bond_x <- EMSanova(grams~country*design,data=bonding_db,type=c("F","F"),nested=c(NA,NA))
bond_desc <- ddply(bonding_db,.(country,design),summarise, strength_mean=mean(grams),ci=1.96*sd(grams)/sqrt(length(grams)))

pd <- position_dodge(width = 0.2)

(ansPlot <-  ggplot(data = bond_desc,aes(x=design,y=strength_mean,group = country))+
    geom_line(aes(linetype = country), position = pd) +
    geom_errorbar(aes(ymin= strength_mean - ci, ymax= strength_mean + ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Country")) +
    labs(title = paste("Mean bond strength depending on",
                       "Machine design and country.",
                       "Error bars represent 95% Confidence Intervals",
                       sep = "\n"),
         x = "Design",
         y = "Grams") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

qqnorm(bonding_db$grams);qqline(bonding_db$grams,col=2)

bond.lm = lm(grams ~ country+design, data=bonding_db) 
bond.res = resid(bond.lm)

plot(bonding_db$grams, bond.res, ylab="Residuals", xlab="Grams", main="Residual vs Fitted values") 
abline(0, 0) 

ols_norm_test(bond.lm)

pf(q=0.147437128, df1=1, df2=24, lower.tail=FALSE)

##problem 32
totals_db <- structure(list(power = structure(c(rep(1L,4),rep(2L,4)), .Label = c('A1','A2'), class = "factor"),
                            SF6 = structure(c(rep(c(rep(1L,2),rep(2L,2)),2)), .Label = c('B1','B2'), class = "factor"),
                            OxyR = structure(c(rep(c(rep(1L,1),rep(2L,1)),4)), .Label = c('C1','C2'), class = "factor"),
                            height  = as.numeric(c(8.45,10.11,11.18,12.42,7.85,9.19,14.76,16.64))),
                       .Names = c("power","SF6","OxyR","height"), class = "data.frame",
                       row.names = c(NA, -8L)) # last number in row.names is just nrow()
##test
totals_db$combi <- c(NA)
for(i in 1:nrow(totals_db)){
  if(totals_db$SF6[i] == 'B1' & totals_db$OxyR[i] == 'C1'){
    totals_db$combi[i] <- c('(-1,-1)')
  }else if(totals_db$SF6[i] == 'B1' & totals_db$OxyR[i] == 'C2'){
    totals_db$combi[i] <- c('(-1,1)')
  }else if(totals_db$SF6[i] == 'B2' & totals_db$OxyR[i] == 'C1'){
    totals_db$combi[i] <- c('(1,-1)')
  }else if(totals_db$SF6[i] == 'B2' & totals_db$OxyR[i] == 'C2'){
    totals_db$combi[i] <- c('(1,1)')
  }
}
totals_db$combi <- as.factor(totals_db$combi)
#model <- aov(measure~combi,data = steel2_db)
#out <- SNK.test(model, 'combi', alpha = 0.05, group=TRUE, main = NULL,console=FALSE)
#out$means$combi <- c('(-1,-1)','(-1,1)','(1,-1)','(1,1)')
#out$means$combi <- as.factor(out$means$combi)
#out$means$ci <- 1.96*out$means$std/sqrt(2)

pd <- position_dodge(width = 0.2)
(ansPlot <-  ggplot(data = totals_db,aes(x=combi,y=height,group=power))+
    geom_line(aes(linetype = power), position = pd) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    labs(title = paste("    Total height of a trench in a 2^3 factorial",
                       "           factorial experiment.",
                       sep = "\n"),
         x = "(B/C)",
         y = "Height") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)
##quantile plot
eff <- c(0.3825,1.2125,0.3925,0.0075,0.02,0.5825,0.06)
q1 <- quantile(eff)
qqnorm(eff);qqline(eff,col=2)
y <- qqnorm(eff)
db <- structure(list(eff=c(eff),
                     x=c(as.numeric(y$x[]))),
                .Names =c("effect","theoretical_quantile"),class = "data.frame",
                row.names=c(NA,-7L))
db$source <- c('A','B','C','AB','AC','BC','ABC')
ggplot(db,aes(theoretical_quantile,effect,label=source))+
  labs(title = '       Normal Q-Q Plot', x = 'Theoretical Quantiles', y = 'Effect') +
  geom_text()+
  theme_bw()


######################
# # # Chapter 10 # # #
######################
# # Example 9.2
sample <- c(125,175,150,250,50,150,250,275,40,50,200,200,0,100,240,350,50,175,325,400)
s_y <- sd(sample)^2
SS_total <- (length(sample)-1)*s_y
## Problem 10.1
pull_db <- structure(list(temperature = structure(c(rep(1L,6),rep(2L,6),rep(3L,6)), .Label = c('Cold','Ambient','Hot'), class = "factor"),
                          humidity = structure(c(rep(c(rep(1L,2),rep(2L,2),rep(3L,2)),3)), .Label = c('50%','70%','90%'), class = "factor"),
                          measure  = as.numeric(c(0.8,2.8,1.0,1.6,2.0,2.2,
                                                  1.5,3.2,1.6,1.8,1.5,0.8,
                                                  2.5,4.2,1.8,1.0,2.5,4.0))),
                     .Names = c("temperature","humidity","measure"), class = "data.frame",
                     row.names = c(NA, -18L)) # last number in row.names is just nrow()

fit <- EMSanova(measure~temperature*humidity,data=pull_db,type=c("F","F"),nested=c(NA,NA))
desc_pull <- ddply(pull_db,.(temperature,humidity),summarise,suma=sum(measure),n=length(measure))

## Problem 10.2
orthogonal_db <- structure(list(treatment = as.character(c('00','01','02','10','11','12','20','21','22')),
                                A_L = as.numeric(c(-1,-1,-1,0,0,0,1,1,1)),
                                A_Q = as.numeric(c(1,1,1,-2,-2,-2,1,1,1)),
                                B_L = as.numeric(c(-1,0,1,-1,0,1,-1,0,1)),
                                B_Q = as.numeric(c(1,-2,1,1,-2,1,1,-2,1)),
                                A_L_B_L = as.numeric(c(1,0,-1,0,0,0,-1,0,1)),
                                A_L_B_Q = as.numeric(c(-1,+2,-1,0,0,0,1,-2,1)),
                                A_Q_B_L = as.numeric(c(-1,0,1,2,0,-2,-1,0,1)),
                                A_Q_B_Q = as.numeric(c(1,-2,1,-2,4,-2,1,-2,1)),
                                totals  = as.numeric(c(desc_pull$suma[]))),
                           .Names = c("treatment","A_L","A_Q","B_L","B_Q","A_L_B_L","A_L_B_Q","A_Q_B_L","A_Q_B_Q","totals"), class = "data.frame",
                           row.names = c(NA, -9L)) # last number in row.names is just nrow()
sum_c1sq <- data.frame(matrix(data = NA,nrow = 8,ncol = 3));colnames(sum_c1sq)[1:3] <- c('treatment','sumc1sq','latex')
sum_c1sq$treatment <- colnames(orthogonal_db)[2:9]
for(i in 1:nrow(sum_c1sq)){
  sum_c1sq$sumc1sq[i] <- sum(orthogonal_db[,i+1]^2)
  for(j in 1:nrow(orthogonal_db)){
    if(j == 1){
      tempa <- c(paste0(as.character(colnames(orthogonal_db)[i+1])," Sigma c_i^2= "))
      tempa <- c(paste0(tempa,'(',orthogonal_db[j,i+1],')^2+'))
    }else if(j > 1 & j < nrow(orthogonal_db)){
      tempa <- c(paste0(tempa,'(',orthogonal_db[j,i+1],')^2+'))
    }else if(j == nrow(orthogonal_db)){
      tempa <- c(paste0(tempa,'(',orthogonal_db[j,i+1],')^2'))
      sum_c1sq$latex[i] <- as.character(tempa)
    }
  }
}
sum_c1sq$latex <- paste0(sum_c1sq$latex[],'=',sum_c1sq$sumc1sq)

contrast <- data.frame(matrix(data = NA,nrow = 8,ncol = 3));colnames(contrast)[1:3] <- c('treatment','contrast','latex')
contrast$treatment <- colnames(orthogonal_db)[2:9]
for(i in 2:9){
  tempo <- orthogonal_db[,i]*orthogonal_db[,10] 
  contrast$contrast[i-1] <- sum(tempo)
  for(j in 1:nrow(orthogonal_db)){
    if(j == 1){
      tempa <- c(paste0("C_{",as.character(colnames(orthogonal_db)[i]),"}="))
      tempa <- c(paste0(tempa,'(',orthogonal_db[j,i],')(',orthogonal_db[j,10],')+'))
    }else if(j > 1 & j < nrow(orthogonal_db)){
      tempa <- c(paste0(tempa,'(',orthogonal_db[j,i],')(',orthogonal_db[j,10],')+'))
    }else if(j == nrow(orthogonal_db)){
      tempa <- c(paste0(tempa,'(',orthogonal_db[j,i],')(',orthogonal_db[j,10],')'))
      contrast$latex[i-1] <- as.character(tempa)
    }
  }
}
contrast$latex <- paste0(contrast$latex[],'=',round(contrast$contrast[],2))

ss <- data.frame(matrix(data = NA,nrow = 8,ncol = 3));colnames(ss)[1:3] <- c('treatment','ss','latex')
ss$treatment <- colnames(orthogonal_db)[2:9]
n <- 2 # define n before, number of observations per treatment combination
for(i in 1:nrow(ss)){
  ss$ss[i] <- round(contrast$contrast[i]^2/(n*sum_c1sq$sumc1sq[i]),3)
  ss$latex[i] <- c(paste0('SS_{',as.character(colnames(orthogonal_db)[i+1]),'}= frac{(',round(contrast$contrast[i],3),')^2}{',n,'(',sum_c1sq$sumc1sq[i],')','}'))
}
ss$latex <- paste0(ss$latex[],'=',ss$ss[])

## Problem 10.4


## Problem 10.5
noxious_db <- structure(list(volume = structure(c(rep(1L,6),rep(2L,6),rep(3L,6)), .Label = c('Low','Medium','High'), class = "factor"),
                             time = structure(c(rep(c(rep(1L,2),rep(2L,2),rep(3L,2)),3)), .Label = c('Short','Medium','Long'), class = "factor"),
                             measure  = as.numeric(c(0.3,0.4,0.1,0.4,0.8,2.0,
                                                     0.1,0.4,0.1,0.2,0.7,1.6,
                                                     0.4,0.4,0.4,1.2,2.5,3.6))),
                        .Names = c("volume","time","measure"), class = "data.frame",
                        row.names = c(NA, -18L)) # last number in row.names is just nrow()

fit <- EMSanova(measure~volume*time,data=noxious_db,type=c("F","F"),nested=c(NA,NA))

##Problem 10.11
concrete_db <- structure(list(surface = structure(c(rep(1L,18),rep(2L,18),rep(3L,18)), .Label = c('Surface3','Surface4','Surface5'), class = "factor"),
                              base = structure(c(rep(c(rep(1L,6),rep(2L,6),rep(3L,6)),3)), .Label = c('Base0','Base3','Base6'), class = "factor"),
                              subbase = structure(c(rep(c(rep(1L,2),rep(2L,2),rep(3L,2)),9)), .Label = c('Sub4','Sub8','Sub12'), class = "factor"),
                              measure  = as.numeric(c(2.8,2.6,4.1,4.4,5.5,5.3,
                                                      4.3,4.5,5.7,5.8,7.0,6.8,
                                                      5.7,5.3,6.9,7.1,8.1,8.3,
                                                      4.1,4.4,5.3,5.1,6.5,6.7,
                                                      5.4,5.5,6.5,6.7,7.7,7.5,
                                                      6.7,6.9,7.7,7.4,8.8,9.1,
                                                      6.0,6.2,6.1,5.8,7.0,7.2,
                                                      6.3,6.5,7.2,7.1,8.0,8.3,
                                                      7.1,6.9,8.1,8.4,9.1,9.0))),
                         .Names = c("surface","base","subbase","measure"), class = "data.frame",
                         row.names = c(NA, -54L)) # last number in row.names is just nrow()

fit <- EMSanova(measure~surface*base*subbase,data=concrete_db,type=c("F","F","F"),nested=c(NA,NA,NA))
desc_concreto0 <- ddply(concrete_db,.(surface,base),summarise,suma=sum(measure))
desc_concreto1 <- ddply(concrete_db,.(surface,subbase),summarise,suma=sum(measure))
desc_concreto2 <- ddply(concrete_db,.(base,subbase),summarise,suma=sum(measure))
desc_concreto3 <- ddply(concrete_db,.(surface,base,subbase),summarise,suma=sum(measure),n=length(measure))

## Problem 10.12
orthogonal_2 <- structure(list(treatment = as.character(c('000','001','002','010','011','012','020','021','022',
                                                          '100','101','102','110','111','112','120','121','122',
                                                          '200','201','202','210','211','212','220','221','222')),
                               A_L = as.numeric(c(rep(-1,9),rep(0,9),rep(1,9))),
                               A_Q = as.numeric(c(rep(1,9),rep(-2,9),rep(1,9))),
                               B_L = as.numeric(c(rep(c(rep(-1,3),rep(0,3),rep(1,3)),3))),
                               B_Q = as.numeric(c(rep(c(rep(1,3),rep(-2,3),rep(1,3)),3))),
                               C_L = as.numeric(c(rep(c(-1,0,1),9))),
                               C_Q = as.numeric(c(rep(c(1,-2,1),9)))),
                          .Names = c("treatment","A_L","A_Q","B_L","B_Q","C_L","C_Q"), class = "data.frame",
                          row.names = c(NA, -27L)) # last number in row.names is just nrow()

orthogonal_2$A_L_B_L <- orthogonal_2$A_L*orthogonal_2$B_L
orthogonal_2$A_L_B_Q <- orthogonal_2$A_L*orthogonal_2$B_Q
orthogonal_2$A_Q_B_L <- orthogonal_2$A_Q*orthogonal_2$B_L
orthogonal_2$A_Q_B_Q <- orthogonal_2$A_Q*orthogonal_2$B_Q

orthogonal_2$A_L_C_L <- orthogonal_2$A_L*orthogonal_2$C_L
orthogonal_2$A_L_C_Q <- orthogonal_2$A_L*orthogonal_2$C_Q
orthogonal_2$A_Q_C_L <- orthogonal_2$A_Q*orthogonal_2$C_L
orthogonal_2$A_Q_C_Q <- orthogonal_2$A_Q*orthogonal_2$C_Q

orthogonal_2$B_L_C_L <- orthogonal_2$B_L*orthogonal_2$C_L
orthogonal_2$B_L_C_Q <- orthogonal_2$B_L*orthogonal_2$C_Q
orthogonal_2$B_Q_C_L <- orthogonal_2$B_Q*orthogonal_2$C_L
orthogonal_2$B_Q_C_Q <- orthogonal_2$B_Q*orthogonal_2$C_Q

orthogonal_2$A_L_B_L_C_L <- orthogonal_2$A_L*orthogonal_2$B_L*orthogonal_2$C_L
orthogonal_2$A_L_B_L_C_Q <- orthogonal_2$A_L*orthogonal_2$B_L*orthogonal_2$C_Q
orthogonal_2$A_L_B_Q_C_L <- orthogonal_2$A_L*orthogonal_2$B_Q*orthogonal_2$C_L
orthogonal_2$A_L_B_Q_C_Q <- orthogonal_2$A_L*orthogonal_2$B_Q*orthogonal_2$C_Q
orthogonal_2$A_Q_B_L_C_L <- orthogonal_2$A_Q*orthogonal_2$B_L*orthogonal_2$C_L
orthogonal_2$A_Q_B_L_C_Q <- orthogonal_2$A_Q*orthogonal_2$B_L*orthogonal_2$C_Q
orthogonal_2$A_Q_B_Q_C_L <- orthogonal_2$A_Q*orthogonal_2$B_Q*orthogonal_2$C_L
orthogonal_2$A_Q_B_Q_C_Q <- orthogonal_2$A_Q*orthogonal_2$B_Q*orthogonal_2$C_Q

orthogonal_2$totals  <-  as.numeric(c(desc_concreto3$suma[]))

x0 <- ncol(orthogonal_2)
x1 <- ncol(orthogonal_2)-1
x2 <- ncol(orthogonal_2)-2
sum_c1sq2 <- data.frame(matrix(data = NA,nrow = x2,ncol = 3));colnames(sum_c1sq2)[1:3] <- c('treatment','sumc1sq','latex')
sum_c1sq2$treatment <- colnames(orthogonal_2)[2:x1]
for(i in 1:nrow(sum_c1sq2)){
  sum_c1sq2$sumc1sq[i] <- sum(orthogonal_2[,i+1]^2)
  for(j in 1:nrow(orthogonal_2)){
    if(j == 1){
      tempa <- c(paste0(as.character(colnames(orthogonal_2)[i+1])," Sigma c_i^2= "))
      tempa <- c(paste0(tempa,'(',orthogonal_2[j,i+1],')^2+'))
    }else if(j > 1 & j < nrow(orthogonal_2)){
      tempa <- c(paste0(tempa,'(',orthogonal_2[j,i+1],')^2+'))
    }else if(j == nrow(orthogonal_2)){
      tempa <- c(paste0(tempa,'(',orthogonal_2[j,i+1],')^2'))
      sum_c1sq2$latex[i] <- as.character(tempa)
    }
  }
}
sum_c1sq2$latex <- paste0(sum_c1sq2$latex[],'=',sum_c1sq2$sumc1sq)

contrast2 <- data.frame(matrix(data = NA,nrow = x2,ncol = 3));colnames(contrast2)[1:3] <- c('treatment','contrast','latex')
contrast2$treatment <- colnames(orthogonal_2)[2:x1]
for(i in 2:x1){
  tempo <- orthogonal_2[,i]*orthogonal_2[,x0] 
  contrast2$contrast[i-1] <- sum(tempo)
  for(j in 1:nrow(orthogonal_2)){
    if(j == 1){
      tempa <- c(paste0("C_{",as.character(colnames(orthogonal_2)[i]),"}="))
      tempa <- c(paste0(tempa,'(',orthogonal_2[j,i],')(',orthogonal_2[j,x0],')+'))
    }else if(j > 1 & j < nrow(orthogonal_2)){
      tempa <- c(paste0(tempa,'(',orthogonal_2[j,i],')(',orthogonal_2[j,x0],')+'))
    }else if(j == nrow(orthogonal_2)){
      tempa <- c(paste0(tempa,'(',orthogonal_2[j,i],')(',orthogonal_2[j,x0],')'))
      contrast2$latex[i-1] <- as.character(tempa)
    }
  }
}
contrast2$latex <- paste0(contrast2$latex[],'=',round(contrast2$contrast[],2))

ss2 <- data.frame(matrix(data = NA,nrow = x2,ncol = 3));colnames(ss2)[1:3] <- c('treatment','ss','latex')
ss2$treatment <- colnames(orthogonal_2)[2:x1]
n <- 2 # define n before, number of observations per treatment combination
for(i in 1:nrow(ss2)){
  ss2$ss[i] <- round(contrast2$contrast[i]^2/(n*sum_c1sq2$sumc1sq[i]),3)
  ss2$latex[i] <- c(paste0('SS_{',as.character(colnames(orthogonal_2)[i+1]),'}= frac{(',round(contrast2$contrast[i],3),')^2}{',n,'(',sum_c1sq2$sumc1sq[i],')','}'))
}
ss2$latex <- paste0(ss2$latex[],'=',ss2$ss[])

setwd("/Users/williamcruz/Desktop/3 STATS")
library('openxlsx')
dss <- read.xlsx('test.xlsx',sheet = 1)
dss$pvalue <- c(NA)
for(i in 1:nrow(dss)){
  dss$pvalue[i] <- round(pf(q=dss$F.value[i], df1=dss$df[i], df2=27, lower.tail=FALSE),3)
}

## Problem 10.15
pd <- position_dodge(width = 0.2)
desc0 <- ddply(concrete_db,.(base,subbase),summarise,totals=sum(measure))
desc1 <- ddply(concrete_db,.(surface,base),summarise,totals=sum(measure))
desc2 <- ddply(concrete_db,.(surface,subbase),summarise,totals=sum(measure))
desc3 <- ddply(concrete_db,.(surface),summarise,totals=sum(measure))

(ansPlot1 <-  ggplot(data = desc1,aes(x=base,y=totals,group = surface))+
    geom_line(aes(linetype = surface), position = pd) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Surface")) +
    labs(title = paste("Cell totals of Surface and Base levels"),
         x = "Base",
         y = "Cell totals") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

(ansPlot2 <-  ggplot(data = desc2,aes(x=subbase,y=totals,group = surface))+
    geom_line(aes(linetype = surface), position = pd) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Surface")) +
    labs(title = paste("Cell totals of Surface and Sub-base levels"),
         x = "Sub-base",
         y = "Cell totals") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

(ansPlot0 <-  ggplot(data = desc0,aes(x=subbase,y=totals,group = base))+
    geom_line(aes(linetype = base), position = pd) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Surface")) +
    labs(title = paste("Cell totals of Base and Sub-base levels"),
         x = "Sub-base",
         y = "Cell totals") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

(ansPlot3 <-  ggplot(data = desc3,aes(x=surface,y=totals))+
    geom_line(aes(linetype = surface),position = pd) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    labs(title = paste("Cell totals of Surface levels"),
         x = "Surface",
         y = "Cell totals") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

concrete_db$interaction <- c(NA)
for(i in 1:nrow(concrete_db)){
  if(concrete_db$base[i] == 'Base0' & concrete_db$subbase[i] == 'Sub4'){
    concrete_db$interaction[i]  <- c('00')
  }else if(concrete_db$base[i] == 'Base0' & concrete_db$subbase[i] == 'Sub8'){
    concrete_db$interaction[i] <- c('01')
  }else if(concrete_db$base[i] == 'Base0' & concrete_db$subbase[i] == 'Sub12'){
    concrete_db$interaction[i] <- c('02')
  }else if(concrete_db$base[i] == 'Base3' & concrete_db$subbase[i] == 'Sub4'){
    concrete_db$interaction[i] <- c('10')
  }else if(concrete_db$base[i] == 'Base3' & concrete_db$subbase[i] == 'Sub8'){
    concrete_db$interaction[i] <- c('11')
  }else if(concrete_db$base[i] == 'Base3' & concrete_db$subbase[i] == 'Sub12'){
    concrete_db$interaction[i] <- c('12')
  }else if(concrete_db$base[i] == 'Base6' & concrete_db$subbase[i] == 'Sub4'){
    concrete_db$interaction[i] <- c('20')
  }else if(concrete_db$base[i] == 'Base6' & concrete_db$subbase[i] == 'Sub8'){
    concrete_db$interaction[i] <- c('21')
  }else if(concrete_db$base[i] == 'Base6' & concrete_db$subbase[i] == 'Sub12'){
    concrete_db$interaction[i] <- c('22')
  }
}
concrete_db$interaction <- as.factor(concrete_db$interaction)

desc3 <- ddply(concrete_db,.(surface,interaction),summarise,mean=mean(measure),measure_ci=1.96*sd(measure)/sqrt(length(measure)))
pd <- position_dodge(width = 0.2)

(ansPlot3 <-  ggplot(data = desc3,aes(x=interaction,y=mean,group = surface))+
    geom_line(aes(linetype = surface), position = pd) +
    geom_errorbar(aes(ymin= mean - measure_ci, ymax= mean + measure_ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Transmission")) +
    labs(title = paste("Mean of Concrete resistance depending on",
                       "Surface and Base-Sub base interaction.",
                       "Error bars represent 95% Confidence Intervals",
                       sep = "\n"),
         x = "B x C interaction",
         y = "Mean") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)
######################
# # # Chapter 11 # # #
######################
## Problem 11.7
color_db <- structure(list(replication = structure(c(rep(1L,20),rep(2L,20),rep(3L,20)), .Label = c('I','II','III'), class = "factor"),
                           operator = structure(c(rep(c(seq(1,4)),15)), .Label = LETTERS[seq(from=1,to=4)], class = "factor"),
                           day = structure(c(rep(c(rep(1L,4),rep(2L,4),rep(3L,4),rep(4L,4),rep(5L,4)),3)), .Label = c('Monday','Tuesday','Wednesday','Thrusday','Friday'), class = "factor"),
                           color  = as.numeric(c(800,760,920,860,950,900,920,940,900,820,840,820,740,740,900,940,880,960,1020,980,
                                                 780,820,740,800,810,940,900,900,880,900,880,800,960,780,840,820,920,820,880,900,
                                                 800,900,800,900,1030,920,1000,980,800,880,920,900,840,920,760,780,800,1000,920,880))),
                      .Names = c("replication","operator","day","color"), class = "data.frame",
                      row.names = c(NA, -60L)) # last number in row.names is just nrow()
fit <- EMSanova(color~replication*day*operator,data=color_db,type=c("R","F","F"),nested=c(NA,NA,NA))

## Problem 11.8
fit2 <- EMSanova(color~replication*day*operator,data=color_db,type=c("R","F","R"),nested=c(NA,NA,NA))

## Problem 11.10
rubber_db <- structure(list(replication = structure(c(rep(1L,27),rep(2L,27),rep(3L,27),rep(4L,27)), .Label = c('I','II','III','IV'), class = "factor"),
                            temperature = structure(c(rep(c(rep(1L,9),rep(2L,9),rep(3L,9)),4)), .Label = c('145','155','165'), class = "factor"),
                            lab = structure(c(rep(c(seq(1,3)),36)), .Label = letters[seq(from=1,to=3)], class = "factor"),
                            mix = structure(c(rep(c(rep(1L,3),rep(2L,3),rep(3L,3)),12)), .Label = LETTERS[seq(from=1,to=3)], class = "factor"),
                            cure  = as.numeric(c(18.6,20,19.7,14.5,18.4,16.3,21.1,22.5,22.7,9.5,11.4,9.3,7.8,10.8,9.1,11.2,13.3,11.3,5.4,6.8,6.7,5.2,6,5.7,6.3,7.7,6.6,
                                                 17,20.1,18.3,15.8,18.1,16.7,20.8,22.7,21.9,9.4,11.5,10.2,8.3,11.1,9.2,10,14,11,5.3,6.9,6,4.9,6.1,5.5,6.4,8,6.5,
                                                 18.7,19.4,16.8,16.5,16.5,14.4,21.8,21.5,19.3,9.5,11.4,9.8,8.9,9.5,8,11.5,12,10.9,5.7,6,5,4.3,5,4.6,5.8,6.6,5.9,
                                                 18.7,20,17.1,17.6,16.7,15.2,21,21.3,19.3,10,11.5,9.5,9.1,9.7,9,11.1,11.5,11.4,5.3,5.7,4.8,5.2,5.2,5.4,5.6,6.3,5.8))),
                       .Names = c("replication","temperature","lab","mix","cure"), class = "data.frame",
                       row.names = c(NA, -108L)) # last number in row.names is just nrow()
fit_rubber <- EMSanova(cure~replication*lab*temperature*mix,data=rubber_db,type=c("R","F","F","F"),nested=c(NA,NA,NA,NA))

desc_rubber <- ddply(rubber_db,.(temperature,mix),summarise,mean=mean(cure),ci=1.96*sd(cure)/sqrt(length(cure)))
(ansPlot <-  ggplot(data = desc_rubber,aes(x=temperature,y=mean,group = mix))+
    geom_line(aes(linetype = mix), position = pd) +
    geom_errorbar(aes(ymin= mean - ci, ymax= mean + ci),
                  width = .2,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Mix")) +
    labs(title = paste("Mean cure rate variation by",
                       "temperature and mix.",
                       "Error bars represent 95% Confidence Intervals",
                       sep = "\n"),
         x = "Temperature",
         y = "Cure mean rate") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

## Problem 11.11
(ansplot2 <- ggplot(db, aes(x=TC, y=cure, fill=TC)) +
    geom_boxplot() +
    labs(title = paste("Treatment Combinations between",
                       "temperature and mix.",
                       sep = "\n"),
         x = "Treatment Combinations",
         y = "Cure rate") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

##Problem 11.12
fit_rubber2 <- EMSanova(cure~replication*lab*temperature*mix,data=rubber_db,type=c("R","R","F","F"),nested=c(NA,NA,NA,NA))

##Problem 11.16
cake_db <- structure(list(replication = structure(c(rep(1L,18),rep(2L,18)), .Label = c('I','II'), class = "factor"),
                          temperature = structure(c(rep(c(rep(1L,3),rep(2L,3),rep(3L,3),rep(4L,3),rep(5L,3),rep(6L,3)),2)), .Label = c('175','185','195','205','215','225'), class = "factor"),
                          mix = structure(c(rep(c(seq(1,3,1)),12)), .Label = letters[seq(from=1,to=3)], class = "factor"),
                          quality  = as.numeric(c(47,35,43,29,46,43,35,47,43,47,39,46,57,52,47,45,61,58,
                                                  26,21,21,28,21,28,32,28,25,25,26,25,37,27,31,33,20,25))),
                     .Names = c("replication","temperature","mix","quality"), class = "data.frame",
                     row.names = c(NA, -36L)) # last number in row.names is just nrow()

fit_cake <- EMSanova(quality~replication*mix*temperature,data=cake_db,type=c("R","F","F"),nested=c(NA,NA,NA))

######################
# # # Chapter 15 # # #
######################
## Example 15.1
dosage_db <- structure(list(dosage = as.numeric(c(rep(0.5,3),rep(1,3),rep(1.5,3),rep(2,3),rep(2.5,3))),
                            time  = as.numeric(c(26,28,29,28,26,30,28,30,31,32,33,31,38,39,38))),
                       .Names = c("dosage","time"), class = "data.frame",
                       row.names = c(NA, -15L)) # last number in row.names is just nrow()

dosage_db$dosage <- as.factor(dosage_db$dosage) # there might be useful to run it without being a factor an converting it into factor
# to estimate the departure from linear SS and MS. Just substract the residual SS from the first analysis
# to that of the second analysis.
dosage_db$dosage2 <- dosage_db$dosage^2
plot(x=dosage_db$dosage,y=dosage_db$time)
scatter.smooth(x=dosage_db$dosage, y=dosage_db$time, main="time ~ dosage")

linearmodel <- lm(time ~ dosage,data = dosage_db)
linearmodel2 <- lm(time ~ dosage+dosage2,data = dosage_db)
summary(linearmodel)
summary(linearmodel2)
dosage_db$dosage2 <- as.factor(dosage_db$dosage2)
summary(aov(time~dosage+dosage2,data = dosage_db))

## Problem 15.5
planting_db <- structure(list(planting = structure(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3)),.Label = c('12','16','20','24','28'), class = "factor"),
                              yield  = as.numeric(c(130.5,129.6,129.9,
                                                    142.5,140.3,143.4,
                                                    145.1,144.8,144.1,
                                                    147.8,146.6,148.4,
                                                    134.8,135.1,136.7))),
                         .Names = c("planting","yield"), class = "data.frame",
                         row.names = c(NA, -15L)) # last number in row.names is just nrow()

#version2
planting2_db <- structure(list(planting = as.numeric(c(rep(12,3),rep(16,3),rep(20,3),rep(24,3),rep(28,3))),
                               yield  = as.numeric(c(130.5,129.6,129.9,
                                                     142.5,140.3,143.4,
                                                     145.1,144.8,144.1,
                                                     147.8,146.6,148.4,
                                                     134.8,135.1,136.7))),
                          .Names = c("planting","yield"), class = "data.frame",
                          row.names = c(NA, -15L)) # last number in row.names is just nrow()

ans_plot <- scatter.smooth(x=as.factor(planting_db$planting), y=planting_db$yield, main="Planting rate and Yield scatter plot",
                           xlab = 'Planting rate', ylab = 'Yield')
scale_y_continuous(limits = c(11,29),breaks = seq(12,28,4)) + 
  theme(
    panel.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white"),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1)
  )

ggplot(planting2_db, aes(x=planting, y=yield)) + 
  geom_point()+
  geom_smooth(method = lm,formula = y ~ splines::bs(x, 3),se = FALSE)
geom_smooth(method = lm,se = FALSE)+
  labs(title = "Scatter plot of Planting rate vs Yield",
       x = "Planting rate",
       y = "Yield") +
  scale_x_continuous(limits = c(11,29),breaks = seq(12,28,4)) + 
  theme(
    panel.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white"),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1)
  )
linearmodel3 <- lm(yield ~ planting,data = planting2_db)
##Problem 15.6
planting2_db$planting2 <- planting2_db$planting^2
planting2_db$planting3 <- planting2_db$planting^3
planting2_db$planting4 <- planting2_db$planting^4
linearmodel3 <- lm(yield ~ planting + planting2 + planting3 + planting4,data = planting2_db)
summary(linearmodel3)
A <- matrix(c(15, 0, 30, # unsuccessful trial, unable to install the matlib package
              0, 30, 0,
              30,  0, 102), 3, 3, byrow=TRUE)
b <- c(2099.6, 49.8, 4055.4)
gaussianElimination(A, b)
gaussianElimination(A, b, verbose=TRUE, fractions=TRUE)
gaussianElimination(A, b, verbose=TRUE, fractions=TRUE, latex=TRUE)

pf(q=4.59, df1=1, df2=10, lower.tail=FALSE)

##Problem 15.7
planting2_db$planting3 <- planting2_db$planting^3
planting2_db$planting4 <- planting2_db$planting^4
linearmodel3 <- lm(yield ~ planting + planting2 + planting3 + planting4,data = planting2_db)
summary(linearmodel3)

##Problem 15.8
#start from plain planting2_db
planting2_db$planting2 <- planting2_db$planting^2 
planting2_db$planting2 <- as.factor(planting2_db$planting2)
model0 <- summary(aov(yield~planting+planting2,data = planting2_db))
linearmodel3 <- lm(yield ~ planting + planting2,data = planting2_db)
summary(linearmodel3)

##Problem 15.18
film_db <- structure(list(resin  = structure(c(rep(1,18),rep(2,18)),.Label = c('1','2'), class = "factor"),
                          weight = structure(c(rep(c(rep(1L,6),rep(2L,6),rep(3L,6)),2)), .Label = c('0.20','0.25','0.30'), class = "factor"),
                          gate   = structure(c(rep(c(rep(1L,2),rep(2L,2),rep(3L,2)),6)), .Label = c('2','4','6'), class = "factor"),
                          thickness  = as.numeric(c(1.6,1.5,2.7,2.7,3.9,4.0,
                                                    1.5,1.3,2.5,2.5,3.6,3.8,
                                                    1.5,1.3,2.4,2.3,3.5,3.4,
                                                    1.5,1.4,2.4,2.3,4.0,4.0,
                                                    1.4,1.3,2.6,2.4,3.7,3.6,
                                                    1.6,1.4,2.2,2.1,3.4,3.3))),
                     .Names = c("resin","weight","gate","thickness"), class = "data.frame",
                     row.names = c(NA, -36L)) # last number in row.names is just nrow()
fit_film <- EMSanova(thickness~resin*gate*weight,data=film_db,type=c("F","F","F"),nested=c(NA,NA,NA))

###Problem 15.18
pf(q=1818.100, df1=2, df2=18, lower.tail=FALSE)

##Problem 15.21
pd <- position_dodge(width = 0.2)
desc0 <- ddply(film_db,.(gate),summarise,totals=sum(thickness))
desc1 <- ddply(film_db,.(resin),summarise,mean=mean(thickness),ci=1.96*sd(thickness)/sqrt(length(thickness)))
desc2 <- ddply(film_db,.(gate,weight),summarise,totals=sum(thickness))
desc3 <- ddply(film_db,.(gate,weight),summarise,mean=mean(thickness),ci=1.96*sd(thickness)/sqrt(length(thickness)))


(ansPlot1 <-  ggplot(data = desc0,aes(x=gate,y=totals))+
    geom_smooth(method = lm,formula = y ~ splines::bs(x, 2),se = FALSE)+
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "black") +
    labs(title = paste("       Cell totals of thickness by Weight"),
         x = "Gate",
         y = "Thickness") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

(ansPlot2 <-  ggplot(data = desc1,aes(x=resin,y=mean))+
    geom_errorbar(aes(ymin= mean - ci, ymax= mean + ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    labs(title = paste("   Resin Types means",
                       "Error bars represent 95% Confidence Intervals",
                       sep = "\n"),
         x = "Resin Types",
         y = "Thickness") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

(ansPlot3 <-  ggplot(data = desc2,aes(x=gate,y=totals,group = weight))+
    geom_line(aes(linetype = weight), position = pd) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Weight")) +
    labs(title = paste("Cell totals of Gate and Weight levels"),
         x = "Gate",
         y = "Cell totals") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

(ansPlot <-  ggplot(data = desc3,aes(x=gate,y=mean,group = weight))+
    geom_line(aes(linetype = weight), position = pd) +
    geom_errorbar(aes(ymin= mean - ci, ymax= mean + ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Weight")) +
    labs(title = paste("         Mean thickness depending on",
                       "             Gate and Weight",
                       "Error bars represent 95% Confidence Intervals",
                       sep = "\n"),
         x = "Gate",
         y = "Thickness") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

##Problem 15.34
pf(q=12.7514, df1=1, df2=46, lower.tail=FALSE)

##Problem 15.36
model_db <- structure(list(model=structure(c(1L,2L,3L,4L),.Label = c('model1','model2','model3','model4'), class="factor"),
                           coefficient=as.numeric(c(0.0907,0.3285,0.4211,0.4548))),
                      .Names = c("model","coefficient"), class = "data.frame",
                      row.names = c(NA, -4L))

(ansPlot1 <-  ggplot(data = model_db,aes(x=model,y=coefficient))+
    geom_smooth(method = lm,formula = y ~ splines::bs(x, 3),se = FALSE)+
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    labs(title = paste("       Models coefficients"),
         x = "model",
         y = "R^2") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)
x <- ((0.4548-0.4211)/(4-3))/((1-0.4548)/(45-4-1))
pf(q=x, df1=1, df2=40, lower.tail=FALSE)

#################
# # # Tests # # #
#################
####test for lines linking factors

geom_smooth(method = lm,formula = y ~ splines::bs(x, 3),se = FALSE)
geom_smooth(method = lm,se = FALSE)+
  pd <- position_dodge(width = 0.2)
mtcars$cyl <- as.factor(mtcars$cyl)
mtcars$am <- factor(mtcars$am,labels = c('automatic','manual'))

hp <- ddply(mtcars,.(cyl,am),summarise,hp_mean=mean(hp),hp_ci=1.96*sd(hp)/sqrt(length(hp)))
(ansPlot <-  ggplot(data = hp,aes(x=cyl,y=hp_mean,group = am))+
    geom_line(aes(linetype = am), position = pd) +
    geom_errorbar(aes(ymin= hp_mean - hp_ci, ymax= hp_mean + hp_ci),
                  width = .1,position = pd, linetype = 1) +
    geom_point(size = 4, position = pd) +
    geom_point(size = 3, position = pd, color = "white") +
    guides(linetype = guide_legend("Transmission")) +
    labs(title = paste("Mean horsepower depending on",
                       "number of cylinders and transmission type.",
                       "Error bars represent 95% Confidence Intervals",
                       sep = "\n"),
         x = "Number of cylinders",
         y = "Gross horsepower") +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black", size = 1),
      axis.line.y = element_line(colour = "black", size = 1)
    )
)

### end of thes figure linking factors
