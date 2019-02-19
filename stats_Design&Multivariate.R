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

# ANOVA for the 'Power Requirement Experiment'
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