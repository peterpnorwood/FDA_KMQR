library(ggplot2)

setwd("~/KMQR")
source("gen_power.R")

## Running for 100 curves, kernel = 1, X2 Method

## Running the function for each delta level
## Change the parameters n, k, x2 for the exact situation you want

## Note: you can also use mclapply here
x2_1_100 <- lapply(delta_vec,complete_process, n=100,k=1,x2=T)

## Formatting the Results
results_x2_1_100 <- matrix(unlist(x2_1_100),ncol = 2,nrow=11,byrow = T)
results_x2_1_100_format <- data.frame(as.numeric(results_x2_1_100),c(seq(0,2, len = 11),seq(0,2, len = 11)),c(rep("QR",11),rep("KMQR",11)))
colnames(results_x2_1_100_format) <- c("pvalue","delta","test")

## Plotting the Results
plot_x2_1_100 <- ggplot(results_x2_1_100_format,aes(delta,pvalue)) + geom_point(aes(color=as.factor(test))) +
  geom_line(aes(group = as.factor(test),color=as.factor(test))) + xlab("Delta") + 
  ylab("Power") + ggtitle("X2, Power Raised = 1, 100 Observations") + 
  scale_colour_discrete(name ="Test")