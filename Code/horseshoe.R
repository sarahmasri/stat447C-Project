library(readr)
library(ramify)
library(glmnet)
library(tidyverse)
library(plotrix)
library(DescTools)
library(rstan)


setwd("/Users/sarahmasri/Desktop/STAT  447C/stat447C-Project/Code/")


Normalized_TrainData <- read.table("../Data/Normalized_TrainData.txt")
Normalized_TestData <- read.table("../Data/Normalized_TestData.txt")
Normalized_TrainData_Th <- read.table("../Data/Normalized_TrainData_Th.txt")
Normalized_TestData_Th <- read.table("../Data/Normalized_TestData_Th.txt")

Normalized_BC <- read.table("../Data/Normalized_BC.txt")



####Classifier for Immune Cells---------------------------------------------
##Class Labels: 1:B Cell, 2:CD4, 3:CD8, 4:Mono, 
##5:Neu, 6:NK, 7:DC, 8:M1, 9:M2, 10:MQ




x <-log2(as.matrix(data.frame(t(Normalized_TrainData[,3:110])))+1)
Test <-log2(as.matrix(data.frame(t(Normalized_TestData[,c(3:62,65:81)])))+1) 

y<-rbind(ones(10, ncol = 1) , 2*ones(10, ncol = 1), 
         3*ones(10, ncol = 1),4*ones(10, ncol = 1),
         5*ones(10, ncol = 1),6*ones(10, ncol = 1),
         1,7*ones(6, ncol = 1),4 ,5 ,6 ,8 ,8, 9, 9,
         2 ,2 ,3 ,3,3,10,10,10,6,ones(3, ncol = 1),
         2*ones(3, ncol = 1) , 3*ones(3, ncol = 1),
         4*ones(3, ncol = 1) , 6*ones(3, ncol = 1),
         7,7,4,4,8,8,9,9,10,10)


y_Test<-rbind( ones(10, ncol = 1), 2*ones(10, ncol = 1),
               3*ones(10, ncol = 1),4*ones(10, ncol = 1),
               5*ones(10, ncol = 1), 6*ones(4, ncol = 1),
               1,7,7,4,5,6,8,9,2,3,3,10,6,1,2,3,4,6,7, 4,
               8,9,10)


##Change Missing Values to -1

x[is.na(x)]<- -1
Test[is.na(Test)]<- -1






## Bayesian regression with Horshoe prior

# Define data
n <- nrow(x)       # Number of observations
p <- ncol(x)       # Number of predictors 
X <- as.matrix(x)  # Predictor matrix
y <- as.vector(y)  # Outcome vector


# Compile Stan model
file <- "./horseshoe.stan"  # Paste the Stan model code here
model <- stan_model(file = file)


# Fit the model
fit <- sampling(model, 
                data = list(n = n, p = p, X = X, y = y), 
                iter = 2000)

coefs <- data.frame(summary(fit, pars = c("beta"))$summary)$mean
horseshoe <- data.frame(coefs)


#write.table(horseshoe, "../Data/horseshoe")