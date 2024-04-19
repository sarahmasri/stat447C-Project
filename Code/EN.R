

library(readr)
library(ramify)
library(glmnet)
library(tidyverse)
library(plotrix)
library(DescTools)


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


##Create Noisy Data to Neglect -1 and Treat It as Missing Values 

#X<-x ; Y<-y
#for (i in 1:100) {
#  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
#  Y<-rbind(Y,y)
#}








##Generation of Classifier (Elastic-Net Logistic Regression)

lambdas <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)


start_time <- Sys.time()
EN <- cv.glmnet(x, y, family="multinomial", alpha=.93, lambda=lambdas, type.multinomial="grouped", nfolds=5)
end_time <- Sys.time()
elapsed_time <- end_time - start_time

(elapsed_time)
##40.89618 secs













##Accuracy of Classifier for Testing and Training Samples

prediction <- data.frame(predict(EN, newx=data.matrix(x), 
                                 interval ="prediction",
                                 type="response", s=0.0001))

colnames(prediction)<-c("B Cell","CD4","CD8","Mono","Neu",
                        "NK","DC","M1","M2","MQ")

which(as.matrix(prediction)==matrixStats::rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(EN, newx=data.matrix(Test), 
                                      interval ="prediction",
                                      type="class", s=0.0001))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy












####Scrambling Data to Creat True Negative Samples and Plot ROC-------------




TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-0*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0


for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  prediction_Test <- data.frame(predict(EN, newx=data.matrix(Test),
                                        interval ="prediction",
                                        type="response", s=i))
  
  prediction_Test_Class <- data.frame(predict(EN, newx=data.matrix(Test),
                                              interval ="prediction",
                                              type="class", s=i))
  
  prediction_TN_Test <- data.frame(predict(EN, newx=data.matrix(TN_Test),
                                           interval ="prediction",
                                           type="response", s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0; FN=0; TN=0; FP=0; R=R+1
    
    for (j in 1:length(TN_y)) {
      
      if (matrixStats::rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
          
        } else {FN=FN+1}
        
      } else {FN=FN+1}
      
      if (matrixStats::rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
        
      } else {FP=FP+1}
    }
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


##Plot ROC Curve



pdf( "../Figures/ROC_New.pdf", width = 5, height =5 )

Color=c("red","orange","yellow","forestgreen","blue","black","violet")
matplot(1-Specificity,Sensitivity, type="l", col=Color,lwd = 1.5, 
        lty=c(1,2,1,4,5,6,2),xlab = "1-Specificity", 
        ylab = "Sensitivity", ylim = c(.5,1))

Legend=c("lambda=1e-1","lambda=5e-2","lambda=1e-2","lambda=5e-3",
         "lambda=1e-3 ","lambda=5e-4","lambda=1e-4")
legend("bottomright", legend = Legend , col = Color,
       lty = c(1,2,1,4,5,6,2), lwd = 1.5, cex = 0.6)

dev.off()


















###################################################
###################################################
################ BREAST CANCER ###################
###################################################
###################################################

####Breast Cancer-------------
#Data was obtained from: (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688)



##Normalization of Breast Cancer Samples------




##Analysis over Breast Cancer Data------


y_Test=rbind(ones(83,ncol=1),2*ones(54,ncol=1),
             11*ones(38,ncol=1),12*ones(3,ncol=1))

Test <-log2(as.matrix(data.frame(t(Normalized_BC[,-c(1,2)])))+1)

Test[is.na(Test)]<- -1

prediction <- data.frame(predict(EN, newx=data.matrix(Test),
                                 interval ="prediction",
                                 type="response", s=0.0001))

colnames(prediction)<-c("B Cell","CD4","CD8","Mono","Neu",
                        "NK","DC","M1","M2","MQ")

Class<-which(as.matrix(prediction)==matrixStats::rowMaxs(as.matrix(prediction)),
             arr.ind = TRUE)

Class<-Class[order(Class[,1]),]

prediction_Class <-predict(EN, newx=data.matrix(Test),
                           interval ="prediction",type="class",
                           s=0.0001)


Correct=0
for (i in min(which(y_Test %in% 1)):max(which(y_Test %in% 1))){
  if(prediction_Class[i,1]==1){
    Correct=Correct+1
  }
}

for (i in min(which(y_Test %in% 2)):max(which(y_Test %in% 2))){
  if(prediction_Class[i,1]==2|prediction_Class[i,1]==3){
    Correct=Correct+1
  }
}

Total_Accuracy<-Correct/max(which(y_Test %in% 2))
Total_Accuracy # 0.3649635

TCell=0
BCell=0
M=0
B=0
j=0
k=0

FN_TCell=c()
FN_BCell=c()
FN_M=c()
FN_B=c()

TP_TCell=c()
TP_BCell=c()
TP_M=c()
TP_B=c()

for (i in min(which(y_Test %in% 2)):max(which(y_Test %in% 2))){
  if(prediction_Class[i,1]==2|prediction_Class[i,1]==3){
    TCell=TCell+1
    k<-k+1
    TP_TCell[k]<-prediction_Class[i,1]
  }else{j<-j+1;FN_TCell[j]<-prediction_Class[i,1]}
}


j=0
k=0

for (i in min(which(y_Test %in% 1)):max(which(y_Test %in% 1))){
  if(prediction_Class[i,1]==1){
    BCell=BCell+1
    k<-k+1
    TP_BCell[k]<-prediction_Class[i,1]
  }else{j<-j+1;FN_BCell[j]<-prediction_Class[i,1]}
}

j=0
k=0

for (i in min(which(y_Test %in% 11)):max(which(y_Test %in% 11))){
  if(prediction_Class[i,1]==0){
    M=M+1
    k<-k+1
    TP_M[k]<-prediction_Class[i,1]
  }else{j<-j+1;FN_M[j]<-prediction_Class[i,1]}
}

j=0
k=0

for (i in min(which(y_Test %in% 12)):max(which(y_Test %in% 12))){
  if(prediction_Class[i,1]==0){
    B=B+1
    k<-k+1
    TP_B[k]<-prediction_Class[i,1]
  }else{j<-j+1;FN_B[j]<-prediction_Class[i,1]}
}


TCell/length(which(y_Test%in%2)) #0.4814815
BCell/length(which(y_Test%in%1)) #0.2891566
M/length(which(y_Test%in%11))    #0
B/length(which(y_Test%in%12))    #0












