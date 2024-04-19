library(rstantools)



## Compare

all_genes <- Normalized_TrainData[,2]
horseshoe$genes <- all_genes



## Combine coefficients
n <- length(EN_coefs)
est_coefs <- data.frame(matrix(nrow=n, ncol=0))
est_coefs$genes <- all_genes
est_coefs$EN <- EN_coefs

Index <-coef(EN,s=0.0001)[["1"]]@i[-1]
EN_genes <- all_genes[Index]





est_coefs$horseshoe <- horseshoe$coefs

threshold <- 0.5
Index <- which(abs(horseshoe$coefs) <= threshold)
est_coefs$horseshoe[Index] <- 0

Index <- which(abs(horseshoe$coefs) >= threshold)
horseshoe_genes <- all_genes[Index]



est_coefs$diff <- est_coefs$EN - est_coefs$horseshoe

( mean(est_coefs$diff) )



## Find selected genes

## EN
EN_coefs <- coef(EN,s=0.0001)[["1"]][-1]

Index <-coef(EN,s=0.0001)[["1"]]@i[-1]
EN_Selected_Genes<-Normalized_TrainData[Index, ]
row.names(EN_Selected_Genes)<-EN_Selected_Genes[,2]



## Horseshoe

threshold <- 0.5
Index <- which(abs(horseshoe$coefs) >= threshold)
horseshoe_Selected_Genes<-Normalized_TrainData[Index, ]
horseshoe_Selected_Genes[,2][which(is.na(horseshoe_Selected_Genes[,2]))] <- "NA"
row.names(horseshoe_Selected_Genes)<-horseshoe_Selected_Genes[,2]

horseshoe_Selected_Genes








########################################
########################################

##Analysis over Breast Cancer Data------


## EN

y_Test=rbind(ones(83,ncol=1),2*ones(54,ncol=1),
             11*ones(38,ncol=1),12*ones(3,ncol=1))

Normalized_BC <- read.table("../Data/Normalized_BC.txt")

Test <-log2(as.matrix(data.frame(t(Normalized_BC[,-c(1,2)])))+1)

Test[is.na(Test)]<- -1

prediction <- data.frame(predict(EN, newx=data.matrix(Test),
                                 interval ="prediction",
                                 type="response", s=0.0001))

colnames(prediction)<-c("B Cell","CD4","CD8","Mono","Neu",
                        "NK","DC","M1","M2","MQ")

Class<-which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)),
             arr.ind = TRUE)

Class<-Class[order(Class[,1]),]

prediction_Class <-predict(Lasso, newx=data.matrix(Test),
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
Total_Accuracy #0.3649635

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











## Horseshoe

y_Test=rbind(ones(83,ncol=1),2*ones(54,ncol=1),
             11*ones(38,ncol=1),12*ones(3,ncol=1))

Normalized_BC <- read.table("../Data/Normalized_BC.txt")

Test <-log2(as.matrix(data.frame(t(Normalized_BC[,-c(1,2)])))+1)

Test[is.na(Test)]<- -1



posterior_pred <- posterior_predict(fit, Test)


prediction <- data.frame(posterior_predict(fit, newx=data.matrix(Test)))

colnames(prediction)<-c("B Cell","CD4","CD8","Mono","Neu",
                        "NK","DC","M1","M2","MQ")

Class<-which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)),
             arr.ind = TRUE)

Class<-Class[order(Class[,1]),]

prediction_Class <-predict(Lasso, newx=data.matrix(Test),
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
Total_Accuracy #0.3649635

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





