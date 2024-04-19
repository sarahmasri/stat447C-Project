library(rstantools)

setwd("/Users/sarahmasri/Desktop/STAT  447C/stat447C-Project/Code/")
horseshoe <- read.table("../Data/horseshoe.txt")

## Compare

#all_genes <- Normalized_TrainData[,2]
#horseshoe$genes <- all_genes



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
EN_Selected_Genes<-Normalized_TrainData[Index, 2]



## Horseshoe

threshold <- 0.1
Index <- which(abs(horseshoe$coefs) >= threshold)
horseshoe_Selected_Genes<-Normalized_TrainData[Index, 2]


## Compare
length(EN_Selected_Genes) ## 110
length(horseshoe_Selected_Genes) ## 1946
length(Normalized_TrainData[, 2]) ## 15453

sum(EN_Selected_Genes %in% horseshoe_Selected_Genes) ## 0




