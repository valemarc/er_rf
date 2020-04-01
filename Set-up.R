###Load devtools
library(devtools)

####Check status in Git
use_git()
use_github()


#install packages
install.packages("caret", dependencies=TRUE)
install.packages("ROCR")
install.packages("kernlab")
install.packages("e1071")
install.packages("rpart")
install.packages("caTools")
install.packages("spdep")
install.packages("ncf")
install.packages("gbm")
install.packages("pROC")
install.packages("randomForest")
install.packages("rpart.plot")
install.packages("ellipse")

install.packages("spdep")
install.packages("here")
install.packages("corrplot")
install.packages("pROC")
install.packages("ranger")
devtools::install_cran("dplyr", force=TRUE)
install.packages("MLeval")

#Open libraries
library(caret)
library(ROCR)
library(kernlab)
library (e1071)
library(rpart)
library(caTools)
library(spdep)
library(ncf)
library(gbm)
library(pROC)
library(randomForest)
library(rpart.plot)
library(ellipse)

library (spdep)
library(here)
library(corrplot)
library(pROC)
library(ranger)
library(dplyr)
library(MLeval)
