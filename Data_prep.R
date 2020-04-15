##################################################################################################
### RANDOM FOREST MODEL FOR PREDICTING EXTINCTION RISK IN DATA DEFICIENT REPTILES ################
##################################################################################################

#Open data
here()
#setwd("insertpath") # insert correct path name
BM<-read.csv(here("Bland_data.csv")) # insert correct path name
dim(BM) # check dataset
names(BM)
summary(BM)
typeof(BM)
str(BM)

####################################################
###Create training and testing datasets#############
####################################################

###Partition the dataset
# Separate non DD species
assessed<- data.frame(subset(BM, RedList != "DD"))
#row.names(assessed)<- assessed$Binomial
nrow (assessed)

###Add binary variable threatened/not threatened
assessed$binary <-"nonThr"
assessed$binary[assessed$RedList == "VU"]="Thr"
assessed$binary[assessed$RedList == "EN"]="Thr"
assessed$binary[assessed$RedList == "CR"]="Thr"

# Setting DD species aside
nonassessedspecies<- data.frame(subset (BM, RedList=="DD"))
#row.names(nonassessedspecies)<- nonassessedspecies$Binomial
#nonassessedspecies<- subset(nonassessedspecies,select= -c(Binomial,RedList, Genus))

# No partitioning # as we separate DD species out rather than setting aside a certain proportion of the dataset chosen at random
trainassessed <- assessed
trainassessed <- trainassessed[complete.cases(trainassessed), ]
#trainB <- trainassessed[,36]
trainB <- trainassessed[,"binary"]  # define classes
#trainB <- as.numeric(trainB)

###Check that the length of TrainB matches the number of rows in trainassessed
dim(trainassessed)
length(trainB)

###If you need to check anything you can write out the files
#write.csv(trainassessed, "trainassessed.csv")
#write.csv(trainB, "trainB.csv")


###Create a smaller sample for testing
trainassessed_small <- sample_n(trainassessed, 200)
summary(trainassessed_small)

#######################################################################################################################
#########################  Basic model to test  #################################
###Run model (calculating how long it takes to run)
start_time <- Sys.time()

model <- train(binary ~ ., data=trainassessed,
               method = "ranger")

end_time <- Sys.time()
end_time - start_time

###Look at results
print(model)
plot(model)

###Predict status of DD species
er_pred <- predict(model, nonassessedspecies)
summary(er_pred)

###########################################################################################################
####Repeat model tuning parameters (takes about 1 min)
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10)
# run a random forest model
start_time <- Sys.time()
model <- train(as.factor(binary) ~ ., 
               data = trainassessed, 
               method = "ranger",
               trControl = model_control)
end_time <- Sys.time()
end_time - start_time

###Look at results
print(model)
plot(model)

###Predict status of DD species
er_pred <- predict(model, nonassessedspecies)
summary(er_pred)


######################################################
####compute ROC and AUC under ROC after training#######
#########################################################
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = TRUE, 
  classProbs = TRUE, )
# run a random forest model
start_time <- Sys.time()
model <- train(as.factor(binary) ~ ., 
               data = trainassessed, 
               method = "ranger",
               trControl = model_control)
end_time <- Sys.time()
end_time - start_time

###Look at results
print(model)
plot(model)

###Predict status of DD species
er_pred <- predict(model, nonassessedspecies)
summary(er_pred)


######################################################
####Same but with method rf#######
#########################################################
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = TRUE, 
  classProbs = TRUE, )
# run a random forest model
start_time <- Sys.time()
model <- train(as.factor(binary) ~ ., 
               data = trainassessed, 
               method = "rf",
               trControl = model_control)
end_time <- Sys.time()
end_time - start_time

###Look at results
print(model)
plot(model)

x <- evalm(model)

## get roc curve plotted in ggplot2

x$roc

## get AUC and other metrics

x$stdres

###Predict status of DD species
er_pred <- predict(model, nonassessedspecies)
summary(er_pred)


#######################################################################################
####set ntree to 500 (although this should be the default value set by rf in caret)####
#######################################################################################
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = TRUE, 
  classProbs = TRUE)
# run a random forest model
start_time <- Sys.time()
model <- train(as.factor(binary) ~ ., 
               data = trainassessed, 
               method = "rf",
               ntree= 500,
               trControl = model_control)
end_time <- Sys.time()
end_time - start_time

print(model)
plot(model)


x <- evalm(model)

## get roc curve plotted in ggplot2

x$roc

## get AUC and other metrics

x$stdres

###Predict status of DD species
er_pred <- predict(model, nonassessedspecies)
summary(er_pred)


#####################################################################################
####Set mtry (number of variables randomly sampled as candidates at each split)######
#####################################################################################
###Takes about 26 minutes
#Generate 15 random values of mtry at each time tunning. 
#We have 15 values because of tunning length is 15.
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = TRUE, 
  classProbs = TRUE)
mtry <- sqrt(ncol(trainassessed[,2:29]))
# run a random forest model
start_time <- Sys.time()
model <- train(as.factor(binary) ~ ., 
               data = trainassessed, 
               method = "rf",
               tunelength=15,
               ntree= 500,
               trControl = model_control)
end_time <- Sys.time()
end_time - start_time

print(model)
plot(model)

###Predict status of DD species
er_pred <- predict(model, nonassessedspecies)
summary(er_pred)

##############################################################################
#mtry: Number of variables randomly sampled as candidates at each split.
###Takes about 9 minutes
####Prediction changes from all threatened to all non threatened (mtry was held constant at 5.29 and
#plot(model) throws an error message - There are no tuning parameters with more than 1 value)
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = TRUE, 
  classProbs = TRUE)
mtry <- sqrt(ncol(trainassessed[,2:29]))
tunegrid <- expand.grid(.mtry=mtry)
# run a random forest model
start_time <- Sys.time()
model <- train(as.factor(binary) ~ ., 
               data = trainassessed, 
               method = "rf",
               tuneGrid=tunegrid,
               ntree= 500,
               trControl = model_control)
end_time <- Sys.time()
end_time - start_time

print(model)
plot(model)


x <- evalm(model)

## get roc curve plotted in ggplot2

x$roc

## get AUC and other metrics

x$stdres

###Predict status of DD species
er_pred <- predict(model, nonassessedspecies)
summary(er_pred)

###############################################################
###Pre-process variables as per Caret package (Kuhn, 2008)#####
###############################################################
###Centre and scale variables
trainassessed_preprocess <- preProcess(select(trainassessed, - c("ID", "binary")), method = c("center", "scale", "nzv"))
trainassessed_preprocess
trainassessed_preprocess$method

###As the function preprocess doesn't actually preprocess the data, we need to do this
trainTransformed <- predict(trainassessed_preprocess, trainassessed)

###########################################################################################################
####Repeat model tuning parameters (takes about 22 min)
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = TRUE, 
  classProbs = TRUE)
mtry <- sqrt(ncol(trainTransformed[,2:29]))
# run a random forest model
start_time <- Sys.time()
model <- train(as.factor(binary) ~ ., 
               data = trainTransformed, 
               method = "rf",
               tunelength=15,
               ntree= 500,
               trControl = model_control)
end_time <- Sys.time()
end_time - start_time

print(model)
plot(model)

###Predict status of DD species
er_pred <- predict(model, nonassessedspecies)
summary(er_pred)

x <- evalm(model)

## get roc curve plotted in ggplot2

x$roc

## get AUC and other metrics

x$stdres

###Predict status of DD species
er_pred <- predict(model, nonassessedspecies)
summary(er_pred)
###All species predicted to be threatened

#####################################################################################################
############FINAL MODEL##############################################################################
#####################################################################################################
###Centering and scaling numeric variables################
###Centre and scale variables
trainassessed_preprocess <- preProcess(BM[,-1], method = c("center", "scale"))
trainassessed_preprocess
trainassessed_preprocess$method

###As the function preprocess doesn't actually preprocess the data, we need to do this
trainTransformed <- predict(trainassessed_preprocess, BM)
head(trainTransformed)
#preProcValues 
#BM_new2_Transformed <- predict(preProcValues, BM_new2)
#testTransformed <- predict(preProcValues, test)

#dim(BM_new2_Transformed)
#names(BM_new2_Transformed)
#summary(BM_new2_Transformed)
#typeof(BM_new2_Transformed)
#str(BM_new2_Transformed)


###Create dataframe with different types of variables
BM_other <- trainTransformed[,c("ID", "Binomial", "Order", "Family", "Genus", "BodySize", "HabitatsIUCN", "EOO", "Latitude", "ElevMin", "Prec", "PrecSeas", "Temp", "TempSeas", "HPD", "HPDMin", "HumanFootprint", "Accessibility", "Afrotropical", "Australasia", "Indo_malayan", "Nearctic", "Neotropical", "Oceania", "Palearctic")]
BM_cat <- trainTransformed[,c("ID", "RedList","ReproductiveMode", "TrophicGroup", "HabitatMode", "Continent")]
BM_numeric <- trainTransformed[,c("ID", "BodySize", "HabitatsIUCN", "EOO", "Latitude", "ElevMin", "Prec", "PrecSeas", "Temp", "TempSeas", "HPD", "HPDMin", "HumanFootprint", "Accessibility", "Afrotropical", "Australasia", "Indo_malayan", "Nearctic", "Neotropical", "Oceania", "Palearctic")]


###Categorical variables transformed to orthogonal dummy variables##############
###Create the dummy variables
dummies <- dummyVars(" ~ .", data = BM_cat, fullRank=T)
head(dummies)
BM_cat_new <- data.frame(predict(dummies, newdata = BM_cat))
summary(BM_cat_new)

###Remove highly correlated variables####################
###Create a correlation matrix
descrCor <-cor(BM_numeric[,-1], y = NULL, use = "complete.obs",
               method = c("pearson", "kendall", "spearman"))
#descrCor <-  cor(BM_numeric[,-1])
summary(descrCor[upper.tri(descrCor)])

###Visualise correlation
corrplot(descrCor)

###Detect highly correlated variables and exclude them (but there are none)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .9) ###as this is zero there is no need to delete any of the columns
#filteredDescr <- BM_numeric[,-highlyCorDescr]
#descrCor2 <- cor(filteredDescr)
#summary(descrCor2[upper.tri(descrCor2)])

####Merge dataframes with variables of different categories back together
BM_new <- merge(BM_other,BM_cat_new, by = "ID")
###Look at the data
dim(BM_new) # check dataset
names(BM_new)
summary(BM_new)
typeof(BM_new)
str(BM_new)

###Removed variables with near-zero variance##############
nearZeroVar(BM_new, saveMetrics = FALSE) # returns a vector of integers corresponding to the column position of the problematic predictors
head(BM_new)
nzv <- nearZeroVar(BM_new)
check <- BM_new[, nzv]
head(check)
BM_new2 <- BM_new[, -nzv]
dim(BM_new2)
names(BM_new2)
summary(BM_new2)
typeof(BM_new2)
str(BM_new2)


###Partition the dataset
# Separate assessed species
assessed<- data.frame(subset(BM_new2, RedList.DD!=1))
assessed <- assessed[complete.cases(assessed), ]
#row.names(assessed)<- assessed$Binomial

assessed$binary<-"nonThr"
assessed$binary[assessed$RedList.VU == 1]="Thr"
assessed$binary[assessed$RedList.EN == 1]="Thr"
assessed$binary[assessed$RedList.CR == 1]="Thr"
dim(assessed)
head(assessed)

#assessed<- subset(assessed, select= -c(Binomial,RedList, binary, Genus))

# Setting nonassessed species aside
nonassessedspecies<- data.frame(subset (BM_new2, RedList.DD==1))
nonassessedspecies <- nonassessedspecies[complete.cases(nonassessedspecies), ]
nonassessedspecies$binary<-NA
head(nonassessedspecies)
dim(nonassessedspecies)

#row.names(nonassessedspecies)<- nonassessedspecies$Binomial
#nonassessedspecies<- subset(nonassessedspecies,select= -c(Binomial,RedList, Genus))

####Repeat model tuning parameters (takes about 25 min)
####Predicts 4 species threatened, 261 non-threatened, uses 63 as final mtry value
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = TRUE, 
  classProbs = TRUE)
#mtry <- sqrt(ncol(assessed[,2:36]))
    # run a random forest model
      start_time <- Sys.time()
      model <- train(as.factor(binary) ~ ., 
                     data = assessed, 
                     method = "rf",
                     tunelength=15,
                     ntree= 500,
                     trControl = model_control)
      end_time <- Sys.time()
      end_time - start_time

print(model)
plot(model)

###Predict status of assessed species
er_pred <- predict(model,assessed)
summary(er_pred)

x <- evalm(model)

## get roc curve plotted in ggplot2

x$roc

## get AUC and other metrics

x$stdres

summary(model$finalModel)


##  Interpreting probabilistic results####
###Predict status of DD species
er_pred <- predict(model, assessed,type="prob")
results_er_pred<-predict(model,assessed,type="prob")[,2]
summary(er_pred)
print(er_pred)


###Create a daframe with predicted and observed status
binary<- assessed$binary # define classes
binary2<- as.vector(binary)
temp <- merge(assessed, BM[, c("ID","RedList")], by = "ID")
RL_status<-temp$RedList
RL_status<-as.character(RL_status)
dataframe<- data.frame(cbind(results_er_pred, RL_status, binary2))
write.csv(dataframe, "Observed_predicted.csv") # insert correct path name to save results as .csv



###Plot probability classes
results_er_pred<- as.vector(results_er_pred)
dataframe<- data.frame(cbind(results_er_pred,binary2))
dataframe$results_er_pred<- as.numeric(as.character(dataframe$results_er_pred)) # as numeric changes the variable value. needs correction
par(mfrow=c(2,2))
hist(dataframe$results_er_pred, xlab="Probability of extinction", main=NULL, cex.main=1, ylim=c(0,50))
abline(v=0.42,col=3,lty=3) # set abline as correct probability threshold from results
hist(dataframe$results_er_pred[dataframe$binary2=="nonThr"], xlab="Probability of extinction", main=NULL, cex.main=1, ylim=c(0,50), xlim=c(0,1))
abline(v=0.42,col=3,lty=3)
hist(dataframe$results_er_pred[dataframe$binary2=="Thr"], xlab="Probability of extinction", main=NULL, cex.main=1, xlim=c(0,1))
abline(v=0.42,col=3,lty=3)

preds<- prediction(predictions=results_er_pred, labels=trainB) #  indicate which model results, e.g. results.rf
print(preds)

# AUC
AUC<- performance(preds, "auc")
AUC
myROC<- performance(preds, "tpr", "fpr")
myROCtable<- data.frame(cbind(myROC@alpha.values[[1]],myROC@y.values[[1]],myROC@x.values[[1]])) # creating dataframe with cutoff, tpr, fpr
# Youden's Index
myROCtable$X4<- myROCtable$X2 - myROCtable$X3
cutoff<- myROCtable$X1[which.max(myROCtable$X4)]
cutoff
Youden<- myROCtable$X4[which.max(myROCtable$X4)]
Youden
# ROC plot
par(mfrow=c(1,1))
plot(myROC, main=NULL, colorize=T, ylim= c(0,1))
# Confusion matrix
score<- ifelse(results_er_pred<=cutoff,"nonThr", "Thr")
confusion<- confusionMatrix(score,trainB, positive="Thr")
confusion







#####Start from here






#######################################################################################################################
############################ Comparing probabilities with IUCN Red List status ########################################

results<- results_er_pred
print(results)

pred.names<- as.vector(row.names(assessed))
assessed2<- data.frame(subset(BM, RL_status!="nonassessed"))
row.names(assessed2)<- assessed2$Genus_species
pred.frame<- assessed2[match(pred.names,assessed2$Genus_species),] # selecting rows for which we have predicted status
pred.frame$RL_status<- droplevels(pred.frame$RL_status) # drop unused levels
pred.frame<- cbind(pred.frame[,1:2], results)
pred.frame$RL_status <- factor(pred.frame$RL_status, levels = c("LC", "NT", "VU", "EN", "CR"), labels = c("LC", "NT", "VU", "EN", "CR")) 
plot(pred.frame$RL_status,pred.frame$results, xlab= "Red List category", cex.lab=0.9,
     ylab= "Predicted probability of threat",
     main=NULL, cex.main=0.9)
box<- boxplot(results~RL_status, labels=row.names(pred.frame),data=pred.frame, id.n=10)
identify(pred.frame$RedList,pred.frame$results.rf,pred.frame$Binomial,pos=F,plot=T)
write.table(pred.order,"insertpath") # insert correct path name to save table



#######################################################################################################################
############################################ Variable Importance ######################################################

random.forest<-randomForest(trainonassessed,trainB,ntreeTry=500, # changed to 500
                            mtry= rf.model$bestTune[1,1],replace=T,importance=TRUE)
rf.importance<- importance(random.forest, type=2)
rf.importance
rf.importance<- data.frame(rf.importance)
rf.importance$names<- row.names(rf.importance)
rf.importance<- rf.importance[order(rf.importance$MeanDecreaseGini, decreasing=FALSE),]
rf.importance

dotchart(rf.importance$MeanDecreaseGini, labels = row.names(rf.importance), xlab= "Mean Decrease in Gini Index", main="Variable importance in random forest model")

#Partial depence plots
par(mfrow=c(3,2))
partialPlot(random.forest,trainnonassessed, x.var= "Range", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="log(Range)", cex.main=1) # n.pt smoothes the graph - number of datapoints from which the function is computed
partialPlot(random.forest,trainnonassessed, x.var= "Isolation", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="Isolation index", cex.main=1) 
partialPlot(random.forest,trainnonassessed, x.var= "Human.Footprint.Index", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="log(Human Footprint Index)", cex.main=1) # n.pt smoothes the graph - number of datapoints from which the function is computed
partialPlot(random.forest,trainnonassessed, x.var= "ForestLoss", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="log(Forest loss)", cex.main=1) 
partialPlot(random.forest,trainnonassessed, x.var= "Population", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="log(Population density)" ,cex.main=1)
partialPlot(random.forest,trainnonassessed, x.var= "AreaProtected", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="log(Percentage prea protected)", cex.main=1)

# getTree(random.forest,k=1) - but can only see split for individual trees, no average.

################################################################################################################################
####################################### Predicting the status of nonassessed species ###########################################

resultsnonassessed<- predict(rf.model$finalModel, nonassessedspecies, type="prob")
resultsnonassessed<- resultsnonassessed[,2]
#print(resultsnonassessed)
scorenonassessed<- ifelse(resultsnonassessed<=cutoff,"nonThr", "Thr")
scorenonassessed<- as.character(scorenonassessed)
resultsnonassessed<-cbind(resultsnonassessed, scorenonassessed)
scorenonassessed<-as.factor(scorenonassessed)
summary(scorenonassessed)

write.csv(resultsnonassessed, file="insertpath") #insert correct path name to save prediction results for non assessed species




####Manually assign different values to tunelength###########################################
#1 - takes XX minutes
start_time <- Sys.time()
model_1 <- train(as.factor(binary) ~ ., 
               data = assessed, 
               method = "rf",
               tunelength=1,
               ntree= 500,
               trControl = model_control)
end_time <- Sys.time()
end_time - start_time

print(model)
plot(model)

###Predict status of DD species
er_pred_1 <- predict(model_1, nonassessedspecies)
summary(er_pred_1)

x_1 <- evalm(model_1)

## get roc curve plotted in ggplot2

x_1$roc

## get AUC and other metrics

x_1$stdres

summary(model$finalModel)




####Manually assign different values to tunelength###########################################
#10 - takes 18 minutes
start_time <- Sys.time()
model <- train(as.factor(binary) ~ ., 
               data = assessed, 
               method = "rf",
               tunelength=10,
               ntree= 500,
               trControl = model_control)
end_time <- Sys.time()
end_time - start_time

print(model)
plot(model)

###Predict status of DD species
er_pred <- predict(model, nonassessedspecies)
summary(er_pred)

x <- evalm(model)

## get roc curve plotted in ggplot2

x$roc

## get AUC and other metrics

x$stdres

summary(model$finalModel)

####Manually assign different values to tunelength###########################################
#20 - takes 19 minutes
start_time <- Sys.time()
model_20 <- train(as.factor(binary) ~ ., 
               data = assessed, 
               method = "rf",
               tunelength=20,
               ntree= 500,
               trControl = model_control)
end_time <- Sys.time()
end_time - start_time

print(model_20)
plot(model_20)

###Predict status of DD species
er_pred_20 <- predict(model_20, nonassessedspecies)
summary(er_pred_20)

x <- evalm(model_20)

## get roc curve plotted in ggplot2

x$roc

## get AUC and other metrics

x$stdres

summary(model_20$finalModel)


####Plot several ROCs together
res <- evalm(list(model,model_1,model_20),gnames='rf')

#######################################################################################################
###Trying to add cut off point based on Youden index##################################################
######################################################################################################
###Example from https://github.com/thie1e/cutpointr
mcp <- multi_cutpointr(suicide, class = suicide, pos_class = "yes", 
                       use_midpoints = TRUE, silent = TRUE) 
summary(mcp)


rocobj <-x$roc
coords(rocobj, "best")
coords(rocobj, x="best", input="threshold", best.method="youden")

data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b)
coords(rocobj, "best")
coords(rocobj, x="best", input="threshold", best.method="youden")


##################################################################################
####Try this to get stats

# NOT RUN {
set.seed(2444)
dat <- twoClassSim(500, intercept = -10)
table(dat$Class)

ctrl <- trainControl(method = "cv", 
                     classProbs = TRUE,
                     savePredictions = "all",
                     summaryFunction = twoClassSummary)

set.seed(2863)
mod <- train(Class ~ ., data = dat, 
             method = "rda",
             tuneLength = 4,
             metric = "ROC",
             trControl = ctrl)

resample_stats <- thresholder(mod, 
                              threshold = seq(.5, 1, by = 0.05), 
                              final = TRUE)

ggplot(resample_stats, aes(x = prob_threshold, y = J)) + 
  geom_point()
ggplot(resample_stats, aes(x = prob_threshold, y = Dist)) + 
  geom_point()
ggplot(resample_stats, aes(x = prob_threshold, y = Sensitivity)) + 
  geom_point() + 
  geom_point(aes(y = Specificity), col = "red")
# }











#######################################################################################################################
#########################  Random Forests Model for IUCN Red Listed and SRLI Species  #################################

# No partitioning #
trainassessed <- assessed
trainassessed <- trainassessed[complete.cases(trainassessed), ]
trainB <- trainassessed[,36]
#trainB <- trainassessed[,"binary"]  # define classes
#trainB <- as.numeric(trainB)

###Check that the length of TrainB matches the number of rows in trainassessed
dim(trainassessed)
length(trainB)

###If you need to check anything you can write out the files
write.csv(trainassessed, "trainassessed.csv")
write.csv(trainB, "trainB.csv")



# Control Parameters #
control.params<- trainControl(summaryFunction= twoClassSummary, selectionFunction= "best", 
                              method="repeatedcv", number=10, repeats=5, 
                              returnData=T,returnResamp="final",classProbs=TRUE)

# The Model #
#rf.model<- train(trainassessed, trainB, method="rf", metric="ROC", tuneLength= ncol(assessed),
#                trControl= control.params, ntrees = 500) # 500 trees works better than 100, or 1,000 trees

rf.model<- train(binary ~ ., data = trainassessed, method="rf", metric="ROC", tuneLength= ncol(assessed),
                 trControl= control.params, ntrees = 500) # 500 trees works better than 100, or 1,000 trees

rf.model<- train(trainassessed,trainB,method="rf", trControl= control.params, ntrees = 500) # 500 trees works better than 100, or 1,000 trees
print(rf.model)

#######################################################################################################################
##################################  Interpreting probabilistic results  ###############################################

results.rf<-predict(rf.model$finalModel,trainnonassessed,type="prob")[,2]
#print(results.rf)

binary2<- as.vector(binary)
RL_status<-as.character(RL_status)
dataframe<- data.frame(cbind(results.rf,RL_status, binary2))
write.csv(dataframe, file="insertpath") # insert correct path name to save results as .csv

results.rf<- as.vector(results.rf)
dataframe<- data.frame(cbind(results.rf,binary2))
dataframe$results.rf<- as.numeric(as.character(dataframe$results.rf)) # as numeric changes the variable value. needs correction
par(mfrow=c(2,2))
hist(dataframe$results.rf, xlab="Probability of extinction", main=NULL, cex.main=1, ylim=c(0,50))
abline(v=0.42,col=3,lty=3) # set abline as correct probability threshold from results
hist(dataframe$results.rf[dataframe$binary2=="nonThr"], xlab="Probability of extinction", main=NULL, cex.main=1, ylim=c(0,50), xlim=c(0,1))
abline(v=0.42,col=3,lty=3)
hist(dataframe$results.rf[dataframe$binary2=="Thr"], xlab="Probability of extinction", main=NULL, cex.main=1, xlim=c(0,1))
abline(v=0.42,col=3,lty=3)

preds<- prediction(predictions=results.rf, labels=trainB) #  indicate which model results, e.g. results.rf
print(preds)

# AUC
AUC<- performance(preds, "auc")
AUC
myROC<- performance(preds, "tpr", "fpr")
myROCtable<- data.frame(cbind(myROC@alpha.values[[1]],myROC@y.values[[1]],myROC@x.values[[1]])) # creating dataframe with cutoff, tpr, fpr
# Youden's Index
myROCtable$X4<- myROCtable$X2 - myROCtable$X3
cutoff<- myROCtable$X1[which.max(myROCtable$X4)]
cutoff
Youden<- myROCtable$X4[which.max(myROCtable$X4)]
Youden
# ROC plot
par(mfrow=c(1,1))
plot(myROC, main=NULL, colorize=T, ylim= c(0,1))
# Confusion matrix
score<- ifelse(results.rf<=cutoff,"nonThr", "Thr")
confusion<- confusionMatrix(score,trainB, positive="Thr")
confusion



#######################################################################################################################
############################ Comparing probabilities with IUCN Red List status ########################################

results<- results.rf
print(results)

pred.names<- as.vector(row.names(trainnonassessed))
assessed2<- data.frame(subset(BM, RL_status!="nonassessed"))
row.names(assessed2)<- assessed2$Genus_species
pred.frame<- assessed2[match(pred.names,assessed2$Genus_species),] # selecting rows for which we have predicted status
pred.frame$RL_status<- droplevels(pred.frame$RL_status) # drop unused levels
pred.frame<- cbind(pred.frame[,1:2], results)
pred.frame$RL_status <- factor(pred.frame$RL_status, levels = c("LC", "NT", "VU", "EN", "CR"), labels = c("LC", "NT", "VU", "EN", "CR")) 
plot(pred.frame$RL_status,pred.frame$results, xlab= "Red List category", cex.lab=0.9,
     ylab= "Predicted probability of threat",
     main=NULL, cex.main=0.9)
box<- boxplot(results~RL_status, labels=row.names(pred.frame),data=pred.frame, id.n=10)
identify(pred.frame$RedList,pred.frame$results.rf,pred.frame$Binomial,pos=F,plot=T)
write.table(pred.order,"insertpath") # insert correct path name to save table



#######################################################################################################################
############################################ Variable Importance ######################################################

random.forest<-randomForest(trainonassessed,trainB,ntreeTry=500, # changed to 500
                            mtry= rf.model$bestTune[1,1],replace=T,importance=TRUE)
rf.importance<- importance(random.forest, type=2)
rf.importance
rf.importance<- data.frame(rf.importance)
rf.importance$names<- row.names(rf.importance)
rf.importance<- rf.importance[order(rf.importance$MeanDecreaseGini, decreasing=FALSE),]
rf.importance

dotchart(rf.importance$MeanDecreaseGini, labels = row.names(rf.importance), xlab= "Mean Decrease in Gini Index", main="Variable importance in random forest model")

#Partial depence plots
par(mfrow=c(3,2))
partialPlot(random.forest,trainnonassessed, x.var= "Range", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="log(Range)", cex.main=1) # n.pt smoothes the graph - number of datapoints from which the function is computed
partialPlot(random.forest,trainnonassessed, x.var= "Isolation", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="Isolation index", cex.main=1) 
partialPlot(random.forest,trainnonassessed, x.var= "Human.Footprint.Index", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="log(Human Footprint Index)", cex.main=1) # n.pt smoothes the graph - number of datapoints from which the function is computed
partialPlot(random.forest,trainnonassessed, x.var= "ForestLoss", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="log(Forest loss)", cex.main=1) 
partialPlot(random.forest,trainnonassessed, x.var= "Population", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="log(Population density)" ,cex.main=1)
partialPlot(random.forest,trainnonassessed, x.var= "AreaProtected", n.pt=10, which.class= "Thr", main=NULL, ylab="Partial dependence", xlab="log(Percentage prea protected)", cex.main=1)

# getTree(random.forest,k=1) - but can only see split for individual trees, no average.

################################################################################################################################
####################################### Predicting the status of nonassessed species ###########################################

resultsnonassessed<- predict(rf.model$finalModel, nonassessedspecies, type="prob")
resultsnonassessed<- resultsnonassessed[,2]
#print(resultsnonassessed)
scorenonassessed<- ifelse(resultsnonassessed<=cutoff,"nonThr", "Thr")
scorenonassessed<- as.character(scorenonassessed)
resultsnonassessed<-cbind(resultsnonassessed, scorenonassessed)
scorenonassessed<-as.factor(scorenonassessed)
summary(scorenonassessed)

write.csv(resultsnonassessed, file="insertpath") #insert correct path name to save prediction results for non assessed species



#Not used
##Example
#install.packages("earth")
#library(earth)
#data(etitanic)
#dummies <- dummyVars(survived ~ ., data = etitanic)
#head(predict(dummies, newdata = etitanic))

#install.packages("dummies")
#library(dummies)
#BM.new <- dummy.data.frame(BM, sep = ".")
#students.new1 <- dummy.data.frame(students, names = c("State","Gender") , sep = ".")







# Test log transformation of all variables for scale on partial dependence plots later on
#BM$Range<- log(BM$Range)
#BM$Human.Footprint.Index<- log(BM$HFI+1)
#BM$ForestLoss<- log(BM$ForestLoss+1)
#BM$AreaProtected<- log(BM$AreaProtected+1)
#BM$Population<- log(BM$Population+1)
#dmy <- dummyVars(" ~ .", data = customers)

RedList<-assessed$RedList
assessed$binary<- factor((assessed$RedList.VU=="VU")+(assessed$RedList=="EN")+(assessed$RedList=="CR"))
levels(assessed$binary)[levels(assessed$binary)=="0"]<- "nonThr"
levels(assessed$binary)[levels(assessed$binary)=="1"]<- "Thr"
binary<- assessed$binary # define classes