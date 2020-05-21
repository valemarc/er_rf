##################################################################################################
### RANDOM FOREST MODEL FOR PREDICTING EXTINCTION RISK IN DATA DEFICIENT REPTILES ################
##################################################################################################

#Open data and look at the dataset
here()
#setwd("insertpath") # insert correct path name
BM <-read.csv(here("Bland_data.csv")) # insert correct path name
# check dataset
dim(BM) 
names(BM)
summary(BM)
typeof(BM)
str(BM)


###Descriptive stats for dataset
skimmed <- skim(BM)
write.csv(skimmed, "descr_stats.csv")

####################################################################################################
####Dataset preparation (complete dataset no imputation)############################################
####################################################################################################

#####Separate out different types of variables based on the transformation that we are going to do
###
BM_other <- BM[,c("ID", "Binomial", "Genus")]
BM_cat <- BM[,c("ID", "Order", "Family",  "Afrotropical", "Australasia", "Indo_malayan", "Nearctic", "Neotropical", "Oceania", "Palearctic", "RedList","ReproductiveMode", "TrophicGroup", "HabitatMode", "Continent")]
BM_numeric <- BM[,c("ID", "HabitatsIUCN", "EOO", "BodySize", "Latitude", "ElevMin", "Prec", "PrecSeas", "Temp", "TempSeas", "HPD", "HPDMin", "HumanFootprint", "Accessibility")]


####Look at variables
boxplot(BM_numeric[,-1], horizontal=TRUE)

###Density plots
par(mfrow=c(3, 5))
colnames <- dimnames(BM_numeric)[[2]]
for (i in 2:14) {
  BM_numeric_complete <- na.omit(BM_numeric)
  d <- density(BM_numeric_complete[,i])
  plot(d, type="n", main=colnames[i])
  polygon(d, col="red", border="gray")
}
dev.off()


###Categorical variables transformed to orthogonal dummy variables##############
###Create the dummy variables
###There are arguments for not following this step (increase in dimensionality) - https://roamanalytics.com/2016/10/28/are-categorical-variables-getting-lost-in-your-random-forests/
dummies <- dummyVars(" ~ .", data = BM_cat, fullRank=T) ###full-rank =T will create only (n-1) columns for a categorical column with n different levels
head(dummies)
BM_cat_new <- data.frame(predict(dummies, newdata = BM_cat))
summary(BM_cat_new)
#write.csv(BM_cat_new, "BM_cat_new.csv")

####Find out which variables have near zero variance in each of these subsets (categorical variables already transformed)
####Based on the settings from Bland and Bohm no variables are excluded

nzv_other <- nearZeroVar(BM_other, freqCut = 999/1, uniqueCut = 0.01, saveMetrics = FALSE)
nzv_other

nzv_cat <- nearZeroVar(BM_cat_new , freqCut = 999/1, uniqueCut = 0.01, saveMetrics = FALSE)
nzv_cat

nzv_numeric <- nearZeroVar(BM_numeric , freqCut = 999/1, uniqueCut = 0.01, saveMetrics = FALSE)
nzv_numeric


#Using caret's default setting s(freqCut = 19 and uniqueCut = 10, which is fairly conservative)
#nzv_other <- nearZeroVar(BM_other, freqCut = 19, uniqueCut = 10, saveMetrics = FALSE)
#nzv_other 
#head(BM_other[, nzv_other]) ####variable 3 (Order) would get excluded but it's a predictor in the final model it can't have been excluded

#nzv_cat <- nearZeroVar(BM_cat_new , freqCut = 19, uniqueCut = 10, saveMetrics = FALSE)
#nzv_cat 
#head(BM_cat_new[, nzv_cat])###excludes ReproductiveMode.Ovoviviparous, TrophicGroup.Herbivorous, TrophicGroup.Omnivorous, HabitatMode.Semi.aquatic, HabitatMode.Semi.arboreal, HabitatMode.Semi.fossorial, HabitatMode.Semi.saxicolous, Continent.Continent.Island

#nzv_numeric <- nearZeroVar(BM_numeric , freqCut = 19, uniqueCut = 10, saveMetrics = FALSE)
#nzv_numeric ###excludes none

###exclude them from dataset (only if chosing the stricter settings above)
#BM_binary_new <- BM_binary[, -nzv_binary]
#BM_cat_new <- BM_cat_new[, -nzv_cat]

###Centering and scaling numeric variables################
###Centre and scale variables
trainassessed_preprocess <- preProcess(BM_numeric[,-1], method = c("center", "scale"))
trainassessed_preprocess
trainassessed_preprocess$method ##check which variables have been transformed, ignored etc.

###As the function preprocess doesn't actually preprocess the data, we need to do this
BM_numeric_new <- predict(trainassessed_preprocess, BM_numeric)
head(BM_numeric_new)


###Density plots
par(mfrow=c(3, 5))
colnames <- dimnames(BM_numeric)[[2]]
for (i in 2:14) {
  BM_numeric_new <- na.omit(BM_numeric_new)
  d <- density(BM_numeric_new[,i])
  plot(d, type="n", main=colnames[i])
  polygon(d, col="red", border="gray")
}
dev.off()

###Remove highly correlated variables####################
###Create a correlation matrix
descrCor <-cor(BM_numeric_new[,-1], y = NULL, use = "complete.obs",
               method = c("pearson", "kendall", "spearman"))
#descrCor <-  cor(BM_numeric[,-1])
summary(descrCor[upper.tri(descrCor)])

###Visualise correlation
corrplot(descrCor)

###Detect highly correlated variables and exclude them (but there are none)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .9)
highlyCorDescr ###as this is zero there is no need to delete any of the columns
#filteredDescr <- BM_numeric[,-highlyCorDescr]
#descrCor2 <- cor(filteredDescr)
#summary(descrCor2[upper.tri(descrCor2)])


####Put the variables back in a single dataset (except for BM_other, as it will throw an error if passes through the impute function)
####Merge dataframes with variables of different categories back together
#BM_final <- merge(BM_binary_new, by = "ID")
###Look at the data
#dim(BM_final) # check dataset

BM_final <- merge(BM_cat_new, BM_numeric_new, by = "ID")
dim(BM_final) # check dataset


####Create dataset used for the model with only complete species
BM_final_complete <- merge(BM_final,BM_other, by = "ID")
dim(BM_final_complete) # check dataset

###Exclude non complete cases
BM_final_complete <- BM_final_complete[complete.cases(BM_final_complete), ]
dim(BM_final_complete)

###Partition the dataset
# Separate assessed species
data_suff_complete <- data.frame(subset(BM_final_complete, RedList.DD!=1))
dim(data_suff_complete)

data_suff_complete$binary<-"nonThr"
data_suff_complete$binary[data_suff_complete$RedList.VU == 1]="Thr"
data_suff_complete$binary[data_suff_complete$RedList.EN == 1]="Thr"
data_suff_complete$binary[data_suff_complete$RedList.CR == 1]="Thr"
dim(data_suff_complete)
head(data_suff_complete)

# Setting DD species aside
data_def_complete <- data.frame(subset (BM_final_complete, RedList.DD==1))
data_def_complete$binary<-NA
head(data_def_complete)
dim(data_def_complete)


##################################################################################################
#Spliting training set into two parts based on outcome: 80% and 20%
index <- createDataPartition(data_suff_complete$binary, p=0.8, list=FALSE)
trainingSet <- data_suff_complete[ index,]
validationtSet <- data_suff_complete[-index,]

###Exclude variables that are not used as predictors (binomial and threat status)
trainingSet_compl_a <- subset(trainingSet, select = -c(ID,Binomial,Genus, RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))
validationtSet_compl_a <- subset(validationtSet, select = -c(ID,Binomial,Genus, RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))
predSet_complete_a <- subset(data_def_complete, select = -c(ID,Binomial,Genus, RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))


#############Run the model#############################################
###Set model parameters
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = "final", 
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

    # run a random forest model
    start_time <- Sys.time()
    model_15_a <- train(binary ~ ., 
                      data = trainingSet_compl_a, 
                      method = "rf",
                      tuneLength=15, ####see https://rpubs.com/phamdinhkhanh/389752
                      ntree= 500,
                      trControl = model_control, 
                      metric="ROC")
    end_time <- Sys.time()
    end_time - start_time


####Print 
print(model_15_a)
plot(model_15_a)
print(model_15_a$finalModel) ####this gives the confusion matrix but doesn't calculate the accuracy 
summary(model_15_a$finalModel)
model_15_a$finalModel$confusion
model_15_a$finalModel$ntree
model_15_a$finalModel$mtry

###Average confusion matrix across 10 folds
confusionMatrix.train(model_15_a) #average across 10 samples #https://stats.stackexchange.com/questions/118568/random-forest-confusion-matrix

# ###NOT WORKING
# #Overall confusion matrix (only works if SavePred=FINAL in train)
# #Predictions are in:
# model_15_a$pred
# #sorted as per CV fols, to sort as in original data frame:
# model_15_a$pred[order(model_15_a$pred$rowIndex),2]
# #to obtain a confusion matrix:
# caret::confusionMatrix(model_15_a$pred[order(model_15_a$pred$rowIndex),2], trainingSet_compl$binary)
# 
# 
# ###Confusion matrix on training data (it doesn't print it out, but doing the counts myself if gives perfect accuracy)
# rfPred <- predict.train(model_15_a, trainingSet_compl, type = "raw")
# confusionMatrix(rfPred, trainingSet_compl$binary) 
# summary(as.factor(trainingSet_compl$binary))
# summary(rfPred)


###Using evalm
###as per https://stackoverflow.com/questions/31138751/roc-curve-from-training-data-in-caret
## run MLeval
res <- evalm(model_15_a)
## get ROC
res$roc
## get calibration curve
res$cc
## get precision recall gain curve
res$prg

###NOT WORKING
###Plot tree
source ("tree_func.R")
tree_num <- which(model_15_a$finalModel$forest$ndbigtree == max(model_15_a$finalModel$forest$ndbigtree))
tree_func(final_model = model_15_a$finalModel, tree_num)

#Warning messages:
#  1: Duplicated aesthetics after name standardisation: na.rm 
#2: Duplicated aesthetics after name standardisation: na.rm 
#3: Duplicated aesthetics after name standardisation: na.rm 
#4: Removed 98 rows containing missing values (geom_text_repel). 
#5: Removed 98 rows containing missing values (geom_label). 
#6: Removed 97 rows containing missing values (geom_label_repel). 


#####ROC curve for final average value
for_lift <- data.frame(binary = model_15_a$pred$obs, rf = model_15_a$pred$Thr)
lift_obj <- lift(binary ~ rf, data = for_lift, class = "Thr")   ###to check

# Plot ROC ----------------------------------------------------------------
ggplot(lift_obj$data) +
  geom_line(aes(1 - Sp, Sn, color = liftModelVar)) +
  scale_color_discrete(guide = guide_legend(title = "method"))#

#Get a ROC curve for each fold#############################################
for_lift <- data.frame(binary = model_15_a$pred$obs, rf = model_15_a$pred$Thr, resample = model_15_a$pred$Resample)
lift_df <-  data.frame()
for (fold in unique(for_lift$resample)) {
  fold_df <- dplyr::filter(for_lift, resample == fold)
  lift_obj_data <- lift(binary ~ rf, data = fold_df, class = "Thr")$data
  lift_obj_data$fold = fold
  lift_df = rbind(lift_df, lift_obj_data)
}

lift_obj <- lift(binary ~ rf, data = for_lift, class = "Thr")

# Plot ROC
ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))


###Overlap
ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))+
  geom_line(data = lift_obj$data, aes(1 - Sp, Sn, color = liftModelVar), lwd=1, colour="black")

###Calculate accuracy within individual samples
ddply(model_15_a$pred, "Resample", summarise,
      accuracy = Accuracy(pred, obs))

###Youden index
#thresholder - This function uses the resampling results from a train object to generate performance statistics
#over a set of probability thresholds for two-class problems.
resample_stats_a <- thresholder(model_15_a,
                              threshold = seq(.1, 1, by = 0.0005),
                              final = TRUE,
                              statistics = "J")
ggplot(resample_stats_a, aes(x = prob_threshold, y = J)) +
  geom_point()

###Extract optimal probability threshold (that maximises Youden)
threshold_a <- resample_stats_a[which.max(resample_stats_a$J),]
threshold_a <- threshold_a$prob_threshold
threshold_a

####To finish
####Recalculate class (Thr/Non threatened) based on new threshold (0.86)
#pred_52 <- subset(model_15$pred, mtry==52)
#label <- ifelse(pred_52$Thr > threshold, 'Thr', 'nonThr')
#summary(label)
#levels(label)
#count(label == "Thr")
#count(label == "nonThr")

# ####To recalculate class (Thr/Non threatened) based on new threshold using the final model
# pred <- predict(model_15, newdata = validationtSet_compl, type="prob")
# label_final <- as.factor(ifelse(pred$Thr > threshold, 'Thr', 'nonThr'))
# summary(label_final)
# ConfusionMatrix(label_final, validationtSet_compl$binary)
# #confusionMatrix(label_final, model_15$final$y)$overall
# 
# 
# ###Not sure why confusionMatrix doesn't work properly
# #Predictions
# #pred <- predict(model_15, newdata = validationtSet_compl,type="raw")
# #table(pred)
# #confusionMatrix(as.factor(pred),as.factor(validationtSet_compl$binary))
# 
# #levels(as.factor(validationtSet_compl$binary))
# #levels(pred)
# 
# ###Confusion matrix on training data (it doesn't print it out)
# #rfPred <- predict.train(model_15, validationtSet_compl, type = "raw")
# #confusionMatrix(rfPred, validationtSet_compl$binary) 
# #summary(as.factor(validationtSet_compl$binary))
# #summary(rfPred)
# 
# ###Compare with Lucie's threshold
# label_final_Bland <- as.factor(ifelse(pred$Thr > 0.614, 'Thr', 'nonThr'))
# summary(label_final_Bland)
# ConfusionMatrix(validationtSet_compl$binary,label_final_Bland)
# 
# ###To predict status of DD species
# pred_DD <- predict(model_15, newdata = predSet_complete, type="prob")
# label_DD <- as.factor(ifelse(pred_DD$Thr > threshold, 'Thr', 'nonThr'))
# summary(label_DD)
# 
# ###Compare with Lucie's threshold
# label_DD_Bland <- as.factor(ifelse(pred_DD$Thr > 0.614, 'Thr', 'nonThr'))
# summary(label_DD_Bland )

####Variable importance
varImp <- varImp(model_15_a, scale = FALSE)
varImp
plot(varImp, top=20)

ggplot(pred_DD, aes(x=Thr)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept=threshold_a, linetype="dotted", color="red")+
  geom_vline(xintercept=0.5, linetype="dotted", color="blue")


#######################################################################################
#############Run the model including genus in predictors###############################
#######################################################################################

##################################################################################################
###Exclude variables that are not used as predictors (binomial and threat status)
trainingSet_compl_b <- subset(trainingSet, select = -c(ID,Binomial,RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))
validationtSet_compl_b <- subset(validationtSet, select = -c(ID,Binomial,RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))
predSet_complete_b <- subset(data_def_complete, select = -c(ID,Binomial,RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))


###Set model parameters
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = "final", 
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

    # run a random forest model
    start_time <- Sys.time()
    model_15_b <- train(binary ~ ., 
                      data = trainingSet_compl_b, 
                      method = "rf",
                      tuneLength=15, ####see https://rpubs.com/phamdinhkhanh/389752
                      ntree= 500,
                      trControl = model_control, 
                      metric="ROC")
    end_time <- Sys.time()
    end_time - start_time


####Print 
print(model_15_b)
plot(model_15_b)
print(model_15_b$finalModel) ####this gives the confusion matrix but doesn't calculate the accuracy 
summary(model_15_b$finalModel)
model_15_b$finalModel$confusion
model_15_b$finalModel$ntree
model_15_b$finalModel$mtry


confusionMatrix.train(model_15_b) #average across 10 samples #https://stats.stackexchange.com/questions/118568/random-forest-confusion-matrix

###Using evalm
###as per https://stackoverflow.com/questions/31138751/roc-curve-from-training-data-in-caret
## run MLeval
res_b <- evalm(model_15_b)
## get ROC
res_b$roc
## get calibration curve
res_b$cc
## get precision recall gain curve
res_b$prg


###Plot tree
source ("tree_func.R")
tree_num <- which(model_15_b$finalModel$forest$ndbigtree == max(model_15_b$finalModel$forest$ndbigtree))
tree_func(final_model = model_15_b$finalModel, tree_num)


#Warning messages:
#  1: Duplicated aesthetics after name standardisation: na.rm 
#2: Duplicated aesthetics after name standardisation: na.rm 
#3: Duplicated aesthetics after name standardisation: na.rm 
#4: Removed 98 rows containing missing values (geom_text_repel). 
#5: Removed 98 rows containing missing values (geom_label). 
#6: Removed 97 rows containing missing values (geom_label_repel). 


#####ROC curve for final average value
for_lift <- data.frame(binary = model_15_b$pred$obs, rf = model_15_b$pred$Thr)
lift_obj <- lift(binary ~ rf, data = for_lift, class = "Thr")   ###to check

# Plot ROC ----------------------------------------------------------------
ggplot(lift_obj$data) +
  geom_line(aes(1 - Sp, Sn, color = liftModelVar)) +
  scale_color_discrete(guide = guide_legend(title = "method"))#

#Get a ROC curve for each fold#############################################
for_lift <- data.frame(binary = model_15_b$pred$obs, rf = model_15_b$pred$Thr, resample = model_15_b$pred$Resample)
lift_df <-  data.frame()
for (fold in unique(for_lift$resample)) {
  fold_df <- dplyr::filter(for_lift, resample == fold)
  lift_obj_data <- lift(binary ~ rf, data = fold_df, class = "Thr")$data
  lift_obj_data$fold = fold
  lift_df = rbind(lift_df, lift_obj_data)
}

lift_obj <- lift(binary ~ rf, data = for_lift, class = "Thr")

# Plot ROC
ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))


###Overlap
ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))+
  geom_line(data = lift_obj$data, aes(1 - Sp, Sn, color = liftModelVar), lwd=1, colour="black")

###Calculate accuracy within individual samples
ddply(model_15_b$pred, "Resample", summarise,
      accuracy = Accuracy(pred, obs))

###Youden index
#thresholder - This function uses the resampling results from a train object to generate performance statistics
#over a set of probability thresholds for two-class problems.
resample_stats_b <- thresholder(model_15_b,
                                threshold = seq(.1, 1, by = 0.0005),
                                final = TRUE,
                                statistics = "J")
ggplot(resample_stats_b, aes(x = prob_threshold, y = J)) +
  geom_point()

###Extract optimal probability threshold (that maximises Youden)
threshold_b <- resample_stats_b[which.max(resample_stats_b$J),]
threshold_b <- threshold_b$prob_threshold
threshold_b

####To finish
####Recalculate class (Thr/Non threatened) based on new threshold (0.86)
#pred_52 <- subset(model_15$pred, mtry==52)
#label <- ifelse(pred_52$Thr > threshold, 'Thr', 'nonThr')
#summary(label)
#levels(label)
#count(label == "Thr")
#count(label == "nonThr")

# ####To recalculate class (Thr/Non threatened) based on new threshold (0.86) using the final model
# pred <- predict(model_15_b, newdata = validationtSet_compl_a, type="prob")
# label_final <- as.factor(ifelse(pred$Thr > threshold, 'Thr', 'nonThr'))
# summary(label_final)
# ConfusionMatrix(label_final, validationtSet_compl_a$binary)
# 
# 
# ###Compare with Lucie's threshold
# label_final_Bland <- as.factor(ifelse(pred$Thr > 0.614, 'Thr', 'nonThr'))
# summary(label_final_Bland)
# ConfusionMatrix(model_15_b$final$y,label_final_Bland)
# 
# ###To predict status of DD species
# pred_DD <- predict(model_15_b, newdata = predSet_complete_a, type="prob")
# label_DD <- as.factor(ifelse(pred_DD$Thr > threshold, 'Thr', 'nonThr'))
# summary(label_DD)
# 
# ###Compare with Lucie's threshold
# label_DD_Bland <- as.factor(ifelse(pred_DD$Thr > 0.614, 'Thr', 'nonThr'))
# summary(label_DD_Bland )

####Variable importance
varImp <- varImp(model_15_b, scale = FALSE)
varImp
plot(varImp, top=20)

ggplot(pred_DD, aes(x=Thr)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept=threshold_b, linetype="dotted", color="red")+
  geom_vline(xintercept=0.5, linetype="dotted", color="blue")




##########################################################################################
####Alternative dataset with imputed values################################################
###########################################################################################
#Open data and look at the dataset
here()
#setwd("insertpath") # insert correct path name
BM <-read.csv(here("Bland_data.csv")) # insert correct path name
# check dataset
dim(BM) 

#####Separate out different types of variables based on the transformation that we are going to do
###
BM_other <- BM[,c("ID", "Binomial", "Genus")]
BM_cat <- BM[,c("ID", "Order", "Family",  "Afrotropical", "Australasia", "Indo_malayan", "Nearctic", "Neotropical", "Oceania", "Palearctic", "RedList","ReproductiveMode", "TrophicGroup", "HabitatMode", "Continent")]
BM_numeric <- BM[,c("ID", "HabitatsIUCN", "EOO", "BodySize", "Latitude", "ElevMin", "Prec", "PrecSeas", "Temp", "TempSeas", "HPD", "HPDMin", "HumanFootprint", "Accessibility")]


###Categorical variables transformed to orthogonal dummy variables##############
###Create the dummy variables
###There are arguments for not following this step (increase in dimensionality) - https://roamanalytics.com/2016/10/28/are-categorical-variables-getting-lost-in-your-random-forests/
dummies <- dummyVars(" ~ .", data = BM_cat, fullRank=T) ###full-rank =T will create only (n-1) columns for a categorical column with n different levels
head(dummies)
BM_cat_new <- data.frame(predict(dummies, newdata = BM_cat))
summary(BM_cat_new)
#write.csv(BM_cat_new, "BM_cat_new.csv")

####Find out which variables have near zero variance in each of these subsets (categorical variables already transformed)
####Based on the settings from Bland and Bohm no variables are excluded

nzv_other <- nearZeroVar(BM_other, freqCut = 999/1, uniqueCut = 0.01, saveMetrics = FALSE)
nzv_other

nzv_cat <- nearZeroVar(BM_cat_new , freqCut = 999/1, uniqueCut = 0.01, saveMetrics = FALSE)
nzv_cat

nzv_numeric <- nearZeroVar(BM_numeric , freqCut = 999/1, uniqueCut = 0.01, saveMetrics = FALSE)
nzv_numeric


#Using caret's default setting s(freqCut = 19 and uniqueCut = 10, which is fairly conservative)
#nzv_other <- nearZeroVar(BM_other, freqCut = 19, uniqueCut = 10, saveMetrics = FALSE)
#nzv_other 
#head(BM_other[, nzv_other]) ####variable 3 (Order) would get excluded but it's a predictor in the final model it can't have been excluded

#nzv_cat <- nearZeroVar(BM_cat_new , freqCut = 19, uniqueCut = 10, saveMetrics = FALSE)
#nzv_cat 
#head(BM_cat_new[, nzv_cat])###excludes ReproductiveMode.Ovoviviparous, TrophicGroup.Herbivorous, TrophicGroup.Omnivorous, HabitatMode.Semi.aquatic, HabitatMode.Semi.arboreal, HabitatMode.Semi.fossorial, HabitatMode.Semi.saxicolous, Continent.Continent.Island

#nzv_numeric <- nearZeroVar(BM_numeric , freqCut = 19, uniqueCut = 10, saveMetrics = FALSE)
#nzv_numeric ###excludes none

###exclude them from dataset (only if chosing the stricter settings above)
#BM_binary_new <- BM_binary[, -nzv_binary]
#BM_cat_new <- BM_cat_new[, -nzv_cat]

###Centering and scaling numeric variables################
###Centre and scale variables
trainassessed_preprocess <- preProcess(BM_numeric[,-1], method = c("center", "scale"))
trainassessed_preprocess
trainassessed_preprocess$method ##check which variables have been transformed, ignored etc.

###As the function preprocess doesn't actually preprocess the data, we need to do this
BM_numeric_new <- predict(trainassessed_preprocess, BM_numeric)
head(BM_numeric_new)


###Remove highly correlated variables####################
###Create a correlation matrix
descrCor <-cor(BM_numeric_new[,-1], y = NULL, use = "complete.obs",
               method = c("pearson", "kendall", "spearman"))
#descrCor <-  cor(BM_numeric[,-1])
summary(descrCor[upper.tri(descrCor)])

###Visualise correlation
corrplot(descrCor)

###Detect highly correlated variables and exclude them (but there are none)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .9)
highlyCorDescr ###as this is zero there is no need to delete any of the columns
#filteredDescr <- BM_numeric[,-highlyCorDescr]
#descrCor2 <- cor(filteredDescr)
#summary(descrCor2[upper.tri(descrCor2)])


####Put the variables back in a single dataset (except for BM_other, as it will throw an error if passes through the impute function)
####Merge dataframes with variables of different categories back together
#BM_final <- merge(BM_binary_new, by = "ID")
###Look at the data
#dim(BM_final) # check dataset
BM_final <- merge(BM_cat_new, BM_numeric_new, by = "ID")
dim(BM_final) # check dataset


####Impute missing values for the model with imputed values
BM_imp <- missForest(BM_final, verbose = TRUE) #https://stats.stackexchange.com/questions/49270/imputation-with-random-forests    
#BM_numeric_new_imp <- missForest(BM_numeric, verbose = TRUE)
BM_final_imputed <- merge(BM_other,BM_imp$ximp, by = "ID")

###Partition the dataset
# Separate assessed species
data_suff_imputed <- data.frame(subset(BM_final_imputed, RedList.DD!=1))
dim(data_suff_imputed)

data_suff_imputed$binary<-"nonThr"
data_suff_imputed$binary[data_suff_imputed$RedList.VU == 1]="Thr"
data_suff_imputed$binary[data_suff_imputed$RedList.EN == 1]="Thr"
data_suff_imputed$binary[data_suff_imputed$RedList.CR == 1]="Thr"
dim(data_suff_imputed)
head(data_suff_imputed)

# Setting nonassessed species aside
data_def_imputed <- data.frame(subset (BM_final_imputed, RedList.DD==1))
data_def_imputed$binary<-NA
head(data_def_imputed)
dim(data_def_imputed)

###Create a smaller sample for testing in case it's useful
#assessed_small <- sample_n(assessed, 200)
#summary(assessed_small)
##################################################################################################
#Spliting training set into two parts based on outcome: 80% and 20%
index <- createDataPartition(data_suff_imputed$binary, p=0.8, list=FALSE)
trainingSet_imp <- data_suff_imputed[ index,]
validationtSet_imp <- data_suff_imputed[-index,]

###Exclude variables that are not used as predictors (binomial and threat status)
trainingSet_imp_c <- subset(trainingSet_imp, select = -c(ID,Binomial,Genus, RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))
validationtSet_imp_c <- subset(validationtSet_imp, select = -c(ID,Binomial,Genus, RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))
predSet_imp_c <- subset(data_def_imputed, select = -c(ID,Binomial,Genus, RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))




############################################################################################
###Set model parameters with savepredictions=="Final"######################################
###########################################################################################
model_control_final <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = "final", 
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

# run a random forest model
start_time <- Sys.time()
model_15_c <- train(binary ~ ., 
                  data = trainingSet_imp_c, 
                  method = "rf",
                  tuneLength=15, ####see https://rpubs.com/phamdinhkhanh/389752
                  ntree= 500,
                  trControl = model_control_final, 
                  metric="ROC")
end_time <- Sys.time()
end_time - start_time


####Print 
print(model_15_c)
plot(model_15_c)
print(model_15_c$finalModel) ####this gives the confusion matrix but doesn't calculate the accuracy 
summary(model_15_c$finalModel)
model_15_c$finalModel$confusion
model_15_c$finalModel$ntree
model_15_c$finalModel$mtry

confusionMatrix.train(model_15_c) #average across 10 samples #https://stats.stackexchange.com/questions/118568/random-forest-confusion-matrix

###Using evalm
###as per https://stackoverflow.com/questions/31138751/roc-curve-from-training-data-in-caret
## run MLeval
res <- evalm(model_15_c)
## get ROC
res$roc
## get calibration curve
res$cc
## get precision recall gain curve
res$prg


###Average confusion matrix across 10 folds
confusionMatrix.train(model_15_c) #average across 10 samples #https://stats.stackexchange.com/questions/118568/random-forest-confusion-matrix

###Plot tree
source ("tree_func.R")
tree_num <- which(model_15_c$finalModel$forest$ndbigtree == max(model_15_c$finalModel$forest$ndbigtree))
tree_func(final_model = model_15_c$finalModel, tree_num)

# Error in dimnames(tree) <- list(1:nrow(tree), c("left daughter", "right daughter",  : 
#                                                   length of 'dimnames' [2] not equal to array extent
#                                                 In addition: Warning messages:
#                                                   1: In if (k > rfobj$ntree) { :
#                                                       the condition has length > 1 and only the first element will be used
#                                                     2: In cbind(rfobj$forest$treemap[, , k], rfobj$forest$bestvar[, k],  :
#                                                                   number of rows of result is not a multiple of vector length (arg 1)
#                                                                 3: In 1:rfobj$forest$ndbigtree[k] :
#                                                                   
#                                                                   Error in dimnames(tree) <- list(1:nrow(tree), c("left daughter", "right daughter",  : 
#                                                                                                                     length of 'dimnames' [2] not equal to array extent 


#Warning messages:
#  1: Duplicated aesthetics after name standardisation: na.rm 
#2: Duplicated aesthetics after name standardisation: na.rm 
#3: Duplicated aesthetics after name standardisation: na.rm 
#4: Removed 98 rows containing missing values (geom_text_repel). 
#5: Removed 98 rows containing missing values (geom_label). 
#6: Removed 97 rows containing missing values (geom_label_repel). 


#####ROC curve for final average value
for_lift <- data.frame(binary = model_15_c$pred$obs, rf = model_15_c$pred$Thr)
lift_obj <- lift(binary ~ rf, data = for_lift, class = "Thr")   ###to check

# Plot ROC ----------------------------------------------------------------
ggplot(lift_obj$data) +
  geom_line(aes(1 - Sp, Sn, color = liftModelVar)) +
  scale_color_discrete(guide = guide_legend(title = "method"))#

#Get a ROC curve for each fold#############################################
for_lift <- data.frame(binary = model_15_c$pred$obs, rf = model_15_c$pred$Thr, resample = model_15_c$pred$Resample)
lift_df <-  data.frame()
for (fold in unique(for_lift$resample)) {
  fold_df <- dplyr::filter(for_lift, resample == fold)
  lift_obj_data <- lift(binary ~ rf, data = fold_df, class = "Thr")$data
  lift_obj_data$fold = fold
  lift_df = rbind(lift_df, lift_obj_data)
}

lift_obj <- lift(binary ~ rf, data = for_lift, class = "Thr")

# Plot ROC
ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))


###Overlap
ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))+
  geom_line(data = lift_obj$data, aes(1 - Sp, Sn, color = liftModelVar), lwd=1, colour="black")

###Calculate accuracy within individual samples
ddply(model_15_c$pred, "Resample", summarise,
      accuracy = Accuracy(pred, obs))

###Youden index
#thresholder - This function uses the resampling results from a train object to generate performance statistics
#over a set of probability thresholds for two-class problems.
resample_stats_c <- thresholder(model_15_c,
                                threshold = seq(.1, 1, by = 0.0005),
                                final = TRUE,
                                statistics = "J")
ggplot(resample_stats_c, aes(x = prob_threshold, y = J)) +
  geom_point()

###Extract optimal probability threshold (that maximises Youden)
threshold_c <- resample_stats_c[which.max(resample_stats_c$J),]
threshold_c <- threshold_c$prob_threshold
threshold_c

####To finish
####Recalculate class (Thr/Non threatened) based on new threshold (0.86)
#pred_52 <- subset(model_15_2$pred, mtry==52)
#label <- ifelse(pred_52$Thr > threshold, 'Thr', 'nonThr')
#summary(label)
#levels(label)
#count(label == "Thr")
#count(label == "nonThr")

# ####To recalculate class (Thr/Non threatened) based on new threshold (0.762) using the final model
# pred <- predict(model_15_c, newdata = validationtSet_imp_c, type="prob")
# label_final <- as.factor(ifelse(pred$Thr > threshold, 'Thr', 'nonThr'))
# summary(label_final)
# ConfusionMatrix(label_final, validationtSet_imp$binary)
# 
# ###Compare with Lucie's threshold
# label_final_Bland <- as.factor(ifelse(pred$Thr > 0.614, 'Thr', 'nonThr'))
# summary(label_final_Bland)
# ConfusionMatrix(validationtSet_imp$binary,label_final_Bland)
# 
# ###To predict status of DD species
# pred_DD <- predict(model_15_c, newdata = predSet_imp, type="prob")
# label_DD <- as.factor(ifelse(pred_DD$Thr > threshold, 'Thr', 'nonThr'))
# summary(label_DD)
# 
# ###Compare with Lucie's threshold
# label_DD_Bland <- as.factor(ifelse(pred_DD$Thr > 0.614, 'Thr', 'nonThr'))
# summary(label_DD_Bland )

####Variable importance
varImp <- varImp(model_15_c, scale = FALSE)
varImp
plot(varImp, top=20)

ggplot(pred_DD, aes(x=Thr)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept=threshold_c, linetype="dotted", color="red")+
  geom_vline(xintercept=0.5, linetype="dotted", color="blue")



##################################################################################################
#####Imputed dataset including genus###############################################################
###################################################################################################
#Spliting training set into two parts based on outcome: 80% and 20%
index <- createDataPartition(data_suff_imputed$binary, p=0.8, list=FALSE)
trainingSet_imp <- data_suff_imputed[ index,]
validationtSet_imp <- data_suff_imputed[-index,]

###Exclude variables that are not used as predictors (binomial and threat status)
trainingSet_imp_d <- subset(trainingSet_imp, select = -c(ID,Binomial,RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))
validationtSet_imp_d <- subset(validationtSet_imp, select = -c(ID,Binomial,RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))
predSet_imp_d <- subset(data_def_imputed, select = -c(ID,Binomial,RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))




  ############################################################################################
  ###Set model parameters with savepredictions=="Final"######################################
  ###########################################################################################
  model_control_final <- trainControl(## 10-fold CV
    method = "cv",
    number = 10,
    savePredictions = "final", 
    classProbs = TRUE,
    summaryFunction = twoClassSummary)
  
  # run a random forest model
  start_time <- Sys.time()
  model_15_d <- train(binary ~ ., 
                      data = trainingSet_imp_d, 
                      method = "rf",
                      tuneLength=15, ####see https://rpubs.com/phamdinhkhanh/389752
                      ntree= 500,
                      trControl = model_control_final, 
                      metric="ROC")
  end_time <- Sys.time()
  end_time - start_time


####Print 
print(model_15_d)
plot(model_15_d)
print(model_15_d$finalModel) ####this gives the confusion matrix but doesn't calculate the accuracy 
summary(model_15_d$finalModel)
model_15_d$finalModel$confusion
model_15_d$finalModel$ntree
model_15_d$finalModel$mtry


###Using evalm
###as per https://stackoverflow.com/questions/31138751/roc-curve-from-training-data-in-caret
## run MLeval
res <- evalm(model_15_d)
## get ROC
res$roc
## get calibration curve
res$cc
## get precision recall gain curve
res$prg


###Average confusion matrix across 10 folds
confusionMatrix.train(model_15_d) #average across 10 samples #https://stats.stackexchange.com/questions/118568/random-forest-confusion-matrix

###Plot tree
source ("tree_func.R")
tree_num <- which(model_15_d$finalModel$forest$ndbigtree == max(model_15_d$finalModel$forest$ndbigtree))
tree_func(final_model = model_15_d$finalModel, tree_num)

# Error in dimnames(tree) <- list(1:nrow(tree), c("left daughter", "right daughter",  : 
#                                                   length of 'dimnames' [2] not equal to array extent
#                                                 In addition: Warning messages:
#                                                   1: In if (k > rfobj$ntree) { :
#                                                       the condition has length > 1 and only the first element will be used
#                                                     2: In cbind(rfobj$forest$treemap[, , k], rfobj$forest$bestvar[, k],  :
#                                                                   number of rows of result is not a multiple of vector length (arg 1)
#                                                                 3: In 1:rfobj$forest$ndbigtree[k] :
#                                                                   
#                                                                   Error in dimnames(tree) <- list(1:nrow(tree), c("left daughter", "right daughter",  : 
#                                                                                                                     length of 'dimnames' [2] not equal to array extent 


#Warning messages:
#  1: Duplicated aesthetics after name standardisation: na.rm 
#2: Duplicated aesthetics after name standardisation: na.rm 
#3: Duplicated aesthetics after name standardisation: na.rm 
#4: Removed 98 rows containing missing values (geom_text_repel). 
#5: Removed 98 rows containing missing values (geom_label). 
#6: Removed 97 rows containing missing values (geom_label_repel). 


#####ROC curve for final average value
for_lift <- data.frame(binary = model_15_d$pred$obs, rf = model_15_d$pred$Thr)
lift_obj <- lift(binary ~ rf, data = for_lift, class = "Thr")   ###to check

# Plot ROC ----------------------------------------------------------------
ggplot(lift_obj$data) +
  geom_line(aes(1 - Sp, Sn, color = liftModelVar)) +
  scale_color_discrete(guide = guide_legend(title = "method"))#

#Get a ROC curve for each fold#############################################
for_lift <- data.frame(binary = model_15_d$pred$obs, rf = model_15_d$pred$Thr, resample = model_15_d$pred$Resample)
lift_df <-  data.frame()
for (fold in unique(for_lift$resample)) {
  fold_df <- dplyr::filter(for_lift, resample == fold)
  lift_obj_data <- lift(binary ~ rf, data = fold_df, class = "Thr")$data
  lift_obj_data$fold = fold
  lift_df = rbind(lift_df, lift_obj_data)
}

lift_obj <- lift(binary ~ rf, data = for_lift, class = "Thr")

# Plot ROC
ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))


###Overlap
ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))+
  geom_line(data = lift_obj$data, aes(1 - Sp, Sn, color = liftModelVar), lwd=1, colour="black")

###Calculate accuracy within individual samples
ddply(model_15_d$pred, "Resample", summarise,
      accuracy = Accuracy(pred, obs))

###Youden index
#thresholder - This function uses the resampling results from a train object to generate performance statistics
#over a set of probability thresholds for two-class problems.
resample_stats_d <- thresholder(model_15_d,
                              threshold = seq(.1, 1, by = 0.0005),
                              final = TRUE,
                              statistics = "J")
ggplot(resample_stats_d, aes(x = prob_threshold_d, y = J)) +
  geom_point()

###Extract optimal probability threshold (that maximises Youden)
threshold_d <- resample_stats_d[which.max(resample_stats_d$J),]
threshold_d <- threshold_d$prob_threshold
threshold_d

####To finish
####Recalculate class (Thr/Non threatened) based on new threshold (0.86)
#pred_52 <- subset(model_15_2$pred, mtry==52)
#label <- ifelse(pred_52$Thr > threshold, 'Thr', 'nonThr')
#summary(label)
#levels(label)
#count(label == "Thr")
#count(label == "nonThr")

# ####To recalculate class (Thr/Non threatened) based on new threshold (0.762) using the final model
# pred <- predict(model_15_d, newdata = validationtSet_imp_d, type="prob")
# label_final <- as.factor(ifelse(pred$Thr > threshold, 'Thr', 'nonThr'))
# summary(label_final)
# ConfusionMatrix(label_final, validationtSet_imp_d$binary)
# 
# ###Compare with Lucie's threshold
# label_final_Bland <- as.factor(ifelse(pred$Thr > 0.614, 'Thr', 'nonThr'))
# summary(label_final_Bland)
# ConfusionMatrix(validationtSet_imp_d$binary,label_final_Bland)
# 
# ###To predict status of DD species
# pred_DD <- predict(model_15_d, newdata = predSet_imp_d, type="prob")
# label_DD <- as.factor(ifelse(pred_DD$Thr > threshold, 'Thr', 'nonThr'))
# summary(label_DD)
# 
# ###Compare with Lucie's threshold
# label_DD_Bland <- as.factor(ifelse(pred_DD$Thr > 0.614, 'Thr', 'nonThr'))
# summary(label_DD_Bland )

####Variable importance
varImp <- varImp(model_15_d, scale = FALSE)
varImp
plot(varImp, top=20)

ggplot(pred_DD, aes(x=Thr)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept=threshold_d, linetype="dotted", color="red")+
  geom_vline(xintercept=0.5, linetype="dotted", color="blue")









####NOT RUN~##########################################################################################  
#############Run the model#############################################
####savePredictions = "final"###########################################
###Exclude variables that are not used as predictors (binomial and threat status)
trainingSet_compl <- subset(trainingSet, select = -c(ID,Binomial,Genus, RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))

###Set model parameters
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = "final", 
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

# run a random forest model
start_time <- Sys.time()
model_15 <- train(binary ~ ., 
                  data = trainingSet_compl, 
                  method = "rf",
                  tuneLength=15, ####see https://rpubs.com/phamdinhkhanh/389752
                  ntree= 500,
                  trControl = model_control, 
                  metric="ROC")
end_time <- Sys.time()
end_time - start_time


####Print 
print(model_15)
plot(model_15)
print(model_15$finalModel) ####this gives the confusion matrix but doesn't calculate the accuracy 
summary(model_15$finalModel)
model_15$finalModel$confusion
model_15$finalModel$ntree
model_15$finalModel$mtry

###Average confusion matrix across 10 folds
confusionMatrix.train(model_15) #average across 10 samples #https://stats.stackexchange.com/questions/118568/random-forest-confusion-matrix

#Overall confusion matrix (only works if SavePred=FINAL in train)
#Predictions are in:
model_15$pred
#sorted as per CV fols, to sort as in original data frame:
model_15$pred[order(model_15$pred$rowIndex),2]
#to obtain a confusion matrix:
caret::confusionMatrix(model_15$pred[order(model_15$pred$rowIndex),2], as.factor(trainingSet_compl$binary))


###Confusion matrix on training data (it doesn't print it out, but doing the counts myself if gives perfect accuracy)
rfPred <- predict.train(model_15, trainingSet_compl, type = "raw")
confusionMatrix(rfPred, trainingSet_compl$binary) 
summary(as.factor(trainingSet_compl$binary))
summary(rfPred)


###Using evalm
###as per https://stackoverflow.com/questions/31138751/roc-curve-from-training-data-in-caret
## run MLeval
res <- evalm(model_15)
## get ROC
res$roc
## get calibration curve
res$cc
## get precision recall gain curve
res$prg


###Plot tree
source ("tree_func.R")
tree_num <- which(model_15$finalModel$forest$ndbigtree == max(model_15$finalModel$forest$ndbigtree))
tree_func(final_model = model_15$finalModel, tree_num)

#Warning messages:
#  1: Duplicated aesthetics after name standardisation: na.rm 
#2: Duplicated aesthetics after name standardisation: na.rm 
#3: Duplicated aesthetics after name standardisation: na.rm 
#4: Removed 98 rows containing missing values (geom_text_repel). 
#5: Removed 98 rows containing missing values (geom_label). 
#6: Removed 97 rows containing missing values (geom_label_repel). 


#####ROC curve for final average value
for_lift <- data.frame(binary = model_15$pred$obs, rf = model_15$pred$Thr)
lift_obj <- lift(binary ~ rf, data = for_lift, class = "Thr")   ###to check

# Plot ROC ----------------------------------------------------------------
ggplot(lift_obj$data) +
  geom_line(aes(1 - Sp, Sn, color = liftModelVar)) +
  scale_color_discrete(guide = guide_legend(title = "method"))#

#Get a ROC curve for each fold#############################################
for_lift <- data.frame(binary = model_15$pred$obs, rf = model_15$pred$Thr, resample = model_15$pred$Resample)
lift_df <-  data.frame()
for (fold in unique(for_lift$resample)) {
  fold_df <- dplyr::filter(for_lift, resample == fold)
  lift_obj_data <- lift(binary ~ rf, data = fold_df, class = "Thr")$data
  lift_obj_data$fold = fold
  lift_df = rbind(lift_df, lift_obj_data)
}

lift_obj <- lift(binary ~ rf, data = for_lift, class = "Thr")

# Plot ROC
ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))


###Overlap
ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))+
  geom_line(data = lift_obj$data, aes(1 - Sp, Sn, color = liftModelVar), lwd=1, colour="black")

###Calculate accuracy within individual samples
ddply(model_15$pred, "Resample", summarise,
      accuracy = Accuracy(pred, obs))


###Youden index
#thresholder - This function uses the resampling results from a train object to generate performance statistics
#over a set of probability thresholds for two-class problems.
resample_stats <- thresholder(model_15,
                              threshold = seq(.1, 1, by = 0.0005),
                              final = TRUE,
                              statistics = "J")
ggplot(resample_stats, aes(x = prob_threshold, y = J)) +
  geom_point()

###Extract optimal probability threshold (that maximises Youden)
threshold <- resample_stats[which.max(resample_stats$J),]
threshold <- threshold$prob_threshold
threshold

####Recalculate class (Thr/Non threatened) based on new threshold (0.86)
#pred_52 <- subset(model_15$pred, mtry==52)
#label <- ifelse(pred_52$Thr > threshold, 'Thr', 'nonThr')
#summary(label)
#levels(label)
#count(label == "Thr")
#count(label == "nonThr")


###Confusion matrix on validation data (it doesn't print it out, but doing the counts myself if gives perfect accuracy)
rfPred <- predict.train(model_15, newdata = validationtSet_compl, type = "raw")
confusionMatrix(rfPred, validationtSet_compl$binary) 
summary(as.factor(validationtSet_compl$binary))
summary(rfPred)


####To recalculate class (Thr/Non threatened) based on new threshold using the final model
pred <- predict.train(model_15, newdata = validationtSet_compl, type="prob")
label_final <- as.factor(ifelse(pred$Thr > threshold, 'Thr', 'nonThr'))
summary(as.factor(validationtSet_compl$binary))
summary(label_final)
ConfusionMatrix(label_final,validationtSet_compl$binary)



###Not sure why confusionMatrix doesn't work properly
#Predictions
#pred <- predict(model_15, newdata = validationtSet_compl,type="raw")
#table(pred)
#confusionMatrix(as.factor(pred),as.factor(validationtSet_compl$binary))

#levels(as.factor(validationtSet_compl$binary))
#levels(pred)

###Confusion matrix on training data (it doesn't print it out)
#rfPred <- predict.train(model_15, validationtSet_compl, type = "raw")
#confusionMatrix(rfPred, validationtSet_compl$binary) 
#summary(as.factor(validationtSet_compl$binary))
#summary(rfPred)

###Compare with Lucie's threshold
label_final_Bland <- as.factor(ifelse(pred$Thr > 0.614, 'Thr', 'nonThr'))
summary(label_final_Bland)
ConfusionMatrix(validationtSet_compl$binary,label_final_Bland)

###To predict status of DD species
pred_DD <- predict(model_15, newdata = predSet_complete, type="prob")
label_DD <- as.factor(ifelse(pred_DD$Thr > threshold, 'Thr', 'nonThr'))
summary(label_DD)

###Compare with Lucie's threshold
label_DD_Bland <- as.factor(ifelse(pred_DD$Thr > 0.614, 'Thr', 'nonThr'))
summary(label_DD_Bland )

####Variable importance
varImp <- varImp(model_15, scale = FALSE)
varImp
plot(varImp, top=20)

ggplot(pred_DD, aes(x=Thr)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept=threshold, linetype="dotted", color="red")+
  geom_vline(xintercept=0.5, linetype="dotted", color="blue")


