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

####################################################################################################
####Dataset preparation#############################################################################
####################################################################################################

#####Separate out different types of variables based on the transformation that we are going to do
###
BM_other <- BM[,c("ID", "Binomial", "Order", "Family", "Genus")]
BM_binary <- BM[,c("ID", "Afrotropical", "Australasia", "Indo_malayan", "Nearctic", "Neotropical", "Oceania", "Palearctic")]
BM_cat <- BM[,c("ID", "RedList","ReproductiveMode", "TrophicGroup", "HabitatMode", "Continent")]
BM_numeric <- BM[,c("ID", "BodySize", "HabitatsIUCN", "EOO", "Latitude", "ElevMin", "Prec", "PrecSeas", "Temp", "TempSeas", "HPD", "HPDMin", "HumanFootprint", "Accessibility")]

###Categorical variables transformed to orthogonal dummy variables##############
###Create the dummy variables
dummies <- dummyVars(" ~ .", data = BM_cat, fullRank=T)
head(dummies)
BM_cat_new <- data.frame(predict(dummies, newdata = BM_cat))
summary(BM_cat_new)


####Find out which variables have near zero variance in each of these subsets (categorical variables already transformed)
####Based on the settings from Bland and Bohm no variables are excluded

nzv_other <- nearZeroVar(BM_other, freqCut = 999/1, uniqueCut = 0.01, saveMetrics = FALSE)
nzv_other

nzv_binary <- nearZeroVar(BM_binary, freqCut = 999/1, uniqueCut = 0.01, saveMetrics = FALSE)
nzv_binary

nzv_cat <- nearZeroVar(BM_cat_new , freqCut = 999/1, uniqueCut = 0.01, saveMetrics = FALSE)
nzv_cat

nzv_numeric <- nearZeroVar(BM_numeric , freqCut = 999/1, uniqueCut = 0.01, saveMetrics = FALSE)
nzv_numeric


#Using caret's default setting s(freqCut = 19 and uniqueCut = 10, which is fairly conservative)
nzv_other <- nearZeroVar(BM_other, freqCut = 19, uniqueCut = 10, saveMetrics = FALSE)
nzv_other 
head(BM_other[, nzv_other]) ####variable 3 (Order) would get excluded but it's a predictor in the final model it can't have been excluded

nzv_binary <- nearZeroVar(BM_binary, freqCut = 19, uniqueCut = 10, saveMetrics = FALSE)
nzv_binary 
head(BM_binary[, nzv_binary]) ###Oceania is excluded

nzv_cat <- nearZeroVar(BM_cat_new , freqCut = 19, uniqueCut = 10, saveMetrics = FALSE)
nzv_cat 
head(BM_cat_new[, nzv_cat])###excludes ReproductiveMode.Ovoviviparous, TrophicGroup.Herbivorous, TrophicGroup.Omnivorous, HabitatMode.Semi.aquatic, HabitatMode.Semi.arboreal, HabitatMode.Semi.fossorial, HabitatMode.Semi.saxicolous, Continent.Continent.Island

nzv_numeric <- nearZeroVar(BM_numeric , freqCut = 19, uniqueCut = 10, saveMetrics = FALSE)
nzv_numeric ###excludes none

###exclude them from dataset (only if chosing the stricter settings above)
BM_binary_new <- BM_binary[, -nzv_binary]
BM_cat_new <- BM_cat_new[, -nzv_cat]

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
BM_final <- merge(BM_binary_new,BM_cat_new, by = "ID")
###Look at the data
dim(BM_final) # check dataset

BM_final <- merge(BM_final,BM_numeric_new, by = "ID")
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

# Setting nonassessed species aside
data_def_complete <- data.frame(subset (BM_final_complete, RedList.DD==1))
data_def_complete$binary<-NA
head(data_def_complete)
dim(data_def_complete)

###########################################################################################

####Impute missing values for the model with imputed values
BM_imp <- missForest(BM_final, verbose = TRUE)
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


#######################################################################################
#############Run the model#############################################
#######################################################################################

###Exclude variables that are not used as predictors (binomial and threat status)
data_suff_compl_a <- subset(data_suff_complete, select = -c(ID,Binomial,RedList.DD, RedList.EN, RedList.LC,RedList.NT,RedList.VU))

###Set model parameters
model_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  savePredictions = TRUE, 
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

# run a random forest model
start_time <- Sys.time()
model_15 <- train(binary ~ ., 
                  data = data_suff_compl_a, 
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


###Using evalm
###as per https://stackoverflow.com/questions/31138751/roc-curve-from-training-data-in-caret
## run MLeval
res <- evalm(model_15)
## get ROC
res_c$roc
## get calibration curve
res_c$cc
## get precision recall gain curve
res_c$prg


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

####To finish
####Recalculate class (Thr/Non threatened) based on new threshold (0.86)
#pred_52 <- subset(model_15$pred, mtry==52)
#label <- ifelse(pred_52$Thr > threshold, 'Thr', 'nonThr')
#summary(label)
#levels(label)
#count(label == "Thr")
#count(label == "nonThr")

####To recalculate class (Thr/Non threatened) based on new threshold (0.86) using the final model
pred <- predict(model_15, newdata = data_suff_complete, type="prob")
label_final <- as.factor(ifelse(pred$Thr > threshold, 'Thr', 'nonThr'))
summary(label_final)
ConfusionMatrix(label_final, model_15$final$y)
confusionMatrix(label_final, model_15$final$y)$overall


###Compare with Lucie's threshold
label_final_Bland <- as.factor(ifelse(pred$Thr > 0.6, 'Thr', 'nonThr'))
summary(label_final_Bland)
ConfusionMatrix(model_15_c$final$y,label_final_Bland)

###To predict status of DD species
pred_DD <- predict(model_15, newdata = data_def_complete, type="prob")
label_DD <- as.factor(ifelse(pred_DD$Thr > threshold, 'Thr', 'nonThr'))
summary(label_DD)

###Compare with Lucie's threshold
label_DD_Bland <- as.factor(ifelse(pred_DD$Thr > 0.6, 'Thr', 'nonThr'))
summary(label_DD_Bland )

####Variable importance
varImp <- varImp(model_15, scale = FALSE)
varImp
plot(varImp, top=20)