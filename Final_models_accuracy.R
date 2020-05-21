##########################################################################################################
###########predictions on the test set####################################################################
##########################################################################################################

####Complete dataset

###No Genus
####To recalculate class (Thr/Non threatened) based on new threshold using the final model
threshold_a <- 0.884
pred <- predict(model_15_a, newdata = validationtSet_compl_a, type="prob")
label_final <- as.factor(ifelse(pred$Thr > threshold_a, 'Thr', 'nonThr'))
summary(label_final)
ConfusionMatrix(label_final, validationtSet_compl_a$binary)

###To predict status of DD species
pred_DD <- predict(model_15_a, newdata = predSet_complete_a, type="prob")
label_DD <- as.factor(ifelse(pred_DD$Thr > threshold_a, 'Thr', 'nonThr'))
summary(label_DD)


############Includes genus
####To recalculate class (Thr/Non threatened) based on new threshold using the final model
threshold_b <- 0.886
pred <- predict(model_15_b, newdata = validationtSet_compl_b, type="prob")
label_final <- as.factor(ifelse(pred$Thr > threshold_b, 'Thr', 'nonThr'))
summary(label_final)
ConfusionMatrix(label_final, validationtSet_compl_b$binary)

###To predict status of DD species
pred_DD <- predict(model_15_b, newdata = predSet_complete_b, type="prob")
label_DD <- as.factor(ifelse(pred_DD$Thr > threshold_b, 'Thr', 'nonThr'))
summary(label_DD)


####Imputed dataset########################################################################################

###Include Genus
####To recalculate class (Thr/Non threatened) based on new threshold using the final model
threshold_c <- 0.83
pred <- predict(model_15_c, newdata = validationtSet_imp_c, type="prob")
label_final <- as.factor(ifelse(pred$Thr > threshold_c, 'Thr', 'nonThr'))
summary(label_final)
ConfusionMatrix(label_final, validationtSet_imp$binary)

###Imputed dataset
###To predict status of DD species
pred_DD <- predict(model_15_c, newdata = predSet_imp_c, type="prob")
label_DD <- as.factor(ifelse(pred_DD$Thr > threshold_c, 'Thr', 'nonThr'))
summary(label_DD)


###########################################################################
###No genus
####To recalculate class (Thr/Non threatened) based on new threshold using the final model
threshold_d <- 0.846
pred <- predict(model_15_d, newdata = validationtSet_imp_d, type="prob")
label_final <- as.factor(ifelse(pred$Thr > threshold_d, 'Thr', 'nonThr'))
summary(label_final)
ConfusionMatrix(label_final, validationtSet_imp_d$binary)


###To predict status of DD species
pred_DD <- predict(model_15_d, newdata = predSet_imp_d, type="prob")
label_DD <- as.factor(ifelse(pred_DD$Thr > threshold_d, 'Thr', 'nonThr'))
summary(label_DD)