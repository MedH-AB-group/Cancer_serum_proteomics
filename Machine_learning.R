## ----libraries, message=F, cache = FALSE, include=FALSE-----------------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(knitr)
library(stringr)
library(doParallel) 
library(glmnet)
library(randomForest)
library(caret)
library(ggvenn)
library(skimr)
library(pROC)




## ----read_data----------------------------------------------------------------------------------------------------------------------------
###getting tabular data and scaling the data
data <- read.csv2("data/protein_screening.csv")
data <- data %>%
        mutate_if(is.character, as.numeric) %>%
        mutate_at("sample_id", as.character) %>%
        arrange(group) 
scaled_data <- cbind(data[1:2], as.data.frame(scale(data[3:length(data)])))


## ----split_data---------------------------------------------------------------------------------------------------------------------------
###splitting the data 80 to 20%, controlling it and defining all necessary objects
scaled_data$index <- as.numeric(rownames(scaled_data))
inTrain <- createDataPartition(y = scaled_data$index, p=0.8, list = FALSE)
training <- scaled_data[inTrain,]
nrow(training)
testing  <- scaled_data[-inTrain,]
nrow(testing)
# skimmed <- skim(training) (takes long with many vars)
outcome_test <- testing %>%
           mutate(group = ifelse(group == '4', T, F)) %>%
           dplyr::select(group)
testing <- as.matrix(testing %>% select(-group, -sample_id))
## Convert data to matrix format (required by glmnet)
ptn_matrix <- as.matrix(training %>% select(-group, -sample_id))
outcome <- training %>%
           mutate(group = ifelse(group == '4', T, F)) %>%
           dplyr::select(group)


## ----optimum_alpha_lambda-----------------------------------------------------------------------------------------------------------------
###Regularization, finding best lambda for each value of alpha and controlling behaviour of alpha
Alpha <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
# # Find optimal lambda for each alpha
# for (j in 1:length(Alpha)) {
#     cv <- cv.glmnet(ptn_matrix, outcome$group, alpha = Alpha[j], nfolds = 10)
#     # Access the optimal lambda (cv$lambda.min) for each alpha
#     # You can also get other performance metrics from cv$cvm
#     print(paste("Optimal lambda for alpha =", Alpha[j], ":", cv$lambda.min))
# }
# Find best alpha based on AIC and BIC values
aic_values <- numeric(length(Alpha))
bic_values <- numeric(length(Alpha))
# Perform cross-validation for each alpha
for (i in seq_along(Alpha)) {
    cv_fit <- cv.glmnet(ptn_matrix, outcome$group, alpha = Alpha[i], nfolds = 10)
    aic_values[i] <- cv_fit$cvm[which.min(cv_fit$cvm)] + 2 * cv_fit$glmnet.fit$df
    bic_values[i] <- cv_fit$cvm[which.min(cv_fit$cvm)] + log(length(outcome$group)) * cv_fit$glmnet.fit$df
}
# Find the optimal alpha based on AIC and BIC
optimal_alpha_aic <- Alpha[which.min(aic_values)]
optimal_alpha_bic <- Alpha[which.min(bic_values)]
cat("Optimal alpha (AIC):", optimal_alpha_aic, "\n")
cat("Optimal alpha (BIC):", optimal_alpha_bic, "\n")


## ----lasso--------------------------------------------------------------------------------------------------------------------------------
# Fit LASSO model using glmnet
###Use same chunk with a different alpha value for elastic net (alpha should be well above 0 to avoid falling in ridge regression)
lasso_model <- glmnet(x = ptn_matrix, y = outcome$group, alpha = 1)
# Cross-validation to select the optimal lambda (regularization parameter)
lasso_result <- cv.glmnet(x = ptn_matrix, y = outcome$group, alpha = 1)
# Get the optimal lambda from cross-validation
optimal_lambda <- lasso_result$lambda.min
# Extract the coefficients for the selected features (non-zero coefficients)
lasso_coefficients <- coef(lasso_model, s = optimal_lambda)
str(lasso_coefficients)
lasso_coefficients <- as.matrix(lasso_coefficients)
lasso_coefficients = as.data.frame(lasso_coefficients) %>%
     filter(s1 != 0) %>%
     rownames_to_column("probes") %>%
     filter(!probes %in% c("index", "(Intercept)"))
# save(lasso_coefficients, file= "results/probes_lasso_608ptn.RData")
# save(lasso_coefficients, file = "results/probes_lasso.RData")
get(load("results/probes_lasso.RData"))


## ----test_lasso---------------------------------------------------------------------------------------------------------------------------
###Testing the model
predictions <- predict(lasso_model, newx = testing, s = 1)
mse <- mean((outcome_test$group - predictions)^2)
rmse <- sqrt(mse)
r_squared <- 1 - mse / var(outcome_test$group)
mse
rmse
r_squared


## ----visualize----------------------------------------------------------------------------------------------------------------------------
###Visualize resutls
plot_data <- data.frame(Actual = outcome_test$group, Predicted = predictions)
# pdf(file = "results/Paper_illustratios/Lasso_test_res2.pdf")
ggplot(plot_data, aes(x = Actual, y = s1)) +
    geom_point() +  # Scatter points
    geom_smooth(method = "lm", color = "red", se = FALSE) +  # Regression line
    labs(title = "LASSO model test results",
         x = "Actual",
         y = "Predicted")
 # while (!is.null(dev.list()))  dev.off()
# pdf(file = "results/Paper_illustratios/Lasso_test_res1.pdf")
plot(outcome_test$group, predictions, main = "LASSO model test results", xlab = "Actual", ylab = "Predicted") +
abline(a = 0, b = 1, col = "red")  # Add a diagonal reference line
# while (!is.null(dev.list()))  dev.off()


## ----lasso_roc----------------------------------------------------------------------------------------------------------------------------
##ROC curve
roc_obj <- roc(outcome_test$group, predictions)
plot(roc_obj, main = "ROC Curve", print.auc = TRUE, legacy.axes = TRUE)




## ----elastic_net--------------------------------------------------------------------------------------------------------------------------
# Fit an Elastic Net model
control <- trainControl(method="loocv")
enet_fit <- train(ptn_matrix[,2:ncol(ptn_matrix)], ptn_matrix[,1], method = 'glmnet',  trControl = control)  ##method = "lm",
enet_fit$results$lambda <- round(enet_fit$results$lambda, digits = 2) 
# pdf(file = "results/Paper_illustratios/ELNET_lambda_cv.pdf")
plot(enet_fit, xvar = "lambda", label = TRUE)
# while (!is.null(dev.list()))  dev.off()
# Cross-validate Elastic Net
cv_fit <- cv.glmnet(x = ptn_matrix, y = outcome$group, alpha = 0.1)
# pdf(file = "results/Paper_illustratios/ELNET_cv_fit.pdf")
plot(cv_fit)
# while (!is.null(dev.list()))  dev.off()


## ----test_EN------------------------------------------------------------------------------------------------------------------------------
##Testing the EN model
predictions <- predict(elasticnet_model, newx = testing, s = 0.05)
mse <- mean((outcome_test$group - predictions)^2)
rmse <- sqrt(mse)
r_squared <- 1 - mse / var(outcome_test$group)
mse
rmse
r_squared


## ----visualize2---------------------------------------------------------------------------------------------------------------------------
###Results visualization
plot_data <- data.frame(Actual = outcome_test$group, Predicted = predictions)
# pdf(file = "results/Paper_illustratios/ELNET_test_res.pdf")
ggplot(plot_data, aes(x = Actual, y = s1)) +
    geom_point() +  # Scatter points
    geom_smooth(method = "lm", color = "red", se = FALSE) +  # Regression line
    labs(title = "Elastic Net model test results",
         x = "Actual",
         y = "Predicted")
 # while (!is.null(dev.list()))  dev.off()
# pdf(file = "results/Paper_illustratios/ELNET_test_res1.pdf")
plot(outcome_test$group, predictions, main = "LASSO model test results", xlab = "Actual", ylab = "Predicted") +
abline(a = 0, b = 1, col = "red")  # Add a diagonal reference line
# while (!is.null(dev.list()))  dev.off()


## ----sig_feat-----------------------------------------------------------------------------------------------------------------------------
####saving the resutls
sig.feat <- rfe_features %>%
                        rownames_to_column("probes") %>%
                        # dplyr::select(probes) %>%
                        mutate(EN = T) %>%
            full_join(lasso_coefficients %>%
                        # select(probes) %>%
                        mutate(Lasso = T),
                      by = "probes") %>%
            replace(is.na(.), FALSE) 
# save(sig.feat, file = "results/GBC_benign.rfe_lasso_sig.RData")


## ----sig_feat_visualization---------------------------------------------------------------------------------------------------------------
###Visualize intersection of the 2 methods
as_tibble(sig.feat)
ggvenn(sig.feat, c("EN", "Lasso"), 
  fill_color = c("#a141a5", "#00FF00"),
  stroke_size = 0.5, set_name_size = 4)


## ----feature_visualization----------------------------------------------------------------------------------------------------------------
###Plot boxes and densitites of the chosen proteins by ML
sig.feat <- sig.feat %>%
           filter(Lasso == T)
featurePlot(x = scaled_data[,c("group", sig.feat$probes)], 
            y = as.factor(scaled_data$group), 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))
featurePlot(x = scaled_data[sig.feat$probes], 
            y = as.factor(scaled_data$group), 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

