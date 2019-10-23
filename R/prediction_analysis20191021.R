library(tidyverse)
library(plyr)
library(igraph)
library(dplyr)
setwd("data_analysis20191015/prediction/identification_table/")

#####################################################################################################
#####data preparation###############################################################################
#####################################################################################################

###This is batch data set
smartd_rplc_batch1 <- 
  readr::read_csv("smartd_rplc_batch1.csv")

smartd_rplc_batch2 <- 
  readr::read_csv("E:/project/smartD/smartD_batch1_2/RPLC/POS_NEG/ms1_data_rplc_batch2.csv")

###combine batch 1 and batch 2 data
##read the batch 1 and batch 2 match information
rplc_pos_batch1_batch2_matched_info <- 
  readr::read_csv("E:/project/smartD/smartD_batch1_2/RPLC/POS/rplc_pos_batch1_batch2_matched_info.csv")

rplc_pos_batch1_batch2_matched_info <- 
  readr::read_csv("E:/project/smartD/smartD_batch1_2/RPLC/POS/rplc_pos_batch1_batch2_matched_info.csv")

rplc_matched_info <- rbind(rplc_pos_batch1_batch2_matched_info[,-1],
                           rplc_neg_batch1_batch2_matched_info)


index1 <- 
  rplc_matched_info$batch1 %in%  smartd_rplc_batch1$name %>% 
  which()

index2 <- 
  rplc_matched_info$batch2 %in%  smartd_rplc_batch2$name %>% 
  which()

index <- intersect(index1, index2)

rplc_matched_info <- 
  rplc_matched_info[index,]


smartd_rplc_batch1 <- 
  smartd_rplc_batch1[match(rplc_matched_info$batch1, smartd_rplc_batch1$name),]

smartd_rplc_batch2 <- 
  smartd_rplc_batch2[match(rplc_matched_info$batch2, smartd_rplc_batch2$name),]

##remove QC and blank form batch 2 dataset
smartd_rplc_batch2 <-
  smartd_rplc_batch2 %>% 
  select(-contains("blk")) %>% 
  select(-contains("QC"))

sample_removed <-
  which(is.na(smartd_rplc_batch2), arr.ind = TRUE)[,2] %>% 
  unique() %>% 
  `[`(colnames(smartd_rplc_batch2), .)

##These samples may have positive or negative samples but don't have, so we should remove them
smartd_rplc_batch2 <- 
  smartd_rplc_batch2 %>% 
  dplyr::select(-sample_removed)

##remove the tags of metabolites
smartd_rplc_batch2 <- 
  smartd_rplc_batch2 %>% 
  dplyr::select(-(name:rt))

###rename P samples in batch2
colnames(smartd_rplc_batch2)[grep("P", colnames(smartd_rplc_batch2))] <-
  colnames(smartd_rplc_batch2)[grep("P", colnames(smartd_rplc_batch2))] %>% 
  stringr::str_replace("P", "") %>% 
  as.numeric() %>% 
  `+`(177) %>% 
  paste("X", ., sep = "")


##remove the tags of metabolites
sample_data1 <- 
  smartd_rplc_batch1 %>% 
  dplyr::select(., -(name:Database))

sample_data2 <- smartd_rplc_batch2

##remove NA samples in sample_data1 and sample_data2
##what samples have NAs
which(is.na(sample_data1), arr.ind = TRUE)[,2] %>% 
  unique()

colnames(sample_data1)[c(146, 147, 148)]

###"SFU65" "SFU74" "SFU91" may have no positive data, so remove them from the dataset
sample_data1 <- 
  sample_data1 %>% 
  dplyr::select(., -c(SFU65, SFU74, SFU91))

###metabolite tags
metabolite_tags <- 
  smartd_rplc_batch1 %>% 
  dplyr::select(., name:Database)

###batch 1 batch 2 data integration
mean1 <- apply(sample_data1, 1, mean)
mean2 <- apply(sample_data2, 1, mean)
ref <- mean1 / mean2

sample_data2 <- 
  sample_data2 * ref
####combine batch 1 and batch 2 data
sample_data <- 
  cbind(sample_data1, sample_data2)

###patient and sampel information
sample_info_191021 <- 
  readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")

match(colnames(sample_data), sample_info_191021$Sample_ID)
##P samples have no GA information

####log 10 and scale
sample_data <- 
  log(sample_data, 10) %>% 
  as.data.frame()

sample_data <-
  apply(sample_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

sample_data <- 
  sample_data %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_ID")

colnames(sample_data)[-1] <- 
  smartd_rplc_batch1$name

###add the patient information to the sample_data
sample_data <- 
  inner_join(x = sample_info_191021, 
             y = sample_data, 
             by = "Sample_ID")


###remnove samples with GA == 0
###sample_data_old is the data that contains PP samples
sampel_data_old <- 
  sample_data

sample_data <- 
  sample_data %>% 
  dplyr::filter(GA != 0)







#####################################################################################################
#####get discovery and validation data###############################################################
#####################################################################################################
###randomly divided into discovery and valdiation datasets
# set.seed(seed = 1)
sample_data$Patient_ID %>% unique()
##we have 36 peple in total
set.seed(1)
dis_patient_id <- 
  sample_data$Patient_ID %>% 
  unique() %>% 
  sample(size = 18, replace = FALSE)

dis_index <-
  (sample_data$Patient_ID %in% dis_patient_id) %>%
  which()

val_index <- setdiff(1:nrow(sample_data), dis_index)
# 
# intersect(dis_index, val_index)
# intersect(sample_data$Patient_ID[dis_index], 
#           sample_data$Patient_ID[val_index])
# load("dis_index")
# load("val_index")
# save(dis_index, file = "dis_index")
# save(val_index, file = "val_index")
sample_data_dis <- 
  sample_data[dis_index,]

sample_data_val <- 
  sample_data[val_index,]




################################################################################
###################################LASSO regression#############################
################################################################################
##using lasso regressio to select variables
setwd("lasso/")
library(glmnet)

###construct dataset for lasso regression
sample_data_dis_x <- 
  sample_data_dis %>% 
  dplyr::select(-(Patient_ID:`Birth Control at Discharge?`)) %>% 
  as.matrix()
sum(is.na(sample_data_dis_x))
sample_data_dis_x[,1]

sample_data_dis_y <- 
  sample_data_dis %>% 
  dplyr::select(GA) %>% 
  as.matrix()

lasso_regression2 <-
  cv.glmnet(
    x = sample_data_dis_x,
    y = sample_data_dis_y,
    family = "gaussian",
    type.measure = "mae",
    nfolds = 7,
    nlambda = 100,
    alpha = 1,
    standardize = FALSE
  )

plotLambdaVScoefficients(object = lasso_regression2$glmnet.fit)

plotLambdaVSerror(object = lasso_regression2)

##construct the best model
best_lasso <-
  glmnet(sample_data_dis_x, 
         sample_data_dis_y,
         family = "gaussian", 
         alpha = 1, 
         lambda = lasso_regression2$lambda.1se, 
         standardize = FALSE)

##markers
which(lasso_regression2$lambda == lasso_regression2$lambda.1se)

marker_lasso <- 
  which(lasso_regression2$lambda == lasso_regression2$lambda.1se) %>% 
  `[`(lasso_regression2$glmnet.fit$beta, ,.) %>% 
  # lasso_regression2$glmnet.fit$beta[,47] %>% 
  tibble(name = rownames(lasso_regression2$glmnet.fit$beta),
         coef = .) %>% 
  filter(coef != 0) 

##markers
marker_lasso %>% 
  left_join(., metabolite_tags, by = c("name")) %>% 
  arrange(coef) %>%
  # arrange(., desc(coef)) %>% 
  ggplot(aes(x = factor(Compound.name, Compound.name), y = coef)) +
  labs(x = "", y = "Coefficents") +
  geom_segment(aes(x = factor(Compound.name, Compound.name), 
                   xend = factor(Compound.name, Compound.name),
                   y = 0, yend = coef), 
               colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))

marker_lasso <- 
left_join(marker_lasso, metabolite_tags, by = "name")

write.csv(marker_lasso, "marker_lasso.csv", row.names = FALSE)
####
##validate in validation dataset
sample_data_val_x <- 
  sample_data_val %>% 
  dplyr::select(-c(Patient_ID:`Birth Control at Discharge?`)) %>% 
  as.matrix()

sum(is.na(sample_data_val_x))
sample_data_val_x[,1]

sample_data_val_y <- 
  sample_data_val %>% 
  dplyr::select(GA) %>% 
  as.matrix()

sample_data_val_x_lasso <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_lasso$name))

sample_data_val_x_lasso <- 
  apply(sample_data_val_x_lasso, 2, as.numeric)

sample_data_dis_x_lasso <- 
  sample_data_dis_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_lasso$name))

sample_data_dis_x_lasso <- 
  apply(sample_data_dis_x_lasso, 2, as.numeric)

###use validation dataset for validation
best_lasso2 <-
  glmnet(as.matrix(sample_data_dis_x_lasso), 
         sample_data_dis_y,
         family = "gaussian", 
         alpha = 1, 
         lambda = lasso_regression2$lambda.1se, 
         standardize = FALSE)

predicted_y <-
  predict(
    object = best_lasso2,
    newx = as.matrix(sample_data_val_x_lasso),
    s = lasso_regression2$lambda.1se
    # type = "response"
  )

prediction_self <- 
  predict(object = best_lasso2, 
          newx = as.matrix(sample_data_dis_x_lasso), 
          s = lasso_regression2$lambda.1se)

plot(sample_data_dis_y[,1], prediction_self)
abline(0, 1)

###from here we can see that why should us linear regression to correct prediction
data.frame(measured = sample_data_dis_y[,1],
           predicted = prediction_self[,1],
           stringsAsFactors = FALSE) %>% 
  mutate(diff = predicted - measured) %>% 
  mutate(colour = ifelse(diff > 0, "pos", "neg")) %>% 
  ggplot(aes(measured, diff)) +
  geom_segment(aes(x = measured, xend = measured, 
                   y = 0, yend = diff, colour = colour), show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2, aes(colour = colour), show.legend = FALSE) +
  geom_smooth(colour = "black", fill = "grey") +
  scale_colour_manual(values = c("pos" = "#800000FF", "neg" = "#155F83FF")) +
  theme_bw() +
  labs(x = "GA_week (measured)", y = "GA error (predicted - measured)") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))


linear_regression <- 
  lm(formula = sample_data_dis_y[,1] ~ prediction_self)

predicted_y2 <- 
  coef(linear_regression)[2] * predicted_y[,1] + coef(linear_regression)[1]

plot(sample_data_val_y[,1], predicted_y)
abline(0, 1)

plot(sample_data_val_y[,1], predicted_y2)
abline(0, 1)

abs(sample_data_val_y[,1] - predicted_y2) %>% 
  mean()

summary(lm(formula = predicted_y2~sample_data_val_y[,1]))


##bootstrap
predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_dis_x_lasso), 
           size = nrow(sample_data_dis_x_lasso), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x_lasso), dis_index)

  lasso_regression_best <- 
    glmnet(as.matrix(sample_data_dis_x_lasso[dis_index,]), 
           sample_data_dis_y[dis_index,],
           family = "gaussian", 
           alpha = 1, 
           lambda = lasso_regression2$lambda.1se)
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = lasso_regression_best,
      newx = as.matrix(sample_data_dis_x_lasso[dis_index,]),
      s = lasso_regression2$lambda.1se
    )[,1]
  
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = sample_data_dis_y[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = lasso_regression_best,
      newx = as.matrix(sample_data_dis_x_lasso[val_index,]),
      s = lasso_regression2$lambda.1se
    )[,1]
  
  temp_predict_y2 <- 
    coef(linear_regression)[2] * temp_predict_y + coef(linear_regression)[1]
  
  predict_y[[i]] <- 
    temp_predict_y2

  y[[i]] <- 
    sample_data_dis_y[val_index,1]
  
}

save(predicted_y, file = "predicted_y")
save(y, file = "y")

prediction <-
  data.frame(y = unlist(y),
             predict = unlist(predict_y),
             stringsAsFactors = FALSE)
temp <- 
  prediction %>% 
  dplyr::arrange(., y) %>% 
  dplyr::group_by(., y) %>% 
  dplyr::summarise(mean = mean(predict), sd = sd(predict))

plot(temp$y, temp$mean)


abline(0,1)

abs(temp$y - temp$mean) %>% 
  mean()

summary(lm(formula = temp$y~temp$mean))

temp %>% 
  mutate(ymax = mean + sd, ymin = mean - sd) %>% 
  ggplot(aes(x = y, y = mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  # scale_x_continuous(limits = c(10, 42)) +
  # scale_y_continuous(limits = c(10, 42)) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))









####predicted time to due date
setwd("data_analysis20191015/prediction/identification_table/lasso/time_to_due_prediction/")
sample_data_dis$Patient_ID

# isDate <- function(x){
#   stringr::str_detect(x, "[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}") %>% 
#     as_tibble() %>% 
#     filter(!is.na(value)) %>% 
#     filter(value) %>% 
#     nrow() > 10
# }

date_info_dis <- 
  sample_data_dis[,c("Patient_ID", "Sample_ID", "GA", "Visit","day_0", "Date.Acquired", "DD", "EDD")] %>% 
  mutate(begin_date = as.Date(day_0,"%m/%d/%Y"), 
         acquired_date = as.Date(Date.Acquired,"%m/%d/%y"),
         due_date = as.Date(DD, "%Y-%m-%d"),
         expected_due_date = as.Date(EDD, "%Y-%m-%d")
  ) %>% 
  select(-c(day_0, Date.Acquired, DD, EDD))


expected_date_remained_dis <-
  as.numeric(date_info_dis$due_date - date_info_dis$acquired_date, units = "weeks") %>% 
  as.matrix()


date_info_val <- 
  sample_data_val[,c("Patient_ID", "Sample_ID", "GA", "Visit","day_0", "Date.Acquired", "DD", "EDD")] %>% 
  mutate(begin_date = as.Date(day_0,"%m/%d/%Y"), 
         acquired_date = as.Date(Date.Acquired,"%m/%d/%y"),
         due_date = as.Date(DD, "%Y-%m-%d"),
         expected_due_date = as.Date(EDD, "%Y-%m-%d")
  ) %>% 
  select(-c(day_0, Date.Acquired, DD, EDD))

?difftime

expected_date_remained_val <-
  as.numeric(date_info_val$due_date - date_info_val$acquired_date, units = "weeks") %>% 
  as.matrix()


###construct dataset for lasso regression
lasso_regression2 <-
  cv.glmnet(
    x = sample_data_dis_x,
    y = expected_date_remained_dis,
    family = "gaussian",
    type.measure = "mae",
    nfolds = 7,
    nlambda = 100,
    alpha = 1,
    standardize = FALSE
  )

plotLambdaVScoefficients(object = lasso_regression2$glmnet.fit)

plotLambdaVSerror(object = lasso_regression2)

##construct the best model
best_lasso <-
  glmnet(sample_data_dis_x, 
         expected_date_remained_dis,
         family = "gaussian", 
         alpha = 1, 
         lambda = lasso_regression2$lambda.1se, 
         standardize = FALSE)

##markers
which(lasso_regression2$lambda == lasso_regression2$lambda.1se)

marker_lasso <- 
  which(lasso_regression2$lambda == lasso_regression2$lambda.1se) %>% 
  `[`(lasso_regression2$glmnet.fit$beta, ,.) %>% 
  # lasso_regression2$glmnet.fit$beta[,47] %>% 
  tibble(name = rownames(lasso_regression2$glmnet.fit$beta),
         coef = .) %>% 
  filter(coef != 0) 

##markers
marker_lasso %>% 
  left_join(., metabolite_tags, by = c("name")) %>% 
  arrange(coef) %>%
  # arrange(., desc(coef)) %>% 
  ggplot(aes(x = factor(Compound.name, Compound.name), y = coef)) +
  labs(x = "", y = "Coefficents") +
  geom_segment(aes(x = factor(Compound.name, Compound.name), 
                   xend = factor(Compound.name, Compound.name),
                   y = 0, yend = coef), 
               colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))

marker_lasso <- 
  left_join(marker_lasso, metabolite_tags, by = "name")

write.csv(marker_lasso, "marker_lasso.csv", row.names = FALSE)

##validate in validation dataset
sample_data_val_x <- 
  sample_data_val %>% 
  dplyr::select(-c(Patient_ID:`Birth Control at Discharge?`)) %>% 
  as.matrix()

sum(is.na(sample_data_val_x))
sample_data_val_x[,1]


sample_data_val_x_lasso <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_lasso$name))

sample_data_val_x_lasso <- 
  apply(sample_data_val_x_lasso, 2, as.numeric)

sample_data_dis_x_lasso <- 
  sample_data_dis_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_lasso$name))

sample_data_dis_x_lasso <- 
  apply(sample_data_dis_x_lasso, 2, as.numeric)

###use validation dataset for validation
best_lasso2 <-
  glmnet(as.matrix(sample_data_dis_x_lasso), 
         expected_date_remained_dis,
         family = "gaussian", 
         alpha = 1, 
         lambda = lasso_regression2$lambda.1se, 
         standardize = FALSE)

predicted_y <-
  predict(
    object = best_lasso2,
    newx = as.matrix(sample_data_val_x_lasso),
    s = lasso_regression2$lambda.1se
    # type = "response"
  )

prediction_self <- 
  predict(object = best_lasso2, 
          newx = as.matrix(sample_data_dis_x_lasso), 
          s = lasso_regression2$lambda.1se)

plot(expected_date_remained_dis[,1], prediction_self)
abline(0, 1)

###from here we can see that why should us linear regression to correct prediction
data.frame(measured = expected_date_remained_dis[,1],
           predicted = prediction_self[,1],
           stringsAsFactors = FALSE) %>% 
  mutate(diff = predicted - measured) %>% 
  mutate(colour = ifelse(diff > 0, "pos", "neg")) %>% 
  ggplot(aes(measured, diff)) +
  geom_segment(aes(x = measured, xend = measured, 
                   y = 0, yend = diff, colour = colour), show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2, aes(colour = colour), show.legend = FALSE) +
  geom_smooth(colour = "black", fill = "grey") +
  scale_colour_manual(values = c("pos" = "#800000FF", "neg" = "#155F83FF")) +
  theme_bw() +
  labs(x = "Time to due (measured)", y = "Predicted error (predicted - measured)") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

linear_regression <- 
  lm(formula = expected_date_remained_dis[,1] ~ prediction_self)

predicted_y2 <- 
  coef(linear_regression)[2] * predicted_y[,1] + coef(linear_regression)[1]

plot(expected_date_remained_val[,1], predicted_y)
abline(0, 1)

plot(expected_date_remained_val[,1], predicted_y2)
abline(0, 1)

abs(expected_date_remained_val[,1] - predicted_y2) %>% 
  mean()

summary(lm(formula = predicted_y2~expected_date_remained_val[,1]))

data.frame(
  measured = expected_date_remained_val[, 1],
  predicted = predicted_y2,
  stringsAsFactors = FALSE
) %>%
  ggplot(aes(measured, predicted)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2) +
  geom_point(size = 2, colour = "#FFA319FF") +
  # scale_x_continuous(limits = c(10, 42)) +
  # scale_y_continuous(limits = c(10, 42)) +
  labs(x = "Time to due (weeks, measured)", y = "Time to due (weeks, predicted)") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

##bootstrap
predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_dis_x_lasso), 
           size = nrow(sample_data_dis_x_lasso), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x_lasso), dis_index)
  
  lasso_regression_best <- 
    glmnet(as.matrix(sample_data_dis_x_lasso[dis_index,]), 
           expected_date_remained_dis[dis_index,],
           family = "gaussian", 
           alpha = 1, 
           lambda = lasso_regression2$lambda.1se)
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = lasso_regression_best,
      newx = as.matrix(sample_data_dis_x_lasso[dis_index,]),
      s = lasso_regression2$lambda.1se
    )[,1]
  
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = expected_date_remained_dis[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = lasso_regression_best,
      newx = as.matrix(sample_data_dis_x_lasso[val_index,]),
      s = lasso_regression2$lambda.1se
    )[,1]
  
  temp_predict_y2 <- 
    coef(linear_regression)[2] * temp_predict_y + coef(linear_regression)[1]
  
  predict_y[[i]] <- 
    temp_predict_y2
  
  y[[i]] <- 
    expected_date_remained_dis[val_index,1]
  
}

save(predicted_y, file = "predicted_y")
save(y, file = "y")

prediction <-
  data.frame(y = unlist(y),
             predict = unlist(predict_y),
             stringsAsFactors = FALSE)
temp <- 
  prediction %>% 
  dplyr::arrange(., y) %>% 
  dplyr::group_by(., y) %>% 
  dplyr::summarise(mean = mean(predict), sd = sd(predict))

plot(temp$y, temp$mean)


abline(0,1)

abs(temp$y - temp$mean) %>% 
  mean()

summary(lm(formula = temp$y~temp$mean))

temp %>% 
  mutate(ymax = mean + sd, ymin = mean - sd) %>% 
  ggplot(aes(x = y, y = mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  # scale_x_continuous(limits = c(10, 42)) +
  # scale_y_continuous(limits = c(10, 42)) +
  labs(x = "Time to due (weeks, measured)", 
       y = "Time to due (weeks, predicted)") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))


##############################################################################
####random forest
#############################################################################

####use other models
setwd("data_analysis20191015/prediction/identification_table/RF/")
library(randomForest)
##use boruta method
library(Boruta)
boruta_test <- 
  Boruta(x = sample_data_dis_x,
         y = sample_data_dis_y, 
         doTrace = 3, 
         holdHistory = TRUE)

plot(boruta_test)

boruta_test

marker_rf <- 
  boruta_test$finalDecision[boruta_test$finalDecision == "Confirmed"] %>% 
  names() %>% 
  sort()

marker_rf

no_infi <- function(x) all(!is.infinite(x))

boruta_test$ImpHistory[1:length(marker_rf),] %>% 
  t() %>% 
  as_tibble() %>% 
  select_if(., no_infi) %>% 
  dplyr::transmute(., mean = apply(., 1, mean),
                   sd = apply(., 1, sd),
                   ymax = mean + sd, 
                   ymin = mean - sd) %>% 
  mutate(name = colnames(boruta_test$ImpHistory)) %>% 
  filter(name %in% marker_rf) %>% 
  left_join(., metabolite_tags, by = "name") %>% 
  arrange(., mean) %>% 
  ggplot(.,aes(x = factor(Compound.name, Compound.name), y = mean)) +
  labs(x = "", y = "Importance") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF", width = 0) +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  coord_flip() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))

marker_rf <- 
  boruta_test$ImpHistory[1:length(marker_rf),] %>% 
  t() %>% 
  as_tibble() %>% 
  select_if(., no_infi) %>% 
  dplyr::transmute(., mean = apply(., 1, mean),
                   sd = apply(., 1, sd),
                   ymax = mean + sd, 
                   ymin = mean - sd) %>% 
  mutate(name = colnames(boruta_test$ImpHistory)) %>% 
  filter(name %in% marker_rf) %>% 
  select(name, everything())
colnames(marker_rf)[-1] <-
  colnames(marker_rf)[-1] %>% 
  paste(., "importance", sep = "_") 
  
marker_rf <- 
  marker_rf %>% 
left_join(metabolite_tags, by = "name")

write.csv(marker_rf, "marker_rf.csv", row.names = FALSE)

####parameter tunning
sample_data_dis_x_rf <- 
  sample_data_dis_x %>% 
  as.data.frame() %>% 
  select(., one_of(marker_rf$name)) %>% 
  as.matrix()

sample_data_dis_x_rf <- 
  apply(sample_data_dis_x_rf, 2, as.numeric)

sample_data_val_x_rf <- 
  sample_data_val_x %>% 
  as.data.frame() %>% 
  select(., one_of(marker_rf$name)) %>% 
  as.matrix()

sample_data_val_x_rf <- 
  apply(sample_data_val_x_rf, 2, as.numeric)

fgl.res <- tuneRF(sample_data_dis_x_rf, 
                  sample_data_dis_y[,1], 
                  mtryStart = 1,
                  stepFactor = 2, 
                  trace = TRUE, 
                  plot = TRUE)


fgl.res %>% 
  as.data.frame() %>% 
  ggplot(aes(mtry, y = OOBError)) +
  geom_line(colour = "#155F83FF") +
  geom_point(size = 4, colour = "#FFA319FF") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))


rf_regression <-
  randomForest(x = sample_data_dis_x_rf, 
               y = sample_data_dis_y[,1], 
               replace = TRUE, 
               importance = TRUE,
               proximity = TRUE, 
               mtry = 4)

##validate in validation dataset
###construct dataset for lasso regression
sample_data_val_x_rf <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_rf$name))

###use validation dataset for validation
predicted_y <-
predict(
  object = rf_regression,
  newdata = sample_data_val_x_rf
  # type = "response"
)

plot(
  sample_data_val_y[,1],
  predicted_y
)

abline(0, 1)

prediction_self <- 
  predict(object = rf_regression, 
          newx = as.matrix(sample_data_dis_x_rf))

plot(sample_data_dis_y[,1], prediction_self)
abline(0, 1)


###from here we can see that why should us linear regression to correct prediction
data.frame(measured = sample_data_dis_y[,1],
           predicted = prediction_self,
           stringsAsFactors = FALSE) %>% 
  mutate(diff = predicted - measured) %>% 
  mutate(colour = ifelse(diff > 0, "pos", "neg")) %>% 
  ggplot(aes(measured, diff)) +
  geom_segment(aes(x = measured, xend = measured, 
                   y = 0, yend = diff, colour = colour), show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2, aes(colour = colour), show.legend = FALSE) +
  geom_smooth(colour = "black", fill = "grey") +
  scale_colour_manual(values = c("pos" = "#800000FF", "neg" = "#155F83FF")) +
  theme_bw() +
  labs(x = "GA_week (measured)", y = "GA_week (predicted - measured)") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))


linear_regression <- 
  lm(formula = sample_data_dis_y[,1] ~ prediction_self)

predicted_y2 <- 
  coef(linear_regression)[2] * predicted_y + coef(linear_regression)[1]


data.frame("measured" = sample_data_val_y[,1], 
           'predicted' = predicted_y2,
           stringsAsFactors = FALSE) %>%  
  ggplot(aes(x = measured, predicted)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  scale_x_continuous(limits = c(10, 42)) +
  scale_y_continuous(limits = c(10, 42)) +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        plot.title = element_text(size = 15))
  


 abs(sample_data_val_y[,1] - predicted_y2) %>% 
  mean()

summary(lm(formula = predicted_y2~sample_data_val_y[,1]))


##bootstrap
predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_dis_x_rf), 
           size = nrow(sample_data_dis_x_rf), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x_rf), dis_index)
  
  rf_regression_temp <-
    randomForest(x = sample_data_dis_x_rf[dis_index,], 
                 y = sample_data_dis_y[dis_index,1], 
                 replace = TRUE, 
                 importance = TRUE,
                 proximity = TRUE,
                 mtry = 4)
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = rf_regression_temp,
      newx = as.matrix(sample_data_dis_x_rf[dis_index,])
    )
  
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = sample_data_dis_y[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = rf_regression_temp,
      newdata = as.matrix(sample_data_dis_x_rf[val_index,])
    )
  
  temp_predict_y2 <- 
    coef(linear_regression)[2] * temp_predict_y + coef(linear_regression)[1]
  
  predict_y[[i]] <- 
    temp_predict_y2
  
  y[[i]] <- 
    sample_data_dis_y[val_index,1]
  
}

save(predict_y, file = "predict_y")
save(y, file = "y")

prediction <-
  data.frame(y = unlist(y),
             predict = unlist(predict_y),
             stringsAsFactors = FALSE)

temp <- 
  prediction %>% 
  arrange(., y) %>% 
  group_by(., y) %>% 
  dplyr::summarise(mean = mean(predict), sd = sd(predict))


abs(temp$y - temp$mean) %>% 
  mean()

summary(lm(formula = temp$y~temp$mean))

temp %>% 
  mutate(ymax = mean + sd, ymin = mean - sd) %>% 
  ggplot(aes(x = y, y = mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  scale_x_continuous(limits = c(10, 42)) +
  scale_y_continuous(limits = c(10, 42)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))


###############################################################################
######for other prediction (time to due)
####use Random Forest
library(randomForest)
##use boruta method
library(Boruta)
boruta_test <- 
  Boruta(x = sample_data_dis_x,
         y = expected_date_remained_dis, 
         doTrace = 3, 
         holdHistory = TRUE)

plot(boruta_test)

boruta_test

marker_rf <- 
  boruta_test$finalDecision[boruta_test$finalDecision == "Confirmed"] %>% 
  names() %>% 
  sort()

marker_rf

no_infi <- function(x) all(!is.infinite(x))

boruta_test$ImpHistory[1:length(marker_rf),] %>% 
  t() %>% 
  as_tibble() %>% 
  select_if(., no_infi) %>% 
  dplyr::transmute(., mean = apply(., 1, mean),
                   sd = apply(., 1, sd),
                   ymax = mean + sd, 
                   ymin = mean - sd) %>% 
  mutate(name = colnames(boruta_test$ImpHistory)) %>% 
  filter(name %in% marker_rf) %>% 
  left_join(., metabolite_tags, by = "name") %>% 
  arrange(., mean) %>% 
  ggplot(.,aes(x = factor(Compound.name, Compound.name), y = mean)) +
  labs(x = "", y = "Importance") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF", width = 0) +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  coord_flip() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))

marker_rf <- 
  boruta_test$ImpHistory[1:length(marker_rf),] %>% 
  t() %>% 
  as_tibble() %>% 
  select_if(., no_infi) %>% 
  dplyr::transmute(., mean = apply(., 1, mean),
                   sd = apply(., 1, sd),
                   ymax = mean + sd, 
                   ymin = mean - sd) %>% 
  mutate(name = colnames(boruta_test$ImpHistory)) %>% 
  filter(name %in% marker_rf) %>% 
  select(name, everything())
colnames(marker_rf)[-1] <-
  colnames(marker_rf)[-1] %>% 
  paste(., "importance", sep = "_") 

marker_rf <- 
  marker_rf %>% 
  left_join(metabolite_tags, by = "name")

write.csv(marker_rf, "marker_rf.csv", row.names = FALSE)

####parameter tunning
sample_data_dis_x_rf <- 
  sample_data_dis_x %>% 
  as.data.frame() %>% 
  select(., one_of(marker_rf$name)) %>% 
  as.matrix()

sample_data_dis_x_rf <- 
  apply(sample_data_dis_x_rf, 2, as.numeric)

sample_data_val_x_rf <- 
  sample_data_val_x %>% 
  as.data.frame() %>% 
  select(., one_of(marker_rf$name)) %>% 
  as.matrix()

sample_data_val_x_rf <- 
  apply(sample_data_val_x_rf, 2, as.numeric)

fgl.res <- tuneRF(sample_data_dis_x_rf, 
                  expected_date_remained_dis[,1], 
                  mtryStart = 1,
                  stepFactor = 2, 
                  trace = TRUE, 
                  plot = TRUE)


fgl.res %>% 
  as.data.frame() %>% 
  ggplot(aes(mtry, y = OOBError)) +
  geom_line(colour = "#155F83FF") +
  geom_point(size = 4, colour = "#FFA319FF") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))


rf_regression <-
  randomForest(x = sample_data_dis_x_rf, 
               y = expected_date_remained_dis[,1], 
               replace = TRUE, 
               importance = TRUE,
               proximity = TRUE, 
               mtry = 4)


##validate in validation dataset
###construct dataset for lasso regression
sample_data_val_x_rf <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_rf))

###use validation dataset for validation
predicted_y <-
  predict(
    object = rf_regression,
    newdata = sample_data_val_x_rf
    # type = "response"
  )


prediction_self <- 
  predict(object = rf_regression, 
          newx = as.matrix(sample_data_dis_x_rf))

plot(expected_date_remained_dis[,1], prediction_self)
abline(0, 1)


###from here we can see that why should us linear regression to correct prediction
data.frame(measured = expected_date_remained_dis[,1],
           predicted = prediction_self,
           stringsAsFactors = FALSE) %>% 
  mutate(diff = predicted - measured) %>% 
  mutate(colour = ifelse(diff > 0, "pos", "neg")) %>% 
  ggplot(aes(measured, diff)) +
  geom_segment(aes(x = measured, xend = measured, 
                   y = 0, yend = diff, colour = colour), show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2, aes(colour = colour), show.legend = FALSE) +
  geom_smooth(colour = "black", fill = "grey") +
  scale_colour_manual(values = c("pos" = "#800000FF", "neg" = "#155F83FF")) +
  theme_bw() +
  labs(x = "Time to due (measured)", y = "Predicted error (predicted - measured)") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))


linear_regression <- 
  lm(formula = expected_date_remained_dis[,1] ~ prediction_self)

predicted_y2 <- 
  coef(linear_regression)[2] * predicted_y + coef(linear_regression)[1]

plot(expected_date_remained_val[,1], predicted_y)
abline(0, 1)

plot(expected_date_remained_val[,1], predicted_y2)
abline(0, 1)

abs(expected_date_remained_val[,1] - predicted_y2) %>% 
  mean()

summary(lm(formula = predicted_y2~expected_date_remained_val[,1]))

##bootstrap
predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_dis_x_rf), 
           size = nrow(sample_data_dis_x_rf), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x_rf), dis_index)
  
  rf_regression_temp <-
    randomForest(x = sample_data_dis_x_rf[dis_index,], 
                 y = expected_date_remained_dis[dis_index,1], 
                 replace = TRUE, 
                 importance = TRUE,
                 proximity = TRUE,
                 mtry = 4)
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = rf_regression_temp,
      newx = as.matrix(sample_data_dis_x_rf[dis_index,])
    )
  
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = expected_date_remained_dis[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = rf_regression_temp,
      newdata = as.matrix(sample_data_dis_x_rf[val_index,])
    )
  
  temp_predict_y2 <- 
    coef(linear_regression)[2] * temp_predict_y + coef(linear_regression)[1]
  
  predict_y[[i]] <- 
    temp_predict_y2
  
  y[[i]] <- 
    expected_date_remained_dis[val_index,1]
  
}

save(predict_y, file = "predict_y")
save(y, file = "y")


prediction <-
  data.frame(y = unlist(y),
             predict = unlist(predict_y),
             stringsAsFactors = FALSE)

temp <- 
  prediction %>% 
  arrange(., y) %>% 
  group_by(., y) %>% 
  dplyr::summarise(mean = mean(predict), sd = sd(predict))

plot(temp$y, temp$mean)


abline(0,1)

abs(temp$y - temp$mean) %>% 
  mean()

summary(lm(formula = temp$y~temp$mean))

temp %>% 
  mutate(ymax = mean + sd, ymin = mean - sd) %>% 
  ggplot(aes(x = y, y = mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  # scale_x_continuous(limits = c(10, 42)) +
  # scale_y_continuous(limits = c(10, 42)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))


























