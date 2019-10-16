
###RPLC data
###align batch 1 and batch 2 data
##positive mode
setwd("smartD_batch1_2/RPLC/POS/")
load("smartd_rplc_pos_batch1_5")
load("smartd_rplc_pos_batch2_5")

library(tidyverse)

batch1_pos_info <- smartd_rplc_pos_batch1_5@ms1.data[[1]] %>% 
  select(., name, mz, rt, contains("QC"))

batch2_pos_info <- smartd_rplc_pos_batch2_5@ms1.data[[1]] %>% 
  select(., name, mz, rt, contains("QC"))

data1 <- 
  batch1_pos_info %>% 
  select(., mz:rt)

data2 <- 
  batch2_pos_info %>% 
  select(., mz:rt)

match_result <-
  sxtTools::sxtMTmatch(data1 = as.matrix(data1), 
                     data2 = as.matrix(data2), 
                     mz.tol = 25, rt.tol = 30, 
                     rt.error.type = "abs")


unique(match_result[,1])

name1 <- batch1_pos_info$name[match_result[,1]]
name2 <- batch2_pos_info$name[match_result[,2]]

rplc_pos_batch1_batch2_matched_info <- 
  data.frame(batch1 = name1, batch2 = name2, 
             stringsAsFactors = FALSE)


write.csv(rplc_pos_batch1_batch2_matched_info, 
          'rplc_pos_batch1_batch2_matched_info.csv')



##negative mode
setwd("smartD_batch1_2/RPLC/NEG/")
load("smartd_rplc_neg_batch1_5")
load("smartd_rplc_neg_batch2_5")

library(tidyverse)

batch1_neg_info <- smartd_rplc_neg_batch1_5@ms1.data[[1]] %>% 
  select(., name, mz, rt, contains("QC"))

batch2_neg_info <- smartd_rplc_neg_batch2_5@ms1.data[[1]] %>% 
  select(., name, mz, rt, contains("QC"))

data1 <- 
  batch1_neg_info %>% 
  select(., mz:rt)

data2 <- 
  batch2_neg_info %>% 
  select(., mz:rt)

match_result <-
  sxtTools::sxtMTmatch(data1 = as.matrix(data1), 
                       data2 = as.matrix(data2), 
                       mz.tol = 25, rt.tol = 30, 
                       rt.error.type = "abs")


unique(match_result[,1])

name1 <- batch1_neg_info$name[match_result[,1]]
name2 <- batch2_neg_info$name[match_result[,2]]

rplc_neg_batch1_batch2_matched_info <- 
  data.frame(batch1 = name1, batch2 = name2, 
             stringsAsFactors = FALSE)


write.csv(rplc_neg_batch1_batch2_matched_info, 
          'rplc_neg_batch1_batch2_matched_info.csv')


####how many metabolites can be found in batch2

setwd("smartD_batch1_2/RPLC/POS_NEG/")
smartd_rplc_batch1 <- readr::read_csv("smartd_rplc_batch1.csv")

name_pos <- paste(rplc_pos_batch1_batch2_matched_info$batch1, "POS", sep = "_")
name_neg <- paste(rplc_neg_batch1_batch2_matched_info$batch1, "NEG", sep = "_")

sum(smartd_rplc_batch1$name %in% c(name_pos, name_neg))













