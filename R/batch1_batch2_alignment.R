
###RPLC data
###align batch 1 and batch 2 data
##positive mode
setwd("smartD_batch1_2/RPLC/POS/")
load("smartd_rplc_pos_batch1_5")
load("smartd_rplc_pos_batch2_5")

library(tidyverse)

batch1_pos_info <- smartd_rplc_pos_batch1_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

batch2_pos_info <- smartd_rplc_pos_batch2_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

data1 <- 
  batch1_pos_info %>% 
  dplyr::select(., mz:rt)

data2 <- 
  batch2_pos_info %>% 
  dplyr::select(., mz:rt)

match_result <-
  sxtTools::sxtMTmatch(data1 = as.matrix(data1), 
                     data2 = as.matrix(data2), 
                     mz.tol = 25, rt.tol = 120, 
                     rt.error.type = "abs")

###remove the 1:multiple
match_result <- 
  match_result %>% 
  as_tibble() %>% 
  group_by(Index1) %>% 
  filter(`rt error` == min(`rt error`)) %>% 
  ungroup()


unique(match_result[,1])

name1 <- batch1_pos_info$name[match_result$Index1]
name2 <- batch2_pos_info$name[match_result$Index2]

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
  dplyr::select(., name, mz, rt, contains("QC"))

batch2_neg_info <- smartd_rplc_neg_batch2_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

data1 <- 
  batch1_neg_info %>% 
  dplyr::select(., mz:rt)

data2 <- 
  batch2_neg_info %>% 
  dplyr::select(., mz:rt)

match_result <-
  sxtTools::sxtMTmatch(data1 = as.matrix(data1), 
                       data2 = as.matrix(data2), 
                       mz.tol = 25, rt.tol = 120, 
                       rt.error.type = "abs")

unique(match_result[,1])

###remove the 1:multiple
match_result <- 
  match_result %>% 
  as_tibble() %>% 
  group_by(Index1) %>% 
  filter(`rt error` == min(`rt error`)) %>% 
  ungroup()

name1 <- batch1_neg_info$name[match_result$Index1]
name2 <- batch2_neg_info$name[match_result$Index2]

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



rplc_pos_batch1_batch2_matched_info$batch1 <- 
  paste(rplc_pos_batch1_batch2_matched_info$batch1, "POS", sep = "_")

rplc_pos_batch1_batch2_matched_info$batch2 <- 
  paste(rplc_pos_batch1_batch2_matched_info$batch2, "POS", sep = "_")


rplc_neg_batch1_batch2_matched_info$batch1 <- 
  paste(rplc_neg_batch1_batch2_matched_info$batch1, "NEG", sep = "_")

rplc_neg_batch1_batch2_matched_info$batch2 <- 
  paste(rplc_neg_batch1_batch2_matched_info$batch2, "NEG", sep = "_")

write.csv(rplc_pos_batch1_batch2_matched_info, 
          'rplc_pos_batch1_batch2_matched_info.csv')

write.csv(rplc_neg_batch1_batch2_matched_info, 
          'rplc_neg_batch1_batch2_matched_info.csv')


#######get the identification table for batch 2
setwd("smartD_batch1_2/RPLC/POS_NEG/")
load("../POS/smartd_rplc_pos_batch2_5")
load("../NEG/smartd_rplc_neg_batch2_5")

ms1_data_rplc_batch2_pos <-
  smartd_rplc_pos_batch2_5@ms1.data[[1]]

ms1_data_rplc_batch2_neg <-
  smartd_rplc_neg_batch2_5@ms1.data[[1]]

ms1_data_rplc_batch2_pos$name <- 
  paste(ms1_data_rplc_batch2_pos$name, "POS", sep = "_")

ms1_data_rplc_batch2_neg$name <- 
  paste(ms1_data_rplc_batch2_neg$name, "NEG", sep = "_")


ms1_data_rplc_batch2 <- 
  dplyr::full_join(x = ms1_data_rplc_batch2_pos, 
                   y = ms1_data_rplc_batch2_neg, 
                   by = intersect(colnames(ms1_data_rplc_batch2_pos), colnames(ms1_data_rplc_batch2_neg)))


write.csv(ms1_data_rplc_batch2, "ms1_data_rplc_batch2.csv", 
          row.names = FALSE)


rplc_pos_batch1_batch2_matched_info
rplc_neg_batch1_batch2_matched_info

dim(rplc_pos_batch1_batch2_matched_info)
dim(rplc_neg_batch1_batch2_matched_info)

ms1_data_rplc_batch2 <- 
ms1_data_rplc_batch2 %>% 
  filter(name %in% c(rplc_pos_batch1_batch2_matched_info$batch2, rplc_neg_batch1_batch2_matched_info$batch2))

write.csv(ms1_data_rplc_batch2, "ms1_data_rplc_batch2.csv", 
          row.names = FALSE)


























