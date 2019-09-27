#-------------------------------------------------------------------------------
##batch 2 RPLC pos
#-------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
setwd("E:/project/smartD/smartD_batch2/RPLC/POS/data cleaning")

###construct sample info
sample_info <- readr::read_csv("sample.info.csv")
ms1_data <- readr::read_csv("Peak.table.csv")


colnames(sample_info) %>%
  head(.,10)

table(sample_info$participant_id)
colnames(ms1_data)
colnames(sample_info)[1:20]

sample_info <- data.frame("sample.name" = colnames(ms1_data),
                          "injection.order" = NA,
                          "class" = NA,
                          "batch" = NA,
                          "group" = NA,
                          stringsAsFactors = FALSE)

head(sample_info)
sample_info <- sample_info[-c(1:3),]
sample_info$injection.order <- 1:nrow(sample_info)
sample_info$class <- stringr::str_replace_all(string = sample_info$sample.name,
                                              pattern = "[0-9]", replacement = "")

sample_info$class[sample_info$class == "X" | sample_info$class == "P"] <- "Subject"
sample_info$class[sample_info$class == "QC."] <- "QC"
sample_info$class[grep("blk", sample_info$class)] <- "Blank"
sample_info$class[grep("QC\\.[A-C]{1}", sample_info$class)] <- "Blank"
sample_info$batch <- 1
sample_info$group <- sample_info$class

patient_info <- readr::read_csv("patient_info.csv")

patient_info <- 
  patient_info %>% 
  select(Sample_ID, `Visit GA`) %>% 
  mutate("Sample_ID2" = stringr::str_replace(Sample_ID, "SFU_B", "X")) %>% 
  select(-Sample_ID, GA = `Visit GA`)
 
# patient_info <- patient_info[,c("Sample_ID", "Visit GA")]

sample_info <-
  left_join(x = sample_info, y = patient_info, by = c("sample.name" = "Sample_ID2"))

sample_info %>% 
  filter(group == "Subject") %>% 
  select(sample.name, GA)

sample_info$GA[is.na(sample_info$GA)] <- 0

write.csv(sample_info, file = "sample_info.csv", row.names = FALSE)

###creat object class
setwd("E:/project/smartD/smartD_batch2/RPLC/POS/data cleaning")
library(metflow2)
(
  smartd_rplc_pos_batch2_1 <-
    creatMetflowObject(
      ms1.data = "Peak.table.csv",
      sample.information = "sample_info.csv",
      path = "."
    )
)

save(smartd_rplc_pos_batch2_1, file = "smartd_rplc_pos_batch2_1")


###remove peaks
(
  smartd_rplc_pos_batch2_2 <-
    filterPeak(
      object = smartd_rplc_pos_batch2_1,
      min.fraction.qc = 0.8,##QC at least more than 80%
      min.fraction = 0.5,
      min.subject.qc.ratio = 1,
      dl.qc.r2.cutoff = 0
    )
)


###remove some samples (outliers)
(
  smartd_rplc_pos_batch2_3 <-
    filterSample(object = smartd_rplc_pos_batch2_2,
                 min.fraction.peak = 0.5)
)

(
  plot <- 
    getMVplot4sample(object = smartd_rplc_pos_batch2_3)
)




###MV imputateion
(
  smartd_rplc_pos_batch2_4 <- 
    imputeMV(object = smartd_rplc_pos_batch2_3, method = 'knn')
)

save(smartd_rplc_pos_batch2_4, file = "smartd_rplc_pos_batch2_4")

# qc_data <- 
#   getData(object = smartd_rplc_pos_batch2_4, slot = "QC")
# 
# qc_data <- t(apply(qc_data, 1, function(x){
#   (x - mean(x))/sd(x)
# }))
# 
# qc_data <- tibble::as_tibble(qc_data)
# 
# qc_data <- 
#   gather(data = qc_data, colnames(qc_data), 
#          key = "QC", value = "Intensity")
# 
# plot <- 
#   ggplot(qc_data, aes(x = QC, y = Intensity)) +
#   geom_boxplot(colour = "#00468B7F", outlier.size = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2, colour = "red") +
#   labs(y = "Intensity (z-score)") +
#   coord_flip() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15))
# 
# plot
# export::graph2ppt(plot, "test", width = 4, heigh = 8)


###PCA analysis
# plot <- 
#   pcaAnalysis(object = smart.d.rplc.pos4, scale.method = "auto")
# (plot <- plot +
#     theme(legend.position = c(0, 1),legend.justification = c(0, 1),
#           legend.background = element_blank())
# )

# export::graph2ppt(plot, "test", width = 6, heigh = 4)


###RSD distribution
qc_rsd <- calRSD(object = smartd_rplc_pos_batch2_4, slot = "QC")

(plot1 <- 
    ggplot(data = qc_rsd, aes(index, rsd)) +
    geom_hline(yintercept = 30, linetype = 2, colour = "red") +
    geom_point() +
    scale_y_continuous(name = "Relative standard deviation (RSD, %)", limits = c(0,100),
                       breaks = c(0, 25, 30, 50, 75, 100), labels = c(c(0, 25, 30, 50, 75, 100)))+
    labs(x = "Peak index") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)



(plot2 <- 
    ggplot(data = qc_rsd, aes(x = rsd)) +
    geom_histogram(binwidth = 10, colour = "white", fill = "#0099B47F", origin = 0) +
    scale_x_continuous(name = "Relative standard deviation (RSD, %)") +
    labs(y = "Peak number") +
    # geom_vline(xintercept = 30, linetype = 2, colour = "red") +
    annotate(geom = "rect", xmin = 0, xmax = 30, ymin = 0, ymax = 380,
             alpha = 0.1, fill = "#ED00007F") +
    theme_bw()+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)


# g1 = ggplot2::ggplotGrob(plot1)
# (plot3 <- plot2 +
#     annotation_custom(grob = g1, xmin = 100, xmax = 500, ymin = 50, ymax = 350)
# )

# export::graph2ppt(plot3, "test", width = 6, heigh = 4)
####data normalization

(
  smartd_rplc_pos_batch2_5 <- 
    normalizeData(object = smartd_rplc_pos_batch2_4, 
                  method = "median")
)


save(smartd_rplc_pos_batch2_5, file = "smartd_rplc_pos_batch2_5")





#-------------------------------------------------------------------------------
##batch 2 RPLC neg
#-------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
setwd("E:/project/smartD/smartD_batch2/RPLC/NEG/data cleaning")

###construct sample info
ms1_data <- readr::read_csv("Peak.table.csv")

sample_info <- data.frame("sample.name" = colnames(ms1_data),
                          "injection.order" = NA,
                          "class" = NA,
                          "batch" = NA,
                          "group" = NA,
                          stringsAsFactors = FALSE)

head(sample_info)
sample_info <- sample_info[-c(1:3),]
sample_info$injection.order <- 1:nrow(sample_info)
sample_info$class <- stringr::str_replace_all(string = sample_info$sample.name,
                                              pattern = "[0-9]", replacement = "")

sample_info$class[sample_info$class == "X" | sample_info$class == "P"] <- "Subject"
sample_info$class[sample_info$class == "QC."] <- "QC"
sample_info$class[grep("blk", sample_info$class)] <- "Blank"
sample_info$class[grep("QC\\.[A-C]{1}", sample_info$class)] <- "Blank"
sample_info$batch <- 1
sample_info$group <- sample_info$class

patient_info <- readr::read_csv("patient_info.csv")

patient_info <- 
  patient_info %>% 
  select(Sample_ID, `Visit GA`) %>% 
  mutate("Sample_ID2" = stringr::str_replace(Sample_ID, "SFU_B", "X")) %>% 
  select(-Sample_ID, GA = `Visit GA`)

sample_info <-
  left_join(x = sample_info, y = patient_info, by = c("sample.name" = "Sample_ID2"))

sample_info %>% 
  filter(group == "Subject") %>% 
  select(sample.name, GA)

sample_info$GA[is.na(sample_info$GA)] <- 0

write.csv(sample_info, file = "sample_info.csv", row.names = FALSE)


###creat object class
setwd("E:/project/smartD/smartD_batch2/RPLC/NEG/data cleaning")
library(metflow2)

(
  smartd_rplc_neg_batch2_1 <-
    creatMetflowObject(
      ms1.data = "Peak.table.csv",
      sample.information = "sample_info.csv",
      path = "."
    )
)

save(smartd_rplc_neg_batch2_1, 
     file = "smartd_rplc_neg_batch2_1")


###remove peaks
(
  smartd_rplc_neg_batch2_2 <-
    filterPeak(
      object = smartd_rplc_neg_batch2_1,
      min.fraction.qc = 0.8,##QC at least more than 80%
      min.fraction = 0.5,
      min.subject.qc.ratio = 1,
      dl.qc.r2.cutoff = 0
    )
)

save(smartd_rplc_neg_batch2_2, 
     file = "smartd_rplc_neg_batch2_2")
###remove some samples (outliers)
(
  smartd_rplc_neg_batch2_3 <-
    filterSample(object = smartd_rplc_neg_batch2_2,
                 min.fraction.peak = 0.5)
)

(
  plot <- 
    getMVplot4sample(object = smartd_rplc_neg_batch2_3)
)

save(smartd_rplc_neg_batch2_3, 
     file = "smartd_rplc_neg_batch2_3")


###MV imputateion
(
  smartd_rplc_neg_batch2_4 <- 
    imputeMV(object = smartd_rplc_neg_batch2_3, method = 'knn')
)

save(smartd_rplc_neg_batch2_4, file = "smartd_rplc_neg_batch2_4")

# qc_data <- 
#   getData(object = smartd_rplc_neg_batch2_4, slot = "QC")
# 
# qc_data <- t(apply(qc_data, 1, function(x){
#   (x - mean(x))/sd(x)
# }))
# 
# qc_data <- tibble::as_tibble(qc_data)
# 
# qc_data <- 
#   gather(data = qc_data, colnames(qc_data), 
#          key = "QC", value = "Intensity")
# 
# plot <- 
#   ggplot(qc_data, aes(x = QC, y = Intensity)) +
#   geom_boxplot(colour = "#00468B7F", outlier.size = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2, colour = "red") +
#   labs(y = "Intensity (z-score)") +
#   coord_flip() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15))
# 
# plot
# export::graph2ppt(plot, "test", width = 4, heigh = 8)


###PCA analysis
# plot <- 
#   pcaAnalysis(object = smart.d.rplc.neg4, scale.method = "auto")
# (plot <- plot +
#     theme(legend.negition = c(0, 1),legend.justification = c(0, 1),
#           legend.background = element_blank())
# )

# export::graph2ppt(plot, "test", width = 6, heigh = 4)


###RSD distribution
qc_rsd <- 
  calRSD(object = smartd_rplc_neg_batch2_4, slot = "QC")

(plot1 <- 
    ggplot(data = qc_rsd, aes(index, rsd)) +
    geom_hline(yintercept = 30, linetype = 2, colour = "red") +
    geom_point() +
    scale_y_continuous(name = "Relative standard deviation (RSD, %)", limits = c(0,100),
                       breaks = c(0, 25, 30, 50, 75, 100), labels = c(c(0, 25, 30, 50, 75, 100)))+
    labs(x = "Peak index") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)



(plot2 <- 
    ggplot(data = qc_rsd, aes(x = rsd)) +
    geom_histogram(binwidth = 10, colour = "white", fill = "#0099B47F", origin = 0) +
    scale_x_continuous(name = "Relative standard deviation (RSD, %)") +
    labs(y = "Peak number") +
    # geom_vline(xintercept = 30, linetype = 2, colour = "red") +
    annotate(geom = "rect", xmin = 0, xmax = 30, ymin = 0, ymax = 380,
             alpha = 0.1, fill = "#ED00007F") +
    theme_bw()+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)


# g1 = ggplot2::ggplotGrob(plot1)
# (plot3 <- plot2 +
#     annotation_custom(grob = g1, xmin = 100, xmax = 500, ymin = 50, ymax = 350)
# )

# export::graph2ppt(plot3, "test", width = 6, heigh = 4)
####data normalization

(
  smartd_rplc_neg_batch2_5 <- 
    normalizeData(object = smartd_rplc_neg_batch2_4, 
                  method = "median")
)


save(smartd_rplc_neg_batch2_5, file = "smartd_rplc_neg_batch2_5")


smartd_rplc_neg_batch2_5
sample.info <- smartd_rplc_neg_batch2_5@sample.info



#-------------------------------------------------------------------------------
##batch 2 HILIC pos
#-------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
setwd("E:/project/smartD/smartD_batch2/HILIC/POS/data cleaning")
###construct sample info
ms1_data <- readr::read_csv("Peak.table.csv")

sample_info <- data.frame("sample.name" = colnames(ms1_data),
                          "injection.order" = NA,
                          "class" = NA,
                          "batch" = NA,
                          "group" = NA,
                          stringsAsFactors = FALSE)

head(sample_info)
sample_info <- sample_info[-c(1:3),]
sample_info$injection.order <- 1:nrow(sample_info)
sample_info$class <- stringr::str_replace_all(string = sample_info$sample.name,
                                              pattern = "[0-9]", replacement = "")


sample_info$class[sample_info$class == "SFU_B" | sample_info$class == "SFU_B_"] <- "Subject"
sample_info$class[sample_info$class == "QC_"] <- "QC"
sample_info$class[sample_info$class == "Batch_QC"] <- "Batch1_QC"

sample_info$class[grep("Blank", sample_info$class)] <- "Blank"
sample_info$class[grep("QC_DL", sample_info$class)] <- "QC_DL"
sample_info$batch <- 1
sample_info$group <- sample_info$class

patient_info <- readr::read_csv("patient_info.csv")
patient_info2 <- readr::read_csv("SmartD_all346urine.csv")

patient_info <- 
  patient_info %>% 
  select(Sample_ID, GA = `Visit GA`) 


patient_info2 <- 
  patient_info2 %>% 
  select(sample_id_HILIC, GA = `Visit.GA`) 


temp <- 
  left_join(x = patient_info, 
            y = patient_info2, 
            by = c("Sample_ID" = "sample_id_HILIC"))


sample_info <-
  left_join(x = sample_info, y = patient_info2, by = c("sample.name" = "sample_id_HILIC"))

sample_info %>% 
  filter(group == "Subject") %>% 
  select(sample.name, GA)

sample_info$GA[is.na(sample_info$GA)] <- 0

write.csv(sample_info, file = "sample_info.csv", row.names = FALSE)


###creat object class
setwd("E:/project/smartD/smartD_batch2/HILIC/POS/data cleaning")
library(metflow2)
(
  smartd_hilic_pos_batch2_1 <-
    creatMetflowObject(
      ms1.data = "Peak.table.csv",
      sample.information = "sample_info.csv",
      path = "."
    )
)

save(smartd_hilic_pos_batch2_1, file = "smartd_hilic_pos_batch2_1")


###remove peaks
(
  smartd_hilic_pos_batch2_2 <-
    filterPeak(
      object = smartd_hilic_pos_batch2_1,
      min.fraction.qc = 0.8,##QC at least more than 80%
      min.fraction = 0.5,
      min.subject.qc.ratio = 1,
      dl.qc.r2.cutoff = 0
    )
)

save(smartd_hilic_pos_batch2_2, file = "smartd_hilic_pos_batch2_2")

###remove some samples (outliers)
(
  smartd_hilic_pos_batch2_3 <-
    filterSample(object = smartd_hilic_pos_batch2_2,
                 min.fraction.peak = 0.5)
)

(
  plot <- 
    getMVplot4sample(object = smartd_hilic_pos_batch2_3)
)


save(smartd_hilic_pos_batch2_3, file = "smartd_hilic_pos_batch2_3")

###MV imputateion
(
  smartd_hilic_pos_batch2_4 <- 
    imputeMV(object = smartd_hilic_pos_batch2_3, method = 'knn')
)

save(smartd_hilic_pos_batch2_4, file = "smartd_hilic_pos_batch2_4")

# qc_data <- 
#   getData(object = smartd_hilic_pos_batch2_4, slot = "QC")
# 
# qc_data <- t(apply(qc_data, 1, function(x){
#   (x - mean(x))/sd(x)
# }))
# 
# qc_data <- tibble::as_tibble(qc_data)
# 
# qc_data <- 
#   gather(data = qc_data, colnames(qc_data), 
#          key = "QC", value = "Intensity")
# 
# plot <- 
#   ggplot(qc_data, aes(x = QC, y = Intensity)) +
#   geom_boxplot(colour = "#00468B7F", outlier.size = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2, colour = "red") +
#   labs(y = "Intensity (z-score)") +
#   coord_flip() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15))
# 
# plot
# export::graph2ppt(plot, "test", width = 4, heigh = 8)


###PCA analysis
# plot <- 
#   pcaAnalysis(object = smart.d.hilic.pos4, scale.method = "auto")
# (plot <- plot +
#     theme(legend.position = c(0, 1),legend.justification = c(0, 1),
#           legend.background = element_blank())
# )

# export::graph2ppt(plot, "test", width = 6, heigh = 4)


###RSD distribution
qc_rsd <- calRSD(object = smartd_hilic_pos_batch2_4, slot = "QC")

(plot1 <- 
    ggplot(data = qc_rsd, aes(index, rsd)) +
    geom_hline(yintercept = 30, linetype = 2, colour = "red") +
    geom_point() +
    scale_y_continuous(name = "Relative standard deviation (RSD, %)", limits = c(0,100),
                       breaks = c(0, 25, 30, 50, 75, 100), labels = c(c(0, 25, 30, 50, 75, 100)))+
    labs(x = "Peak index") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)



(plot2 <- 
    ggplot(data = qc_rsd, aes(x = rsd)) +
    geom_histogram(binwidth = 10, colour = "white", fill = "#0099B47F", origin = 0) +
    scale_x_continuous(name = "Relative standard deviation (RSD, %)") +
    labs(y = "Peak number") +
    # geom_vline(xintercept = 30, linetype = 2, colour = "red") +
    annotate(geom = "rect", xmin = 0, xmax = 30, ymin = 0, ymax = 380,
             alpha = 0.1, fill = "#ED00007F") +
    theme_bw()+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)


# g1 = ggplot2::ggplotGrob(plot1)
# (plot3 <- plot2 +
#     annotation_custom(grob = g1, xmin = 100, xmax = 500, ymin = 50, ymax = 350)
# )

# export::graph2ppt(plot3, "test", width = 6, heigh = 4)
####data normalization

(
  smartd_hilic_pos_batch2_5 <- 
    normalizeData(object = smartd_hilic_pos_batch2_4, 
                  method = "median")
)


save(smartd_hilic_pos_batch2_5, file = "smartd_hilic_pos_batch2_5")


qc_rsd <- calRSD(object = smartd_hilic_pos_batch2_5, slot = "QC")

(plot1 <- 
    ggplot(data = qc_rsd, aes(index, rsd)) +
    geom_hline(yintercept = 30, linetype = 2, colour = "red") +
    geom_point() +
    scale_y_continuous(name = "Relative standard deviation (RSD, %)", limits = c(0,100),
                       breaks = c(0, 25, 30, 50, 75, 100), labels = c(c(0, 25, 30, 50, 75, 100)))+
    labs(x = "Peak index") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)


#-------------------------------------------------------------------------------
##batch 2 HILIC neg
#-------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
setwd("E:/project/smartD/smartD_batch2/HILIC/NEG/data cleaning")

###construct sample info
ms1_data <- readr::read_csv("Peak.table.csv")

sample_info <- data.frame("sample.name" = colnames(ms1_data),
                          "injection.order" = NA,
                          "class" = NA,
                          "batch" = NA,
                          "group" = NA,
                          stringsAsFactors = FALSE)

head(sample_info)
sample_info <- sample_info[-c(1:3),]
sample_info$injection.order <- 1:nrow(sample_info)
sample_info$class <- stringr::str_replace_all(string = sample_info$sample.name,
                                              pattern = "[0-9]", replacement = "")


sample_info$class[sample_info$class == "SFU_B" | sample_info$class == "SFU_B_"] <- "Subject"
sample_info$class[sample_info$class == "QC_"] <- "QC"
sample_info$class[sample_info$class == "Batch_QC"] <- "Batch1_QC"

sample_info$class[grep("Blank", sample_info$class)] <- "Blank"
sample_info$class[grep("QC_DL", sample_info$class)] <- "QC_DL"
sample_info$batch <- 1
sample_info$group <- sample_info$class

patient_info <- readr::read_csv("patient_info.csv")

patient_info <- 
  patient_info %>% 
  select(Sample_ID, GA = `Visit GA`) 

sample_info <-
  left_join(x = sample_info, y = patient_info, by = c("sample.name" = "Sample_ID"))

sample_info %>% 
  filter(group == "Subject") %>% 
  select(sample.name, GA)

sample_info$GA[is.na(sample_info$GA)] <- 0

write.csv(sample_info, file = "sample_info.csv", row.names = FALSE)


###creat object class
setwd("E:/project/smartD/smartD_batch2/HILIC/NEG/data cleaning")
library(metflow2)

(
  smartd_hilic_neg_batch2_1 <-
    creatMetflowObject(
      ms1.data = "Peak.table.csv",
      sample.information = "sample_info.csv",
      path = "."
    )
)

save(smartd_hilic_neg_batch2_1, 
     file = "smartd_hilic_neg_batch2_1")


###remove peaks
(
  smartd_hilic_neg_batch2_2 <-
    filterPeak(
      object = smartd_hilic_neg_batch2_1,
      min.fraction.qc = 0.8,##QC at least more than 80%
      min.fraction = 0.5,
      min.subject.qc.ratio = 1,
      dl.qc.r2.cutoff = 0
    )
)

save(smartd_hilic_neg_batch2_2, 
     file = "smartd_hilic_neg_batch2_2")
###remove some samples (outliers)
(
  smartd_hilic_neg_batch2_3 <-
    filterSample(object = smartd_hilic_neg_batch2_2,
                 min.fraction.peak = 0.5)
)

(
  plot <- 
    getMVplot4sample(object = smartd_hilic_neg_batch2_3)
)

save(smartd_hilic_neg_batch2_3, 
     file = "smartd_hilic_neg_batch2_3")


###MV imputateion
(
  smartd_hilic_neg_batch2_4 <- 
    imputeMV(object = smartd_hilic_neg_batch2_3, method = 'knn')
)

save(smartd_hilic_neg_batch2_4, file = "smartd_hilic_neg_batch2_4")

# qc_data <- 
#   getData(object = smartd_hilic_neg_batch2_4, slot = "QC")
# 
# qc_data <- t(apply(qc_data, 1, function(x){
#   (x - mean(x))/sd(x)
# }))
# 
# qc_data <- tibble::as_tibble(qc_data)
# 
# qc_data <- 
#   gather(data = qc_data, colnames(qc_data), 
#          key = "QC", value = "Intensity")
# 
# plot <- 
#   ggplot(qc_data, aes(x = QC, y = Intensity)) +
#   geom_boxplot(colour = "#00468B7F", outlier.size = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2, colour = "red") +
#   labs(y = "Intensity (z-score)") +
#   coord_flip() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15))
# 
# plot
# export::graph2ppt(plot, "test", width = 4, heigh = 8)


###PCA analysis
# plot <- 
#   pcaAnalysis(object = smart.d.hilic.neg4, scale.method = "auto")
# (plot <- plot +
#     theme(legend.negition = c(0, 1),legend.justification = c(0, 1),
#           legend.background = element_blank())
# )

# export::graph2ppt(plot, "test", width = 6, heigh = 4)


###RSD distribution
qc_rsd <- 
  calRSD(object = smartd_hilic_neg_batch2_4, slot = "QC")

(plot1 <- 
    ggplot(data = qc_rsd, aes(index, rsd)) +
    geom_hline(yintercept = 30, linetype = 2, colour = "red") +
    geom_point() +
    scale_y_continuous(name = "Relative standard deviation (RSD, %)", limits = c(0,100),
                       breaks = c(0, 25, 30, 50, 75, 100), labels = c(c(0, 25, 30, 50, 75, 100)))+
    labs(x = "Peak index") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)



(plot2 <- 
    ggplot(data = qc_rsd, aes(x = rsd)) +
    geom_histogram(binwidth = 10, colour = "white", fill = "#0099B47F", origin = 0) +
    scale_x_continuous(name = "Relative standard deviation (RSD, %)") +
    labs(y = "Peak number") +
    # geom_vline(xintercept = 30, linetype = 2, colour = "red") +
    annotate(geom = "rect", xmin = 0, xmax = 30, ymin = 0, ymax = 380,
             alpha = 0.1, fill = "#ED00007F") +
    theme_bw()+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)


# g1 = ggplot2::ggplotGrob(plot1)
# (plot3 <- plot2 +
#     annotation_custom(grob = g1, xmin = 100, xmax = 500, ymin = 50, ymax = 350)
# )

# export::graph2ppt(plot3, "test", width = 6, heigh = 4)
####data normalization

(
  smartd_hilic_neg_batch2_5 <- 
    normalizeData(object = smartd_hilic_neg_batch2_4, 
                  method = "median")
)


save(smartd_hilic_neg_batch2_5, file = "smartd_hilic_neg_batch2_5")


smartd_hilic_neg_batch2_5
sample.info <- smartd_hilic_neg_batch2_5@sample.info

