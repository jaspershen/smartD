###This script is for data cleaning for batch 1 and batch 2. Batch 1 and batch 2 data are processed using XCMS
##together
setwd("smartD_batch1_2/RPLC_xcms/POS/data_cleaning/")

###construct sample information
# sample_info1 <-
#   readr::read_csv("sample_info1.csv")
# 
# sample_info2 <-
#   readr::read_csv("sample_info2.csv")
# 
# 
# sample_info1$sample.name
# sample_info2$sample.name
# 
# peak_table <- 
#   readr::read_csv("Peak.table.csv")
# 
# 
# intersect(c(sample_info1$sample.name, sample_info2$sample.name),
#           colnames(peak_table))
# 
# 
# setdiff(c(sample_info1$sample.name, sample_info2$sample.name),
#         colnames(peak_table))
# 
# setdiff(colnames(peak_table),
#         c(sample_info1$sample.name, 
#           sample_info2$sample.name)
# )
# 
# sample_info2$batch <- 2
# 
# ###remove QC_DL samples
# 
# sample_info1 <-
#   sample_info1 %>% 
#   filter(., !stringr::str_detect(sample.name, "DL"))
# 
# sample_info2 <-
#   sample_info2 %>% 
#   filter(., !stringr::str_detect(sample.name, "DL"))
# 
# sample_info <- 
#   rbind(sample_info1, sample_info2)
# 
# 
# peak_table <- 
#   peak_table %>% 
#   select(., -contains("DL"))
# 
# setdiff(sample_info$sample.name,
#         colnames(peak_table))
# 
# setdiff(colnames(peak_table),
#         sample_info$sample.name
# )
# 
# ##remove "QC2.4" and "QC2.5" from peak_table
# peak_table <- 
#   peak_table %>% 
#   select(-c(QC2.4, QC2.5))
# 
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# write.csv(sample_info, "sample_info.csv", row.names = FALSE)


##RPLC positive
###creat object class
setwd("smartD_batch1_2/RPLC_xcms/POS/data_cleaning/")
library(metflow2)

(
  smartd_rplc_pos_1 <-
    creatMetflowObject(
      ms1.data = "peak_table.csv",
      sample.information = "sample_info.csv",
      path = "."
    )
)

save(smartd_rplc_pos_1, file = "smartd_rplc_pos_1")


###remove peaks
(
  smartd_rplc_pos_2 <-
    filterPeak(
      object = smartd_rplc_pos_1,
      min.fraction.qc = 0.8,##QC at least more than 80%
      min.fraction = 0.5,
      min.subject.qc.ratio = 1,
      dl.qc.r2.cutoff = 0
    )
)


###remove some samples (outliers)
(
  smartd_rplc_pos_3 <-
    filterSample(object = smartd_rplc_pos_2,
                 min.fraction.peak = 0.5)
)

(
  plot <- 
    getMVplot4sample(object = smartd_rplc_pos_3)
)




###MV imputateion
(
  smartd_rplc_pos_4 <- 
    imputeMV(object = smartd_rplc_pos_3, method = 'knn')
)

save(smartd_rplc_pos_4, file = "smartd_rplc_pos_4")

# qc_data <- 
#   getData(object = smartd_rplc_pos_batch1_4, slot = "QC")
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
qc_rsd <- calRSD(object = smartd_rplc_pos_4, slot = "QC")

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
  smartd_rplc_pos_5 <- 
    normalizeData(object = smartd_rplc_pos_4, 
                  method = "mean")
)


###batch intergration

(
  smartd_rplc_pos_5 <- 
    metflow2:::integrateData(
      object = smartd_rplc_pos_5, 
      method = "subject.mean")
)

qc_rsd <- calRSD(object = smartd_rplc_pos_5, 
                 slot = "QC")

plot(qc_rsd$rsd)
save(smartd_rplc_pos_5, file = "smartd_rplc_pos_5")









##RPLC negative
###creat object class
# setwd("smartD_batch1_2/RPLC_xcms/NEG/data_cleaning/")
# sample_info <- 
#   readr::read_csv("sample_info.csv")
# 
# peak_table <- 
#   readr::read_csv("Peak.table.csv")
# 
# setdiff(colnames(peak_table),
#         sample_info$sample.name)
# 
# setdiff(sample_info$sample.name,
#         colnames(peak_table)
# )
# 
# 
# ###remove the samples which are not in sample_info or 
# name1 <-
#   setdiff(colnames(peak_table),
#           sample_info$sample.name) 
# 
# name1 <-
#   name1[!name1 %in% c("name", "mz", "rt")]
# 
# peak_table <- 
#   peak_table %>% 
#   select(-name1)
# 
# 
# name2 <-
#   setdiff(sample_info$sample.name,
#           colnames(peak_table)
#   )
#   
# sample_info <- 
#   sample_info %>% 
#   filter(!sample.name%in%name2)
# 
# 
# write.csv(sample_info, "sample_info.csv", row.names = FALSE)  
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)


library(metflow2)

(
  smartd_rplc_neg_1 <-
    creatMetflowObject(
      ms1.data = "peak_table.csv",
      sample.information = "sample_info.csv",
      path = "."
    )
)

save(smartd_rplc_neg_1, file = "smartd_rplc_neg_1")


###remove peaks
(
  smartd_rplc_neg_2 <-
    filterPeak(
      object = smartd_rplc_neg_1,
      min.fraction.qc = 0.8,##QC at least more than 80%
      min.fraction = 0.5,
      min.subject.qc.ratio = 1,
      dl.qc.r2.cutoff = 0
    )
)


###remove some samples (outliers)
(
  smartd_rplc_neg_3 <-
    filterSample(object = smartd_rplc_neg_2,
                 min.fraction.peak = 0.5)
)

(
  plot <- 
    getMVplot4sample(object = smartd_rplc_neg_3)
)




###MV imputateion
(
  smartd_rplc_neg_4 <- 
    imputeMV(object = smartd_rplc_neg_3, method = 'knn')
)

save(smartd_rplc_neg_4, file = "smartd_rplc_neg_4")

# qc_data <- 
#   getData(object = smartd_rplc_neg_batch1_4, slot = "QC")
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
qc_rsd <- calRSD(object = smartd_rplc_neg_4, slot = "QC")

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
  smartd_rplc_neg_5 <- 
    normalizeData(object = smartd_rplc_neg_4, 
                  method = "mean")
)


###batch intergration

(
  smartd_rplc_neg_5 <- 
    metflow2:::integrateData(
      object = smartd_rplc_neg_5, 
      method = "subject.mean")
)

qc_rsd <- calRSD(object = smartd_rplc_neg_5, 
                 slot = "QC")
plot(qc_rsd$rsd)

save(smartd_rplc_neg_5, file = "smartd_rplc_neg_5")





