#-------------------------------------------------------------------------------
###RPLC pos


setwd("E:/project/smartD/RPLC/POS/data cleaning")
# # sample_info <- readr::read_csv("sample.info.csv")
# # ms1_data <- readr::read_csv("Peak.table.csv")
# # library(magrittr)
# # colnames(sample.info)
# # colnames(sample_info) %>%
# #   head(.,10)
# #
# # table(sample_info$participant_id)
# # colnames(ms1_data)
# # colnames(sample_info)[1:20]
# 
# sample_info <- data.frame("sample.name" = colnames(ms1_data),
#                           "injection.order" = NA,
#                           "class" = NA,
#                           "batch" = NA,
#                           "group" = NA,
#                           stringsAsFactors = FALSE)
# 
# head(sample_info)
# sample_info <- sample_info[-c(1:3),]
# sample_info$injection.order <- 1:nrow(sample_info)
# sample_info$class <- stringr::str_replace_all(string = sample_info$sample.name,
#                                               pattern = "[0-9]", replacement = "")
# 
# sample_info$class[sample_info$class == "X" | sample_info$class == "P"] <- "Subject"
# sample_info$class[sample_info$class == "QC."] <- "QC"
# sample_info$class[grep("blk", sample_info$class)] <- "Blank"
# sample_info$class[grep("QC\\.[A-C]{1}", sample_info$class)] <- "Blank"
# sample_info$batch <- 1
# sample_info$group <- sample_info$class
# 
# write.csv(sample_info, file = "sample_info.csv", row.names = FALSE)

###creat object class
setwd("E:/project/smartD/RPLC/POS/data cleaning")
library(metflow2)
(smart.d.rplc.pos <- creatMetflowObject(ms1.data = "Peak.table.csv",
                                        sample.information = "sample_info.csv",
                                        path = ".")
)

qc_data <- getData(object = smart.d.rplc.pos, slot = "QC")
na.fraction <- apply(qc_data, 1, function(x){
  sum(is.na(x) * 100/ncol(qc_data))
})

class <- rep("NO", length(na.fraction))
class[which(na.fraction > 20)] <- "YES"
na.fraction <- data.frame(index = 1:length(na.fraction),
                          na.fraction,
                          class,
                          stringsAsFactors = FALSE)
require(ggplot2)
plot <- ggplot(data = na.fraction) +
  geom_point(aes(x = index, y = na.fraction, colour = class), size = 2) +
  labs(x = "Peak index", y = "Missing value ratio (%)") +
  scale_colour_manual(values = c("#00A087FF", "#E64B35FF")) +
  guides(colour = FALSE) +
  geom_hline(yintercept = 20, color = "red", linetype = 2) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))

plot
# export::graph2ppt(plot, "test")


# y1 <- as.numeric(y1)/max(as.numeric(y1))
# y2 <- as.numeric(y2)/max(as.numeric(y2))
# y3 <- as.numeric(y3)/max(as.numeric(y3))
# y4 <- as.numeric(y4)/max(as.numeric(y4))
# 
# y <- as.numeric(c(y1, y2, y3, y4))
# class <- rep(c("Peak1", "Peak2", "Peak3", "Peak4"),4)
# class <- sort(class)
# x <- rep(c(1,2,4,8), 4)
# 
# temp <- data.frame(x, y, class, stringsAsFactors = FALSE)
# 
# 
# plot <- ggplot(temp, aes(x, y)) +
#   geom_point(colour = "black", size = 2) +
#   geom_smooth(method = "lm", se = FALSE, colour = "#E64B35FF", linetype = 2) +
#   labs(x = 'Times of dilution', y = 'Relative intensity') +
#   facet_wrap(~class, dir = "v") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15))
# 
# export::graph2ppt(plot, "test")
# 
# 
# lm(formula = y~x, data = dplyr::filter(.data = temp, class == "Peak4")) %>%
#   summary(.) %>%
#   `[[`(., "r.squared")


###remove peaks
smart.d.rplc.pos2 <- filterPeak(object = smart.d.rplc.pos,
                                min.fraction.qc = 0.8,
                                min.fraction = 0.5,
                                min.subject.qc.ratio = 1,
                                dl.qc.r2.cutoff = 0.7)

smart.d.rplc.pos2

smart.d.rplc.pos3 <- filterSample(object = smart.d.rplc.pos2, 
                                  min.fraction.peak = 0.5)

(plot <- getMVplot4sample(object = smart.d.rplc.pos3))

# export::graph2ppt(plot, "text", width = 8, height = 6)



###MV imputateion
smart.d.rplc.pos4 <- imputeMV(object = smart.d.rplc.pos3, method = 'knn')
smart.d.rplc.pos4

qc_data <- getData(object = smart.d.rplc.pos4, slot = "QC")

qc_data <- t(apply(qc_data, 1, function(x){
  (x - mean(x))/sd(x)
}))

qc_data <- tibble::as_tibble(qc_data)

qc_data <- gather(data = qc_data, colnames(qc_data), key = "QC", value = "Intensity")

plot <- ggplot(qc_data, aes(x = QC, y = Intensity)) +
  geom_boxplot(colour = "#00468B7F", outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  labs(y = "Intensity (z-score)") +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))

plot
# export::graph2ppt(plot, "test", width = 4, heigh = 8)


###PCA analysis
plot <- pcaAnalysis(object = smart.d.rplc.pos4, scale.method = "auto")
(plot <- plot +
    theme(legend.position = c(0, 1),legend.justification = c(0, 1),
          legend.background = element_blank())
)



# export::graph2ppt(plot, "test", width = 6, heigh = 4)


###RSD distribution
qc_rsd <- calRSD(object = smart.d.rplc.pos4, slot = "QC")

(plot1 <- ggplot(data = qc_rsd, aes(index, rsd)) +
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



(plot2 <- ggplot(data = qc_rsd, aes(x = rsd)) +
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


g1 = ggplot2::ggplotGrob(plot1)
(plot3 <- plot2 +
    annotation_custom(grob = g1, xmin = 100, xmax = 500, ymin = 50, ymax = 350)
)

# export::graph2ppt(plot3, "test", width = 6, heigh = 4)
####data normalization
(smart.d.rplc.pos5 <- normalizeData(object = smart.d.rplc.pos4, method = "median"))


qc_data <- getData(object = smart.d.rplc.pos5, slot = "QC")

qc_data <- t(apply(qc_data, 1, function(x){
  (x - mean(x))/sd(x)
}))

qc_data <- tibble::as_tibble(qc_data)

qc_data <- tidyr::gather(data = qc_data, colnames(qc_data), key = "QC", value = "Intensity")
require(ggplot2)
(plot <- ggplot(qc_data, aes(x = QC, y = Intensity)) +
    geom_boxplot(colour = "#00468B7F", outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = 2, colour = "red") +
    labs(y = "Intensity (z-score)") +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)

# export::graph2ppt(plot, "test", width = 4, heigh = 8)



###PCA analysis
plot <- pcaAnalysis(object = smart.d.rplc.pos5, scale.method = "auto")
(plot <- plot +
    theme(legend.position = c(0, 1),legend.justification = c(0, 1),
          legend.background = element_blank())
)


# export::graph2ppt(plot, "test", width = 6, heigh = 4)


###RSD distribution
qc_rsd <- calRSD(object = smart.d.rplc.pos5, slot = "QC")
qc_rsd <- dplyr::mutate(qc_rsd, index = 1:nrow(qc_rsd))

(plot1 <- ggplot(data = qc_rsd, aes(index, rsd)) +
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



(plot2 <- ggplot(data = qc_rsd, aes(x = rsd)) +
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


g1 = ggplotGrob(plot1)
(plot3 <- plot2 +
    annotation_custom(grob = g1, xmin = 100, xmax = 500, ymin = 50, ymax = 350)
)


# export::graph2ppt(plot3, "test", width = 6, heigh = 4)



qc_rsd1 <- calRSD(object = smart.d.rplc.pos4, slot = 'QC')
qc_rsd2 <- calRSD(object = smart.d.rplc.pos5, slot = 'QC')
colnames(qc_rsd1)[3] <- "before"
colnames(qc_rsd2)[3] <- "after"

qc_rsd <- dplyr::left_join(x = qc_rsd1, qc_rsd2)

qc_rsd <-
  qc_rsd %>%
  mutate(., class = ifelse(after - before >= 0, "Increase", "Decrease"))


(plot <-
    ggplot(qc_rsd, aes(after, before)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "black") +
    geom_point(aes(colour = class), alpha = 0.5, size = 2) +
    scale_color_manual(values = c("#42B5404C", "#ED00004C")) +
    theme_bw()+
    scale_x_continuous(name = "After data normalization (RSD, %)", limits = c(0,300)) +
    scale_y_continuous(name = "Before data normalization (RSD, %)", limits = c(0,300)) +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15),
          legend.position = c(0, 1), legend.justification = c(0,1),
          legend.background = element_blank())
)



# export::graph2ppt(plot, "test", width = 10, heigh = 4)




qc_data1 <- getData(object = smart.d.rplc.pos4, slot = "QC")

qc_data1 <- t(apply(qc_data1, 1, function(x){
  (x - mean(x))/sd(x)
}))

qc_data1 <- tibble::as_tibble(qc_data1)

qc_data1 <- gather(data = qc_data1, colnames(qc_data1), key = "QC", value = "Intensity")

qc_data1 <-
  qc_data1 %>%
  mutate(., Class = "Before")


qc_data2 <- getData(object = smart.d.rplc.pos5, slot = "QC")

qc_data2 <- t(apply(qc_data2, 1, function(x){
  (x - mean(x))/sd(x)
}))

qc_data2 <- tibble::as_tibble(qc_data2)

qc_data2 <- gather(data = qc_data2, colnames(qc_data2), key = "QC", value = "Intensity")

qc_data2 <-
  qc_data2 %>%
  mutate(., Class = "After")

qc_data <- rbind(qc_data1, qc_data2)
qc_data$Class <- factor(qc_data$Class, levels = c("Before", "After"))

(plot <- ggplot(qc_data, aes(x = QC, y = Intensity)) +
    geom_boxplot(mapping = aes(, fill = Class),colour = "#00468B7F", outlier.size = 0.5) +
    scale_fill_manual(values = c("#ED0000B2", "#42B540B2")) +
    geom_hline(yintercept = 0, linetype = 2, colour = "black") +
    labs(y = "Intensity (z-score)") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text.x = element_text(size = 10, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_blank(),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)


export::graph2ppt(plot, "test", width = 10, heigh = 4)




####metabolite identificatio
library(metIdentify)

#####RPLC positive 25 V
setwd("E:/project/smartD/RPLC/POS/metabolite identification/NCE25")
parameter1 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "msDatabase_rplc0.0.1",
                               threads = 5)

parameter2 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "hmdbDatabase0.0.1",
                               threads = 5)

parameter3 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "massbankDatabase0.0.1",
                               threads = 5)

parameter4 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "metlinDatabase0.0.1",
                               threads = 5)

parameter5 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "monaDatabase0.0.1",
                               threads = 5)

parameter6 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "nistDatabase0.0.1",
                               threads = 5)

parameter7 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "orbitrapDatabase0.0.1",
                               threads = 5)

parameter8 <- mzIdentifyParam(ms1.match.ppm = 25,
                              polarity = "positive",
                              column = "rp",
                              candidate.num = 3,
                              database = "HMDB.metabolite.data",
                              threads = 5)
ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mzXML", dir(), value = TRUE)
result.pRPLC.nce25 <- metIdentify4all(ms1.data = ms1.data,
                                      ms2.data = ms2.data,
                                      parameter.list = c(parameter1,
                                                         parameter2,
                                                         parameter3,
                                                         parameter4,
                                                         parameter5,
                                                         parameter6,
                                                         parameter7,
                                                         parameter8),
                                      path = ".")

save(result.pRPLC.nce25, file = "result.pRPLC.nce25")
print(1)


#####RPLC positive 50 V
setwd("E:/project/smartD/RPLC/POS/metabolite identification/NCE50")

ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mzXML", dir(), value = TRUE)
result.pRPLC.nce50 <- metIdentify4all(ms1.data = ms1.data,
                                      ms2.data = ms2.data,
                                      parameter.list = c(parameter1,
                                                         parameter2,
                                                         parameter3,
                                                         parameter4,
                                                         parameter5,
                                                         parameter6,
                                                         parameter7,
                                                         parameter8),
                                      path = ".")

save(result.pRPLC.nce50, file = "result.pRPLC.nce50")
print(1)





identification.table1 <- getIdentificationTable(result.pRPLC.nce25[[1]],
                                                result.pRPLC.nce50[[1]],
                                                candidate.num = 3,
                                                type = "old")

identification.table1 <-
  identification.table1[which(!is.na(identification.table1$Identification)),]
identification.table1 <- data.frame(identification.table1, Level = 1, stringsAsFactors = FALSE)
dim(identification.table1)

identification.table2 <- getIdentificationTable(result.pRPLC.nce25[[2]],
                                                result.pRPLC.nce25[[3]],
                                                result.pRPLC.nce25[[4]],
                                                result.pRPLC.nce25[[5]],
                                                result.pRPLC.nce25[[6]],
                                                result.pRPLC.nce25[[7]],
                                                result.pRPLC.nce50[[2]],
                                                result.pRPLC.nce50[[3]],
                                                result.pRPLC.nce50[[4]],
                                                result.pRPLC.nce50[[5]],
                                                result.pRPLC.nce50[[6]],
                                                result.pRPLC.nce50[[7]],
                                                candidate.num = 3,
                                                type = "old")

identification.table2 <- dplyr::filter(identification.table2, !is.na(Identification))
require(tidyverse)
identification.table2 <- identification.table2 %>%
  dplyr::filter(., !(name %in% identification.table1$name))

identification.table2 <- data.frame(identification.table2,
                                    Level = 2, stringsAsFactors = FALSE)
dim(identification.table2)


identification.table3 <- getIdentificationTable2(result.pRPLC.nce25[[8]],
                                                 candidate.num = 3,
                                                 type = "old")

identification.table3 <- dplyr::filter(identification.table3, !is.na(Identification))

identification.table3 <- identification.table3 %>%
  dplyr::filter(., !(name %in% c(identification.table1$name, identification.table2$name)))


identification.table3 <- data.frame(identification.table3,
                                    Level = 3, stringsAsFactors = FALSE)
dim(identification.table3)


identification.table <- dplyr::full_join(rbind(identification.table1, identification.table2),
                                         identification.table3,
                                         by = colnames(identification.table3))

write.csv(identification.table, "identification.table.csv", row.names = FALSE)

identification.table.new <- trans2newStyle(identification.table = identification.table)

peak_table <- readr::read_csv("Peak.table.csv")
temp <- left_join(x = peak_table, y = identification.table.new, by = "name")
temp <- temp %>%
  select(., -one_of(c("mz.x", "rt.x")))
colnames(temp)[2:3] <- c("mz", "rt")
identification.table.new <- temp
write.csv(temp, "identification.table.new.csv", row.names = FALSE)

level <- temp$Level
level[is.na(level)] <- "Unknown"
level <- data.frame(level, stringsAsFactors = FALSE)


plot <- ggplot(level, aes(x = level)) +
  geom_bar(width = 0.7, aes(fill = level)) +
  scale_fill_manual(values = c("#00468BB2", "#ED0000B2", "#42B540B2", "#0099B4B2"), guide = FALSE) +
  theme_bw() +
  labs(x = "Annotation level",  y = 'Count') +
  annotate(geom = "text", x = 1:4, y = table(level$level), label = table(level$level), size = 5) +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))


export::graph2ppt(plot, width = 8, height = 6)

ms1.data <- smart.d.rplc.pos5@ms1.data[[1]]
temp <-dplyr::left_join(x = ms1.data, y = identification.table.new)

ms1.data <- temp
ms1.data <- list(ms1.data)
smart.d.rplc.pos5@ms1.data <- ms1.data
save(smart.d.rplc.pos5, file = "smart.d.rplc.pos5")
#####RPLC negative
setwd("/home/shenxt/plasma/RPLC-NEG-MS2/CE50")
parameter1 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "msDatabase_rplc0.0.1",
                               threads = 10)

parameter2 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "hmdbDatabase0.0.1",
                               threads = 10)

parameter3 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "massbankDatabase0.0.1",
                               threads = 10)

parameter4 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "metlinDatabase0.0.1",
                               threads = 10)

parameter5 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "monaDatabase0.0.1",
                               threads = 10)

parameter6 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "nistDatabase0.0.1",
                               threads = 10)

parameter7 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "rp",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "orbitrapDatabase0.0.1",
                               threads = 10)

parameter8 <- mzIdentifyParam(ms1.match.ppm = 25,
                              polarity = "negative",
                              column = "rp",
                              candidate.num = 3,
                              database = "HMDB.metabolite.data",
                              threads = 10)
ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.nRPLC.nce50 <- metIdentify4all(ms1.data = ms1.data,
                                      ms2.data = ms2.data,
                                      parameter.list = c(parameter1,
                                                         parameter2,
                                                         parameter3,
                                                         parameter4,
                                                         parameter5,
                                                         parameter6,
                                                         parameter7,
                                                         parameter8),
                                      path = ".")

save(result.nRPLC.nce50, file = "result.nRPLC.nce50")
print(1)


#####HILIC positive
setwd("/home/shenxt/plasma/HILIC-POS-MS2/CE50")
library(metIdentify)
parameter1 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "msDatabase_hilic0.0.1",
                               threads = 10)

parameter2 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "hmdbDatabase0.0.1",
                               threads = 10)

parameter3 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "massbankDatabase0.0.1",
                               threads = 10)

parameter4 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "metlinDatabase0.0.1",
                               threads = 10)

parameter5 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "monaDatabase0.0.1",
                               threads = 10)

parameter6 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "nistDatabase0.0.1",
                               threads = 10)

parameter7 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "positive",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "orbitrapDatabase0.0.1",
                               threads = 10)

parameter8 <- mzIdentifyParam(ms1.match.ppm = 25,
                              polarity = "positive",
                              column = "hilic",
                              candidate.num = 3,
                              database = "HMDB.metabolite.data",
                              threads = 10)
ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.pHILIC.nce50 <- metIdentify4all(ms1.data = ms1.data,
                                       ms2.data = ms2.data,
                                       parameter.list = c(parameter1,
                                                          parameter2,
                                                          parameter3,
                                                          parameter4,
                                                          parameter5,
                                                          parameter6,
                                                          parameter7,
                                                          parameter8),
                                       path = ".")

save(result.pHILIC.nce50, file = "result.pHILIC.nce50")
print(1)














#####HILIC negative
library(metIdentify)
setwd("/home/shenxt/plasma/HILIC-NEG-MS2/CE25")
parameter1 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "msDatabase_hilic0.0.1",
                               threads = 10)

parameter2 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "hmdbDatabase0.0.1",
                               threads = 10)

parameter3 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "massbankDatabase0.0.1",
                               threads = 10)

parameter4 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "metlinDatabase0.0.1",
                               threads = 10)

parameter5 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "monaDatabase0.0.1",
                               threads = 10)

parameter6 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "nistDatabase0.0.1",
                               threads = 10)

parameter7 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60,
                               ms1.match.ppm = 25,
                               ms2.match.tol = 0.4,
                               rt.match.tol = 90,
                               polarity = "negative",
                               ce = "all",
                               column = "hilic",
                               total.score.tol = 0,
                               candidate.num = 3,
                               database = "orbitrapDatabase0.0.1",
                               threads = 10)

parameter8 <- mzIdentifyParam(ms1.match.ppm = 25,
                              polarity = "negative",
                              column = "hilic",
                              candidate.num = 3,
                              database = "HMDB.metabolite.data",
                              threads = 10)
ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.nHILIC.nce25 <- metIdentify4all(ms1.data = ms1.data,
                                       ms2.data = ms2.data,
                                       parameter.list = c(parameter1,
                                                          parameter2,
                                                          parameter3,
                                                          parameter4,
                                                          parameter5,
                                                          parameter6,
                                                          parameter7,
                                                          parameter8),
                                       path = ".")

save(result.nHILIC.nce25, file = "result.nHILIC.nce25")
print(1)



####
setwd("E:/project/NGLY1/NGLY1-plasma-id-20190709/RPLC-NEG-MS2")
load("CE25/result.nRPLC.nce25")
load("CE50/result.nRPLC.nce50")


load("smart.d.rplc.pos5")
smart.d.rplc.pos5

sample.info <- smart.d.rplc.pos5@sample.info


temp <- 
  patient_info %>% 
  select(Sample_ID, `Visit GA`)

temp$Sample_ID <- 
  paste("X",stringr::str_extract(temp$'Sample_ID', "[0-9]{1,4}"), sep = "")

colnames(temp)[1] <- "sample.name"
sample.info <-
  left_join(sample.info, temp)


match(sample.info$sample.name, temp$Sample_ID)

subject_data <- getData(object = smart.d.rplc.pos5, slot = "Subject")
sample.info <- 
  sample.info %>% 
  filter(!is.na(`Visit GA`))

subject_data <- 
subject_data %>% 
  select(one_of(sample.info$sample.name))


subject_data <- 
  t(
    apply(subject_data, 1, function(x){
  (x - mean(x))/sd(x)
})
)


subject_data <-
  subject_data %>% 
  as_tibble()


subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample.info$`Visit GA`)

pca_object <- prcomp(x = select(subject_data2, -GA))


library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, GA = subject_data2$GA, stringsAsFactors = FALSE)

ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_point() +
  scale_colour_gradientn(colours = c('skyblue', "black", "red")) +
  theme_bw()


(plot <- ggplot2::autoplot(pca_object, data = subject_data2,
                           colour = 'GA',
                          # colour = 'class',
                          frame = FALSE, frame.type = "norm") +
  # scale_y_continuous(limits = c(0, 0.015)) +
  # scale_x_continuous(limits = c(-0.1, 0.1)) +
  theme_bw() +
  scale_colour_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))
)

plot
