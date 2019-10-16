### PCA analysis to check the data quality
##load dataset

setwd("data_analysis20191015/all_feature/")
load("smartd_rplc_pos_5")
load("smartd_rplc_neg_5")

##Let us see if the samples distrubite in PCA score plot.
#### RPLC pos batch 1&2
library(metflow2)
library(ggplot2)
library(tidyverse)

sample.info <- 
  smartd_rplc_pos_5@sample.info

subject_data <- 
  getData(object = smartd_rplc_pos_5, slot = "Subject")

sample.info <- 
  sample.info %>% 
  filter(GA != 0)

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))

##auto scale
subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
# rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample.info$GA, batch = sample.info$batch)

##PCA analysis
pca_object <- 
  prcomp(x = select(subject_data2, -c(GA,batch)))

library(ggfortify)
x <- pca_object$x
x <- x[,1:2]
x <- 
  data.frame(x, GA = subject_data2$GA, 
             batch = subject_data2$batch,
             stringsAsFactors = FALSE)
x <- data.frame(x, name = colnames(subject_data), 
                stringsAsFactors = FALSE)

pca_score_plot <- 
  ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 3) +
  scale_colour_gradientn(guide = guide_colorbar(title = "GA (week)"),
                         colours = c(
                           alpha("royalblue", 1),
                           alpha("royalblue", 0.3),
                           alpha("red", 0.3),
                           alpha("red", 1)
                         )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))
pca_score_plot
ggsave(pca_score_plot, filename = "PCA_score_plot_batch1_2_RPLC_POS.pdf", 
       width = 8, height = 6)


##batch1
sample.info <- 
  smartd_rplc_pos_5@sample.info

subject_data <- 
  getData(object = smartd_rplc_pos_5, slot = "Subject")

sample.info <- 
  sample.info %>% 
  filter(GA != 0)

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))

##get batch 1 data
sample.info <- 
  sample.info %>% 
  filter(batch == 1)

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))

##auto scale
subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
# rownames(subject_data2)
subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample.info$GA, batch = sample.info$batch)

##PCA analysis
pca_object <- 
  prcomp(x = select(subject_data2, -c(GA,batch)))

library(ggfortify)
x <- pca_object$x
x <- x[,1:2]
x <- 
  data.frame(x, GA = subject_data2$GA, 
             batch = subject_data2$batch,
             stringsAsFactors = FALSE)
x <- data.frame(x, name = colnames(subject_data), 
                stringsAsFactors = FALSE)

pca_score_plot <- 
  ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 3) +
  scale_colour_gradientn(guide = guide_colorbar(title = "GA (week)"),
                         colours = c(
                           alpha("royalblue", 1),
                           alpha("royalblue", 0.3),
                           alpha("red", 0.3),
                           alpha("red", 1)
                         )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))
pca_score_plot
ggsave(pca_score_plot, filename = "PCA_score_plot_batch1_RPLC_POS.pdf", 
       width = 8, height = 6)



##batch2
sample.info <- 
  smartd_rplc_pos_5@sample.info

subject_data <- 
  getData(object = smartd_rplc_pos_5, slot = "Subject")

sample.info <- 
  sample.info %>% 
  filter(GA != 0)

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))

##get batch 2 data
sample.info <- 
  sample.info %>% 
  filter(batch == 2)

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))

##auto scale
subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
# rownames(subject_data2)
subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample.info$GA, batch = sample.info$batch)

##PCA analysis
pca_object <- 
  prcomp(x = select(subject_data2, -c(GA,batch)))

library(ggfortify)
x <- pca_object$x
x <- x[,1:2]
x <- 
  data.frame(x, GA = subject_data2$GA, 
             batch = subject_data2$batch,
             stringsAsFactors = FALSE)
x <- data.frame(x, name = colnames(subject_data), 
                stringsAsFactors = FALSE)

pca_score_plot <- 
  ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 3) +
  scale_colour_gradientn(guide = guide_colorbar(title = "GA (week)"),
                         colours = c(
                           alpha("royalblue", 1),
                           alpha("royalblue", 0.3),
                           alpha("red", 0.3),
                           alpha("red", 1)
                         )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))
pca_score_plot
ggsave(pca_score_plot, filename = "PCA_score_plot_batch2_RPLC_POS.pdf", 
       width = 8, height = 6)









#### RPLC neg batch 1&2
library(metflow2)
library(ggplot2)
library(tidyverse)

sample.info <- 
  smartd_rplc_neg_5@sample.info

subject_data <- 
  getData(object = smartd_rplc_neg_5, slot = "Subject")

sample.info <- 
  sample.info %>% 
  filter(GA != 0)

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))

##auto scale
subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
# rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample.info$GA, batch = sample.info$batch)

##PCA analysis
pca_object <- 
  prcomp(x = select(subject_data2, -c(GA,batch)))

library(ggfortify)
x <- pca_object$x
x <- x[,1:2]
x <- 
  data.frame(x, GA = subject_data2$GA, 
             batch = subject_data2$batch,
             stringsAsFactors = FALSE)
x <- data.frame(x, name = colnames(subject_data), 
                stringsAsFactors = FALSE)

pca_score_plot <- 
  ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 3) +
  scale_colour_gradientn(guide = guide_colorbar(title = "GA (week)"),
                         colours = c(
                           alpha("royalblue", 1),
                           alpha("royalblue", 0.3),
                           alpha("red", 0.3),
                           alpha("red", 1)
                         )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))
pca_score_plot
ggsave(pca_score_plot, filename = "PCA_score_plot_batch1_2_RPLC_NEG.pdf", 
       width = 8, height = 6)


##batch1
sample.info <- 
  smartd_rplc_neg_5@sample.info

subject_data <- 
  getData(object = smartd_rplc_neg_5, slot = "Subject")

sample.info <- 
  sample.info %>% 
  filter(GA != 0)

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))

##get batch 1 data
sample.info <- 
  sample.info %>% 
  filter(batch == 1)

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))

##auto scale
subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
# rownames(subject_data2)
subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample.info$GA, batch = sample.info$batch)

##PCA analysis
pca_object <- 
  prcomp(x = select(subject_data2, -c(GA,batch)))

library(ggfortify)
x <- pca_object$x
x <- x[,1:2]
x <- 
  data.frame(x, GA = subject_data2$GA, 
             batch = subject_data2$batch,
             stringsAsFactors = FALSE)
x <- data.frame(x, name = colnames(subject_data), 
                stringsAsFactors = FALSE)

pca_score_plot <- 
  ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 3) +
  scale_colour_gradientn(guide = guide_colorbar(title = "GA (week)"),
                         colours = c(
                           alpha("royalblue", 1),
                           alpha("royalblue", 0.3),
                           alpha("red", 0.3),
                           alpha("red", 1)
                         )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))
pca_score_plot
ggsave(pca_score_plot, filename = "PCA_score_plot_batch1_RPLC_NEG.pdf", 
       width = 8, height = 6)



##batch2
sample.info <- 
  smartd_rplc_neg_5@sample.info

subject_data <- 
  getData(object = smartd_rplc_neg_5, slot = "Subject")

sample.info <- 
  sample.info %>% 
  filter(GA != 0)

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))

##get batch 2 data
sample.info <- 
  sample.info %>% 
  filter(batch == 2)

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))

##auto scale
subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
# rownames(subject_data2)
subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample.info$GA, batch = sample.info$batch)

##PCA analysis
pca_object <- 
  prcomp(x = select(subject_data2, -c(GA,batch)))

library(ggfortify)
x <- pca_object$x
x <- x[,1:2]
x <- 
  data.frame(x, GA = subject_data2$GA, 
             batch = subject_data2$batch,
             stringsAsFactors = FALSE)
x <- data.frame(x, name = colnames(subject_data), 
                stringsAsFactors = FALSE)

pca_score_plot <- 
  ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 3) +
  scale_colour_gradientn(guide = guide_colorbar(title = "GA (week)"),
                         colours = c(
                           alpha("royalblue", 1),
                           alpha("royalblue", 0.3),
                           alpha("red", 0.3),
                           alpha("red", 1)
                         )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))
pca_score_plot
ggsave(pca_score_plot, filename = "PCA_score_plot_batch2_RPLC_NEG.pdf", 
       width = 8, height = 6)



### Patient information
##All the patient information are in 'E:/project/smartD/data_analysis20190828/patient_info'
##This is the file contains all the information of all the samples
all346_urine_info <- 
  readr::read_csv("E:/project/smartD/data_analysis20190828/patient_info/SmartD_all346urine.csv")

##PTID is the ID of patients
all346_urine_info %>% 
  group_by(PTID) %>% 
  dplyr::summarise(n = n())
dim(all346_urine_info)

##Ther are 36 patients and 346 samples in total. 

all346_urine_info %>% 
  group_by(PTID) %>% 
  dplyr::summarise(n = n()) %>% 
  arrange(n) %>% 
  mutate(PTID = factor(PTID, levels = PTID)) %>% 
  ggplot(.,aes(x = PTID, y = n)) +
  geom_point(colour = "red") +
  geom_segment(aes(x = 1:36, y = 0, xend = 1:36, yend = n), colour = "grey") +
  coord_flip() +
  labs(x = "Patient ID", y = "Sample number") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

#You can see that different people have different number of samples. 
#Now we want to get the distribution of sample collection.
temp_data <- 
  all346_urine_info %>% 
  select(PTID, Visit.GA)

plot1 <- 
  ggplot(temp_data, aes(x = Visit.GA, y = PTID)) +
  geom_point() +
  labs(x = "GA (weeks)", y = "Patient ID") +
  theme_bw() +
  scale_x_continuous(limits = c(10, 43)) +
  theme(panel.grid.major.y = element_line(colour = "#8A9045FF"),
        plot.margin = margin(t = 0,r = 0, b = 0, l = 0))

plot2 <-
  ggplot(temp_data, aes(x = Visit.GA)) +
  geom_bar(binwidth = 0.5) +
  theme_bw() +
  scale_x_continuous(limits = c(10, 43), 
                     name = NULL, labels = NULL, breaks = NULL) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(0,0,0,0))


plot3 <- 
  ggplot(temp_data, aes(x = PTID)) +
  geom_bar(width = 0.8) +
  theme_bw() +
  scale_x_discrete(name = NULL, label = NULL, breaks = NULL) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
  coord_flip() +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(0,0,0,0))

library(patchwork)

plot <- 
  {plot2 + plot_spacer() + plot_layout(ncol = 2, widths = c(3, 1))} -
  {plot1 + plot3 + plot_layout(ncol = 2, widths = c(3, 1))} +
  plot_layout(ncol = 1, heights = c(1,3))

plot

ggsave(plot, 
       filename = "Sample_colection_distribution.pdf", 
       width = 8, height = 6)

plot




##we use Batch 1 as discovery dataset and batch 2 as validation dataset
##So we should get the demographics characterisitics of discovery and validation datasets.

