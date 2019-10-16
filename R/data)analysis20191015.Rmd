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
  ggplot(x, aes(PC1, PC2, colour = batch)) +
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