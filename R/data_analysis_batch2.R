#-------------------------------------------------------------------------------
##batch 2 RPLC pos
#-------------------------------------------------------------------------------
setwd("E:/project/smartD/smartD_batch2/RPLC/POS/data_analysis")
load("smartd_rplc_pos_batch2_5")
sample.info <- smartd_rplc_pos_batch2_5@sample.info

subject_data <- getData(object = smartd_rplc_pos_batch2_5, slot = "Subject")
sample.info <- 
  sample.info %>% 
  filter(GA != 0)

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
  mutate(GA = sample.info$GA)

##PCA analysis
pca_object <- prcomp(x = select(subject_data2, -GA))

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, GA = subject_data2$GA, stringsAsFactors = FALSE)

ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("red", 0.1), 
  #                       high = alpha("red", 1)) +
  scale_colour_gradientn(colours = c(
    alpha("blue", 1),
    alpha("blue", 0.1),
    alpha("red", 0.1),
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

###tsne analysis
tsne_object <- Rtsne::Rtsne(
  X = t(as.matrix(subject_data)),
  dims = 2,
  perplexity = 30,
  verbose = TRUE
)

Y <- tsne_object$Y
Y <-
  data.frame(Y, "GA" = sample.info$`GA`, stringsAsFactors = FALSE)

(
  plot <- ggplot(Y, aes(X1, X2, colour = GA)) +
    geom_point(size = 3) +
    labs(x = "Dimension 1",
         y = "Dimension 2") +
    theme_bw() +
    scale_colour_gradientn(colours = c(
      alpha("blue", 1),
      alpha("blue", 0.1),
      alpha("red", 0.1),
      alpha("red", 1)
    )) +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      strip.background = element_rect(fill = "#0099B47F"),
      strip.text = element_text(color = "white", size = 15)
    )
)







#-------------------------------------------------------------------------------
##batch 2 RPLC neg
#-------------------------------------------------------------------------------
setwd("E:/project/smartD/smartD_batch2/RPLC/NEG/data_analysis")
load("smartd_rplc_neg_batch2_5")

subject_data <- 
  getData(object = smartd_rplc_neg_batch2_5, slot = "Subject")

sample.info <-
  smartd_rplc_neg_batch2_5@sample.info

sample.info <- 
  sample.info %>% 
  filter(GA != 0)

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
  mutate(GA = sample.info$GA)

rownames(subject_data2) <-
  sample.info$sample.name

##PCA analysis
pca_object <- prcomp(x = select(subject_data2, -GA))

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, GA = subject_data2$GA, stringsAsFactors = FALSE)

plot <- 
ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("red", 0.1), 
  #                       high = alpha("red", 1)) +
  scale_colour_gradientn(colours = c(
    alpha("blue", 1),
    alpha("blue", 0.1),
    alpha("red", 0.1),
    alpha("red", 1)
  )) +
  # xlim(-20,20) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))


plot
x2 <- 
  x %>% 
  rownames_to_column(var = "sample.name")
plot +
  ggrepel::geom_text_repel(aes(x = PC1, y = PC2, 
                               label = sample.name), 
                           x2, colour = "black")

##remove some outliers
sample.info <- 
  sample.info %>% 
  filter(!sample.name %in% c("X65", "X82", "X68", "X66", "X75", "X67", "X85", "X62",
                             "X77", "X84", "X76", "X83", "X81", "X72", "X79", "X89",
                             "X89", "X78", "X70", "X71", "X80", "X73", "X74"))

subject_data <- 
  subject_data %>% 
  select(one_of(sample.info$sample.name))



###tsne analysis
tsne_object <- Rtsne::Rtsne(
  X = t(as.matrix(subject_data)),
  dims = 2,
  perplexity = 30,
  verbose = TRUE
)

Y <- tsne_object$Y
Y <-
  data.frame(Y, "GA" = sample.info$`GA`, stringsAsFactors = FALSE)

(
  plot <- ggplot(Y, aes(X1, X2, colour = GA)) +
    geom_point(size = 3) +
    labs(x = "Dimension 1",
         y = "Dimension 2") +
    theme_bw() +
    scale_colour_gradientn(colours = c(
      alpha("blue", 1),
      alpha("blue", 0.1),
      alpha("red", 0.1),
      alpha("red", 1)
    )) +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      strip.background = element_rect(fill = "#0099B47F"),
      strip.text = element_text(color = "white", size = 15)
    )
)





#-------------------------------------------------------------------------------
##batch 2 HILIC pos
#-------------------------------------------------------------------------------
setwd("E:/project/smartD/smartD_batch2/HILIC/POS/data_analysis")
load("smartd_hilic_pos_batch2_5")
sample.info <- smartd_hilic_pos_batch2_5@sample.info

subject_data <- getData(object = smartd_hilic_pos_batch2_5, 
                        slot = "Subject")
sample.info <- 
  sample.info %>% 
  filter(GA != 0)

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
  mutate(GA = sample.info$GA)

##PCA analysis
pca_object <- prcomp(x = select(subject_data2, -GA))

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, GA = subject_data2$GA, stringsAsFactors = FALSE)

ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_point() +
  # scale_colour_gradient(low = alpha("red", 0.1), 
  #                       high = alpha("red", 1)) +
  scale_colour_gradientn(colours = c(
    alpha("blue", 1),
    alpha("blue", 0.1),
    alpha("red", 0.1),
    alpha("red", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  ylim(-10, 20) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))

###tsne analysis
tsne_object <- Rtsne::Rtsne(
  X = t(as.matrix(subject_data)),
  dims = 2,
  perplexity = 30,
  verbose = TRUE
)

Y <- tsne_object$Y
Y <-
  data.frame(Y, "GA" = sample.info$`GA`, stringsAsFactors = FALSE)

(
  plot <- ggplot(Y, aes(X1, X2, colour = GA)) +
    geom_point(size = 3) +
    labs(x = "Dimension 1",
         y = "Dimension 2") +
    theme_bw() +
    scale_colour_gradientn(colours = c(
      alpha("blue", 1),
      alpha("blue", 0.1),
      alpha("red", 0.1),
      alpha("red", 1)
    )) +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      strip.background = element_rect(fill = "#0099B47F"),
      strip.text = element_text(color = "white", size = 15)
    )
)


#-------------------------------------------------------------------------------
##batch 2 HILIC neg
#-------------------------------------------------------------------------------
setwd("E:/project/smartD/smartD_batch2/HILIC/NEG/data_analysis")
load("smartd_hilic_neg_batch2_5")
sample.info <- smartd_hilic_neg_batch2_5@sample.info

subject_data <- getData(object = smartd_hilic_neg_batch2_5, slot = "Subject")
sample.info <- 
  sample.info %>% 
  filter(GA != 0)

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
  mutate(GA = sample.info$GA)

##PCA analysis
pca_object <- prcomp(x = select(subject_data2, -GA))

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, GA = subject_data2$GA, stringsAsFactors = FALSE)

ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_point() +
  # scale_colour_gradient(low = alpha("red", 0.1), 
  #                       high = alpha("red", 1)) +
  scale_colour_gradientn(colours = c(
    alpha("blue", 1),
    alpha("blue", 0.1),
    alpha("red", 0.1),
    alpha("red", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  # ylim(-10, 20) +
  xlim(0, 10) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))

###tsne analysis
tsne_object <- Rtsne::Rtsne(
  X = t(as.matrix(subject_data)),
  dims = 2,
  perplexity = 30,
  verbose = TRUE
)

Y <- tsne_object$Y
Y <-
  data.frame(Y, "GA" = sample.info$`GA`, stringsAsFactors = FALSE)

(
  plot <- ggplot(Y, aes(X1, X2, colour = GA)) +
    geom_point(size = 3) +
    labs(x = "Dimension 1",
         y = "Dimension 2") +
    theme_bw() +
    scale_colour_gradientn(colours = c(
      alpha("blue", 1),
      alpha("blue", 0.1),
      alpha("red", 0.1),
      alpha("red", 1)
    )) +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      strip.background = element_rect(fill = "#0099B47F"),
      strip.text = element_text(color = "white", size = 15)
    )
)





