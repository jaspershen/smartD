###negative
setwd("E:/project/smartD/RPLC/NEG/data cleaning")
load("smart.d.rplc.neg5")
smart.d.rplc.neg5

sample.info <- smart.d.rplc.neg5@sample.info

temp <- 
  patient_info %>% 
  select(Sample_ID, `Visit GA`)

temp$Sample_ID <- 
  paste("X",stringr::str_extract(temp$'Sample_ID', "[0-9]{1,4}"), sep = "")

colnames(temp)[1] <- "sample.name"
sample.info <-
  left_join(sample.info, temp)

subject_data <- getData(object = smart.d.rplc.neg5, slot = "Subject")
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

##PCA analysis
pca_object <- prcomp(x = select(subject_data2, -GA))

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, GA = subject_data2$GA, stringsAsFactors = FALSE)

ggplot(x, aes(PC1, PC2, colour = GA)) +
  geom_point() +
  scale_colour_gradient(low = alpha("red", 0.1), 
                         # mid = alpha("red", 0.7), 
                         high = alpha("red", 1)) +
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
  data.frame(Y, "GA" = sample.info$`Visit GA`, stringsAsFactors = FALSE)

(
plot <- ggplot(Y, aes(X1, X2, colour = GA)) +
  geom_point() +
  labs(x = "Dimension 1",
       y = "Dimension 2") +
  theme_bw() +
    scale_colour_gradient(low = alpha("red", 0.1), 
                          # mid = alpha("red", 0.7), 
                          high = alpha("red", 1)) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15)
  )
)






###positive
setwd("E:/project/smartD/RPLC/POS/data cleaning")
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
  data.frame(Y, "GA" = sample.info$`Visit GA`, stringsAsFactors = FALSE)

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









