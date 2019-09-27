###20190918 use all the identified metabolites (Level 1 and Level 2)
setwd("E:/project/smartD/data_analysis20190828/correlation_network_analysis")
setwd("identification_table/")

identification_table_pos <- readr::read_csv("POS/identification.table.new.csv")
identification_table_neg <- readr::read_csv("NEG/identification.table.new.csv")

load("POS/smartd_rplc_pos_batch1_5")
ms1_data_pos <- smartd_rplc_pos_batch1_5@ms1.data[[1]]

load("NEG/smartd_rplc_neg_batch1_5")
ms1_data_neg <- smartd_rplc_neg_batch1_5@ms1.data[[1]]

library(tidyverse)

###we only use the level 1 and level 2 identifications
identification_table_pos <- 
  identification_table_pos %>% 
  filter(Level == 1 | Level == 2)

identification_table_neg <- 
  identification_table_neg %>% 
  filter(Level == 1 | Level == 2)

identification_table_pos <- 
  inner_join(identification_table_pos, ms1_data_pos, by = c("name", "mz", "rt"))

identification_table_neg <- 
  inner_join(identification_table_neg, ms1_data_neg, by = c("name", "mz", "rt"))

identification_table_pos <- 
  identification_table_pos %>% 
  mutate(., name = paste(name, "POS", sep = "_"))

identification_table_neg <- 
  identification_table_neg %>% 
  mutate(., name = paste(name, "NEG", sep = "_"))



write.csv(identification_table_pos, file = "./POS_NEG/identification_table_pos.csv", row.names = FALSE)

write.csv(identification_table_neg, file = "./POS_NEG/identification_table_neg.csv", row.names = FALSE)

### combine pos and neg table together
setwd("E:/project/smartD/data_analysis20190828/correlation_network_analysis/identification_table/POS_NEG")
identification_table_pos <- readr::read_csv("identification_table_pos.csv")
identification_table_neg <- readr::read_csv("identification_table_neg.csv")

##some samples are not different in pos and neg
smartd_rplc_pos_batch1 <- identification_table_pos
smartd_rplc_neg_batch1 <- identification_table_neg

smartd_rplc_batch1 <- full_join(smartd_rplc_pos_batch1, 
                                smartd_rplc_neg_batch1, 
                                by = intersect(colnames(smartd_rplc_pos_batch1), 
                                               colnames(smartd_rplc_neg_batch1)))

##remove QC_DL and Blank Samples
smartd_rplc_batch1 <- 
  smartd_rplc_batch1 %>% 
  select(., -contains("BLK")) %>% 
  select(., -contains("QC")) %>% 
  select(name:rt, MS2.spectrum.name: Database, everything())

#####smartd_rplc_batch dataset is used to analysis for correlation networks analysis
#### step 1 remove the duplicated metabolites according to score
smartd_rplc_batch1 <- 
  smartd_rplc_batch1 %>% 
  group_by(., Compound.name) %>% 
  filter(., SS == max(SS)) %>% 
  ungroup()




################cross-sectional correlation network
cross_temp <- 
  smartd_rplc_batch1 %>% 
  select(-(name:Database))

table(which(is.na(cross_temp), arr.ind = TRUE)[,2])
colnames(cross_temp)[c(146, 147, 148)]
sum(is.na(cross_temp))

##SFU65, SFU74, SFU91 have no positive data, so remove it
cross_temp <- 
  cross_temp %>% 
  select(., -c(SFU65, SFU74, SFU91))

write.csv(smartd_rplc_batch1, file = "smartd_rplc_batch1.csv", row.names = FALSE)

###calculate correlation matrix
###cross sectional network
setwd("E:/project/smartD/data_analysis20190828/correlation_network_analysis/identification_table/POS_NEG/cross_sectional_netowork")
sfu1_148 <- 
  readr::read_csv("E:/project/smartD/patient information/SFU1-148_GA.csv")

patient_info <- 
  readr::read_csv("E:/project/smartD/data_analysis20190828/patient_info/patient_info.csv")

sfu1_148 <- 
  sfu1_148 %>% 
  mutate(subject_id = as.character(subject_id), visit = as.character(visit)) %>% 
  arrange(subject_id)

patient_info <- 
  patient_info %>% 
  mutate(Patient_ID = as.character(Patient_ID), Visit = as.character(Visit)) %>% 
  arrange(Patient_ID)

match(colnames(cross_temp), sfu1_148$sample_id)

cross_temp <- 
  log(cross_temp, 10) %>% 
  as.data.frame(cross_temp)

####log 10 and scale
cross_temp <-
  apply(cross_temp, 1, function(x){
    (x - mean(x))/sd(x)
  })

cross_temp <- 
  cross_temp %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_name")

colnames(cross_temp)[-1] <- 
  smartd_rplc_batch1$name

cross_temp <- 
  inner_join(x = sfu1_148[,-1], cross_temp, 
             by = c("sample_id" = "Sample_name"))


###calculate correlation

cross_temp2 <- 
  cross_temp %>% 
  plyr::dlply(., .(subject_id)) %>% 
  lapply(., function(x){
    x <- 
      x %>% 
      select(., -c(subject_id:GA_week)) %>% 
      apply(., 2, mean)
  }) %>% 
  do.call(rbind, .) %>% 
  as_tibble()
  
cross_correlation_matrix <- 
  cross_temp2 %>% 
  as.matrix(.) %>% 
  cor(., method = "spearman")


cross_correlation_data <- 
  cross_correlation_matrix 

rownames(cross_correlation_data) <- 
  colnames(cross_correlation_data)

##set lower tri as o
cross_correlation_data[lower.tri(cross_correlation_data)] <- 2

cross_correlation_data <- 
  cross_correlation_data %>% 
  as_tibble()

rownames(cross_correlation_data) <- 
  colnames(cross_correlation_data)

cross_correlation_data <- 
  cross_correlation_data %>% 
  rownames_to_column(., var = "Peak.name1") %>% 
  gather(., key = "Peak.name2", value = "Correlation", -Peak.name1) %>% 
  distinct() %>% 
  filter(., Correlation != 1 & Correlation != 2) %>% 
  arrange(., desc(abs(Correlation))) %>% 
  filter(abs(Correlation) > 0.6)

cross_correlation_p <-
  pbapply::pbapply(cross_correlation_data, 1, function(x) {
    peak1 <- as.character(x[1])
    peak2 <- as.character(x[2])
    int1 <-
      pull(cross_temp2, var = peak1)
    int2 <-
      pull(cross_temp2, var = peak2)
    p <-  cor.test(int1, int2, alternative = "two.sided",
                   method = "spearman")$p.value
    p
  })


###p adjustment
cross_correlation_p2 <- 
  p.adjust(cross_correlation_p, method = "BH")

cross_correlation_data <- 
  cross_correlation_data %>% 
  mutate(P.adjust = cross_correlation_p2) %>% 
  filter(P.adjust < 0.05 & abs(Correlation) > 0.6)

###add compound name
cross_correlation_data <- 
  cross_correlation_data %>% 
  left_join(., smartd_rplc_batch1[,c("name", "Compound.name")], 
            by = c("Peak.name1" = "name")) %>% 
  dplyr::rename(Compound.name1 = Compound.name) %>% 
  left_join(., smartd_rplc_batch1[,c("name", "Compound.name")], 
            by = c("Peak.name2" = "name")) %>% 
  dplyr::rename(Compound.name2 = Compound.name)


######circos figure
library(circlize)
cross_correlation_data2 <- 
  cross_correlation_data[1:100,]

ggplot(data = cross_correlation_data2, aes(x = Peak.name1, y = Peak.name2)) +
  geom_tile(aes(fill = Correlation)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw()

Metabolite <- 
  unique(c(cross_correlation_data2$Compound.name1, cross_correlation_data2$Compound.name2))
cross_temp_data <- data.frame(factor = "Metabolomics", 
                        x = 1:length(Metabolite),
                        y = 1,
                        Metabolite = Metabolite)


library(circlize)
par(mar = c(2,2,2,2))
circos.par(start.degree = 75, clock.wise = TRUE, gap.after = 6.3)
circos.initialize(factors = cross_temp_data$factor, x = cross_temp_data$x)
circos.track(factors = cross_temp_data$factor, 
             x = cross_temp_data$x, 
             y = cross_temp_data$y - 1,

             ylim = c(0, 1),
             bg.border = NA,
             bg.col = NA,
             track.height = 0.4,
             # track.margin = c(0,0,0,0),
             panel.fun = function(x, y){
               circos.points(x = x, y = y, pch = 19, col = "black")
               circos.text(x = x, y = y + uy(2, "mm"), 
                           niceFacing = TRUE, 
                           labels = cross_temp_data$Metabolite,
                           facing = "clockwise", 
                           adj = c(0, 0.5),
                           col = "black", cex = 0.8)
  
})

col_fun = colorRamp2(c(-1, 0, 1), c("#3C5488FF", "white", "#E64B35FF"))

for(i in 1:nrow(cross_correlation_data2)){
  cat(i, " ")
  point1 <- 
    match(cross_correlation_data[i,]$Compound.name1, cross_temp_data$Metabolite)
  point2 <- 
    match(cross_correlation_data[i,]$Compound.name2, cross_temp_data$Metabolite)
  correlation <- cross_correlation_data[i,]$Correlation
    col <- col_fun(x = correlation)
  circos.link(sector.index1 = "Metabolomics", point1 = point1, 
              sector.index2 = "Metabolomics", point2 = point2, 
              h = 1, col = col, lwd = 2)
}

circos.clear()



######network community analysis
##https://github.com/IOR-Bioinformatics/PCSF/
##modularity定义是指网络中连接社区结构内部顶点的边所占的比例与另外一个随机网络中丽娜姐社区结构内部顶部
#的边所占比例的期望值相减得到的差值.
library(igraph)

###construct igraph project of correlation data
cross_correlation_data3 <- 
  cross_correlation_data %>% 
  select(., from = Peak.name1, to = Peak.name2, 
         Correlation, P.adjust) %>% 
  as.data.frame()

node.attr <- smartd_rplc_batch1[,c(1:18)] %>% 
  filter(., name %in% (c(cross_correlation_data3$from, cross_correlation_data3$to)))

cross_correlation_data_igraph <- 
  graph_from_data_frame(d = cross_correlation_data3, 
                        directed = FALSE, 
                        vertices = node.attr
                        )
  
# coords <- layout.davidson.harel(cross_correlation_data_igraph)
lay <- layout.auto(graph = cross_correlation_data_igraph)
plot(cross_correlation_data_igraph, 
     vertex.shape = "circle",
     vertex.label = NA,
     vertex.color = "#FFA319FF",
     vertex.size = vertex_attr(cross_correlation_data_igraph, "SS") * 10,
     # edge.color = "#767676FF",
     edge.color = ifelse(edge_attr(cross_correlation_data_igraph, "Correlation") > 0,
                         "#80000099", "#155F8399"),
     edge.width = abs(edge_attr(cross_correlation_data_igraph, "Correlation")) * 5,
     layout = lay
     )

##legend
# size_legend <-  
# quantile(x = vertex_attr(cross_correlation_data_igraph, "SS") * 10 %>% 
#              sort(), 
#            probs = c(0, 0.5, 1)) %>% 
#   round()
#   
# legend("topleft", 
#        legend = size_legend, 
#        bty = "n", 
#        pt.cex = size_legend/3,
#        pch = 19, 
#        # col = "black",
#        col = "#FFA319FF")



cross_subnetwork <- 
  cluster_edge_betweenness(graph = cross_correlation_data_igraph, 
                           weights = abs(edge_attr(cross_correlation_data_igraph, "Correlation")))

par(mar = c(5,5,4,2))

ggplot(
  data.frame(index = 1:length(cross_subnetwork$modularity),
             modu = cross_subnetwork$modularity, stringsAsFactors = FALSE),
  aes(index, modu) 
) +
  geom_vline(xintercept = which.max(cross_subnetwork$modularity), 
             linetype = 2, colour = "#800000B2") + 
  labs(x = "Community analysis iteration", y = "Modularity") +
  geom_line() +
  # geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))


membership(cross_subnetwork)
modularity(cross_subnetwork)
length(cross_subnetwork)
sizes(cross_subnetwork)
algorithm(cross_subnetwork)
merges(cross_subnetwork)
crossing(communities = cross_subnetwork, graph = cross_correlation_data_igraph)
code_len(cross_subnetwork)
# is_hierarchical(cross_subnetwork)
# plot(as.dendrogram(cross_subnetwork))
# plot(as.hclust(cross_subnetwork))

which(sort(table(membership(communities = cross_subnetwork)), decreasing = TRUE) > 2)
which(sort(table(membership(communities = cross_subnetwork)), decreasing = TRUE) > 5)




temp_id <-
  sizes(cross_subnetwork) %>% 
  as_tibble() %>% 
  dplyr::rename(., subnetwork_ID = `Community sizes`, size = n) %>% 
  arrange(., desc(size)) %>% 
  filter(., size > 2) %>% 
  pull(., subnetwork_ID)

temp_subnetwork <-
lapply(temp_id, function(x){
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == x])
})

###node number
node_number <- 
lapply(temp_subnetwork, function(x){
vertex_attr(graph = x, name = "name") %>% 
    length()
}) %>% 
  unlist() 
  # sum()


var <- temp_id  # the categorical data 
## Prep data (nothing to change here)
nrows <- 10
df <- expand.grid(y = 1:nrows, x = 1:nrows)
categ_table <- round(table(var) * ((nrows*nrows)/(length(var))))
categ_table

df$category <- factor(rep(names(categ_table), categ_table))  
# NOTE: if sum(categ_table) is not 100 (i.e. nrows^2), it will need adjustment to make the sum to 100.

## Plot
ggplot(df, aes(x = x, y = y, fill = category)) + 
        geom_tile(color = "black", size = 0.5) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0), trans = 'reverse') +
        scale_fill_brewer(palette = "Set3") +
        labs(title="Waffle Chart", subtitle="'Class' of vehicles",
             caption="Source: mpg") + 
        theme(
          # panel.border = element_rect(size = 1),
              plot.title = element_text(size = rel(1.2)),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              legend.title = element_blank(),
              legend.position = "right")

###edge number
lapply(temp_subnetwork, function(x){
as_edgelist(x) %>% nrow()
}) %>% 
  unlist() %>% 
  mean()





sizes(cross_subnetwork) %>% 
  as_tibble() %>% 
  dplyr::rename(., subnetwork_ID = `Community sizes`, size = n) %>% 
  arrange(., desc(size)) %>% 
  filter(., size > 2) 


sizes(cross_subnetwork) %>% 
  as_tibble() %>% 
  dplyr::rename(., subnetwork_ID = `Community sizes`, size = n) %>% 
  arrange(., desc(size)) %>% 
  filter(., size > 5)

###there are 15 cross_subnetworks that with node bigger than 5.
##there are "5"  "26" "53" "9"  "33" "37" "43" "41" "55" "71" "29" "31" "49" "84" "96"

cross_subnetwork5 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 5])

cross_subnetwork26 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 26])

cross_subnetwork53 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 53])

cross_subnetwork9 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 9])

cross_subnetwork33 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 33])

cross_subnetwork33 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 33])

cross_subnetwork37 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 37])

cross_subnetwork43 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 43])

cross_subnetwork41 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 41])

cross_subnetwork55 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 55])

cross_subnetwork71 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 71])

cross_subnetwork29 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 29])

cross_subnetwork31 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 31])

cross_subnetwork49 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 49])

cross_subnetwork84 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 84])

cross_subnetwork96 <- 
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == 96])

save(cross_subnetwork5, file = "cross_subnetwork5")
save(cross_subnetwork26, file = "cross_subnetwork26")
save(cross_subnetwork53, file = "cross_subnetwork53")
save(cross_subnetwork9, file = "cross_subnetwork9")
save(cross_subnetwork33, file = "cross_subnetwork33")
save(cross_subnetwork37, file = "cross_subnetwork37")
save(cross_subnetwork43, file = "cross_subnetwork43")
save(cross_subnetwork41, file = "cross_subnetwork41")
save(cross_subnetwork55, file = "cross_subnetwork55")
save(cross_subnetwork71, file = "cross_subnetwork71")
save(cross_subnetwork29, file = "cross_subnetwork29")
save(cross_subnetwork31, file = "cross_subnetwork31")
save(cross_subnetwork49, file = "cross_subnetwork49")
save(cross_subnetwork84, file = "cross_subnetwork84")
save(cross_subnetwork96, file = "cross_subnetwork96")


lay <- layout.auto(graph = cross_subnetwork5)



cross_subnetwork5 %>% 
  plot(., 
     vertex.label = NA,
     vertex.color = "#FFA319FF",
     vertex.size = vertex_attr(., "SS") * 10,
     # vertex.size = igraph::degree(graph = .)/5,
     # edge.color = "#767676FF",
     edge.color = ifelse(edge_attr(., "Correlation") > 0,
                         "#80000099", "#155F8399"),
     edge.width = abs(edge_attr(., "Correlation")) * 5,
     vertex.label.dist = 0.5,
     layout = lay)




#####pathway enrichment for each cross_subnetwork
cross_subnetwork5_node_info <-
  vertex_attr(cross_subnetwork5) %>% 
  bind_cols()


write.csv(cross_subnetwork5_node_info,
          file = "cross_subnetwork5_node_info.csv",
          row.names = FALSE)

###cross_subnetwork5
cross_subnetwork5_ID_trans <-
  read.table(
    "cross_subnetwork5_ID_trans.txt",
    sep = ",",
    header = TRUE,
    stringsAsFactors = FALSE
  )


cross_subnetwork5_cross_temp <- 
cross_subnetwork5_node_info %>% 
  select(., name, Compound.name, HMDB.ID, KEGG.ID) %>% 
  left_join(., cross_subnetwork5_ID_trans, by = c("Compound.name" = "Query"))


# cross_subnetwork5_cross_temp <- 
cross_subnetwork5_cross_temp$KEGG[!is.na(cross_subnetwork5_cross_temp$KEGG) & cross_subnetwork5_cross_temp$KEGG == ''] <- NA

cross_subnetwork5_cross_temp$HMDB[!is.na(cross_subnetwork5_cross_temp$HMDB) & cross_subnetwork5_cross_temp$HMDB == ''] <- NA


cross_subnetwork5_cross_temp$KEGG

HMDB.KEGG <-
  apply(cross_subnetwork5_cross_temp, 1, function(x){
  HMDB1 <- x[3]
  KEGG1 <- x[4]
  HMDB2 <- x[6]
  KEGG2 <- x[8]
  
  if(is.na(HMDB1)){
    HMDB <- HMDB2
  }else{
    HMDB <- HMDB1
  }
  
  if(is.na(KEGG1)){
    KEGG <- KEGG2
  }else{
    KEGG <- KEGG1
  }
  
  c(HMDB, KEGG)
  
}) %>% 
  t()


colnames(HMDB.KEGG) <- 
  c("HMDB.ID", "KEGG.ID")

cbind(cross_subnetwork5_node_info$HMDB.ID, HMDB.KEGG[,1])

cross_subnetwork5_node_info[,c("HMDB.ID", 'KEGG.ID')] <-
  HMDB.KEGG



write.csv(cross_subnetwork5_node_info, 
          "cross_subnetwork5_node_info.csv", row.names = FALSE)

id <- cross_subnetwork5_node_info$KEGG.ID
id <- id[!is.na(id)]


cross_pathway5 <- 
  enrichPathway(id = id, database = hsa.kegg.pathway)

save(cross_pathway5, file = "cross_pathway5")


cross_pathway5 %>% 
  mutate(Marker = ifelse(p.value < 0.05, "Yes", "No")) %>%
ggplot(., aes(Overlap.frac * 100, -log(p.value, 10))) +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2, colour = "#FFA319FF") +
  geom_point(aes(size = Pathway.length, colour = Marker), shape = 16) +
  scale_colour_manual(values = c("Yes" = "#FFA319FF", "No" = "#767676FF")) +
  guides(colour = FALSE, 
         size = guide_legend(title = "Pathway size",
                             title.theme = element_text(size = 15),
                             label.theme = element_text(size = 13))) +
  labs(x = "Overlap (%)", y = "-log10(P value)") +
  theme_bw() +
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 13),
        legend.position = c(1, 0), 
        legend.justification = c(1,0),
        legend.background = element_blank()) +
  ggrepel::geom_label_repel(aes(label = Pathway.name), 
                            data = cross_pathway5 %>% filter(., p.value < 0.05))


cross_pathway5 %>% 
  mutate(log.p = -log(p.value, 10)) %>% 
  mutate(enrichment = ifelse(p.value < 0.05, "Yes", "No")) %>% 
  `[`(., 1:15, ) %>% 
  # top_n(., 11, wt = "p.value") %>% 
  arrange(., log.p) %>% 
  ggplot(., aes(x = factor(Pathway.name, levels = Pathway.name), 
                y = log.p)) +  
geom_bar(stat = "identity", width = 0.8, 
         aes(fill = enrichment), show.legend = FALSE,
         # fill = "#FFA319FF"
         ) +
  scale_fill_manual(values = c("Yes" = "#FFA319FF", "No" = "#767676FF")) +
  labs(x = "", y = "-log10(P value)") +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))

  
###calculate the betweness and degree and all other
cross_betweenness5 <- betweenness(graph = cross_subnetwork5)
cross_degree5 <- igraph::degree(graph = cross_subnetwork5)
cross_closeness5 <- closeness(graph = cross_subnetwork5)

cross_importance5 <- 
tibble(cross_betweenness5, cross_degree5, cross_closeness5) %>% 
  mutate_all(., function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  mutate(.,mean = apply(., 1 ,mean)) %>% 
  mutate(peak.name = vertex_attr(graph = cross_subnetwork5, name = "name")) %>% 
  left_join(., cross_subnetwork5_node_info[,c("name", "Compound.name", "KEGG.ID")], 
            by = c('peak.name' = "name")) %>%  
  select(peak.name, Compound.name, everything()) %>% 
  arrange(., desc(mean))


# grep("Cysteine and methionine metabolism", names(hsa.kegg.pathway)) %>% 
#   `[[`(hsa.kegg.pathway, .) %>% 
#   intersect(., cross_importance5[c(1:50),]$KEGG.ID)
# 
# c('C00021', 'C00022')


cross_importance5[c(1:30),] %>% 
  arrange(., mean) %>% 
  mutate(Compound.name = factor(Compound.name, levels = Compound.name)) %>% 
  # mutate(node = ifelse(is.element(KEGG.ID, c('C00021', 'C00022')), 'Yes', "No")) %>% 
  ggplot(.,aes(x = Compound.name, y = mean)) +
  # geom_point(size = 2, aes(colour = node)) +
  geom_point(size = 2, colour = "#155F83FF") +
  # scale_y_continuous(expand = c(0,0.1)) +
  labs(x = "", y = "Mean (Betweenness, Degree, Closeness)") +
  geom_segment(aes(x = Compound.name, y = 0, 
                   xend = Compound.name, yend = mean),
               colour = "#155F83FF") +
  scale_colour_manual(values = c("Yes" = "#FFA319FF", "No" = "#767676FF")) +
  theme_bw() +
  coord_flip() +
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 10), legend.position = "none", 
        axis.text.x = element_text())


###betweenness, degree and closness
# test <- 
cross_importance5[c(1:30),] %>% 
  select(., Compound.name:mean) %>% 
  arrange(., mean) %>% 
  gather(., key = "Class", value = "Value", -Compound.name) %>% 
  mutate(Compound.name = factor(Compound.name, levels = Compound.name[1:30])) %>% 
  filter(., Class != "mean") %>% 
  ggplot(.,aes(x = Compound.name, 
               y = Value)) +
  geom_point(size = 2, 
             aes(colour = Class)
             # colour = "#155F83FF"
             ) +
  geom_segment(aes(x = Compound.name, y = 0, 
                   xend = Compound.name, 
                   yend = Value,
                   colour = Class)
               # colour = "#155F83FF"
               ) +
  scale_color_manual(values = c(
    "cross_betweenness5" = "#FFA319FF",
    "cross_closeness5" = "#8A9045FF",
    "cross_degree5" = "#C16622FF"
  )) +
  facet_grid(. ~ Class, 
             labeller = as_labeller(
               c("cross_betweenness5" = "Betweenness",
                 "cross_closeness5" = "Closeness",
                 "cross_degree5" = "Degree"
                 )
             )
             ) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 10), legend.position = "none", 
        axis.text.x = element_text(), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 13)) +
  coord_flip()

###find the nodes with top 3 imporance

cross_importance5[1,]$peak.name

cross_temp_subgraph <- 
neighbors(graph = cross_subnetwork5, 
          v = cross_importance5[1,]$peak.name) %>% 
  names() %>% 
  `c`(., cross_importance5[1,]$peak.name) %>% 
  unique() %>% 
  subgraph(cross_subnetwork5, v = .)


cross_temp_subgraph <- 
set_vertex_attr(graph = cross_temp_subgraph, name = "label", 
                value = ifelse(vertex_attr(cross_temp_subgraph, "name") == cross_importance5[1,]$peak.name, "Yes", "No"))

##get top 20 correlations
E(cross_temp_subgraph)
cross_temp_subgraph <- 
as_long_data_frame(cross_temp_subgraph) %>% 
  select(., from_name, to_name, Correlation) %>%
  filter(., from_name == cross_importance5[1,]$peak.name | to_name == cross_importance5[1,]$peak.name) %>% 
  mutate(Correlation = abs(Correlation)) %>% 
  arrange(., desc(Correlation)) %>% 
  top_n(., 20)
  
  
cross_temp_subgraph <- 
  subgraph(graph = cross_subnetwork5, 
           v = unique(c(cross_temp_subgraph$from_name, cross_temp_subgraph$to_name)))


cross_temp_subgraph <- 
  set_vertex_attr(graph = cross_temp_subgraph, name = "label", 
                  value = ifelse(vertex_attr(cross_temp_subgraph, "name") == cross_importance5[1,]$peak.name, "Yes", "No"))

which(vertex_attr(cross_temp_subgraph, "name") == cross_importance5[1,]$peak.name)

lay <- layout_as_tree(cross_temp_subgraph, root = 13)

  plot(cross_temp_subgraph,
       vertex.color = ifelse(vertex_attr(cross_temp_subgraph, "label")=="Yes", "#FFA31999", "#8A904599"),
       vertex.label = vertex_attr(cross_temp_subgraph, "Compound.name"),
       vertex.label.color = "black",
       vertex.size = vertex_attr(cross_temp_subgraph, "SS") * 20,
       vertex.label.dist = 1,
       
       # edge.color = "#767676FF",
       edge.color = ifelse(edge_attr(cross_temp_subgraph, "Correlation") > 0,
                           "#80000099", "#155F8399"),
       edge.width = abs(edge_attr(cross_temp_subgraph, "Correlation")) * 5,
       layout = -lay[, 2:1]
       )
  
 
  
  
  
  
  
  
  
  
  
  






  








