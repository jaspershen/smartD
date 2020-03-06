sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/metabolites/RF/GA_prediction/")
marker1 <- readr::read_csv("marker_rf_final.csv")

sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/metabolites/RF/time_to_due_prediction/remove_cs/")
marker2 <- readr::read_csv("marker_rf_final.csv")


marker <- rbind(marker1, marker2) %>% 
  distinct(name, .keep_all = TRUE)

sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/metabolites/RF/all_marker_ms2_shape/peak_shape/")
xlsx::write.xlsx2(x = marker %>% select(name, mz, rt) %>% as.data.frame(), 
            file = "marker.xlsx", row.names = FALSE)


peak_data_pos <- 
  metflow2::extractPeaks(path = ".", ppm = 15, 
                         threads = 4, rt.tolerance = 20,
                         is.table = "marker.xlsx")



plot <- vector(mode = "list", length = nrow(marker))
for(i in 1:nrow(marker)){
  cat(i, " ")
  plot[[i]] <- metflow2::showPeak(object = peak_data_pos, 
                                  peak.index = i, 
                                  alpha = 0.5, 
                                  interactive = FALSE)
  
  plot[[i]] <- 
    plot[[i]] +
    labs(title = marker$name[i]) +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 6),
          axis.ticks = element_blank()) 
}


library(cowplot)

plot_grid(plot[[1]], plot[[2]],
          plot[[3]], plot[[4]],
          plot[[5]], plot[[6]],
          plot[[7]], plot[[8]],
          plot[[9]], plot[[10]],
          plot[[11]], plot[[12]],
          plot[[13]], plot[[14]],
          plot[[15]], plot[[16]],
          plot[[17]], plot[[18]],
          plot[[19]], plot[[20]],
          plot[[21]], plot[[22]],
          plot[[23]], plot[[24]],
          labels = letters[1:24], 
          ncol = 6)



###out put MS2 plot
sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/metabolites/RF/all_marker_ms2_shape/MS2_plot/")
load("result.pRPLC.nce25")
load("result.pRPLC.nce50")
load("HMDB.metabolite.data")
load("hmdbDatabase0.0.1")
load("massbankDatabase0.0.1")
load("metlinDatabase0.0.1")
load("monaDatabase0.0.1")
load("msDatabase_hilic0.0.1")
load("msDatabase_rplc0.0.1")
load("nistDatabase0.0.1")
load("orbitrapDatabase0.0.1")

i <- 24
temp_name <- 
marker$name[i]
temp_compund <- marker$Compound.name[i]
temp_compund
temp_database <- marker$Database[i]
temp_idx <- which(names(result.pRPLC.nce25) == temp_database)

plot <- 
metID::ms2plot(object = result.pRPLC.nce50[[temp_idx]], 
               database = get(temp_database), 
               which.peak = temp_name)

plot <- plot +
  labs(title = paste(temp_name, temp_database, temp_compund, sep = "/")) +
  theme(plot.title = element_text(size = 15))

plot
ggsave(filename = paste(temp_name, "pdf", sep = "."), width = 8, height = 6)




