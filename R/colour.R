library(ggsci)

library(scales)

test_colour <- ggsci::pal_uchicago(alpha = 0.8)(10)
test_colour <- ggsci::pal_uchicago(alpha = 1)(10)
library(scales)
show_col(test_colour)

test_colour


my_theme <- 
    ggplot2::theme(axis.title = element_text(size = 1.5),
                   axis.text = element_text(size = 1.3))  
  
  
