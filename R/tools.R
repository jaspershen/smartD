plotLambdaVSdeviation <- function(object,
                                  xlab = "Log lambda",
                                  ylab = "Deviation ratio (%)") {
  data.frame(
    lambda = object$lambda,
    dev.ratio = object$dev.ratio,
    stringsAsFactors = FALSE
  ) %>%
    ggplot(aes(log(lambda), dev.ratio*100)) +
    labs(x = xlab, y = ylab) +
    geom_point(size = 2) +
    geom_line() +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 13))
  
}


plotLambdaVScoefficients <- 
  function(
    object,
    xlab = "Log lambda", 
    ylab = "Coefficients"
  ){
    
    beta <-
      object$beta %>%
      as.matrix() %>%
      t() %>%
      as_tibble() %>%
      mutate(lambda = lasso_regression$lambda) %>%
      tidyr::gather(., key = "feature", value = "coef", -lambda)
    
    label_index <- seq(range(log(beta$lambda))[1], 
                 range(log(beta$lambda))[2], by = 1)
    label <- 
    lapply(label_index, function(x){
      c(log(object$lambda) - x) %>% 
        abs() %>% 
        which.min() %>% 
        `[`(object$df + 1, .)
    }) %>% 
      unlist()
    
    
    beta %>% 
      ggplot(., aes(log(lambda), coef)) +
      geom_line(aes(colour = feature), show.legend = FALSE) +
      scale_x_continuous(position = "bottom",
                         sec.axis = sec_axis(~., name = "", 
                                             breaks = label_index,
                                             labels = label)) +
      scale_colour_manual(values = colorRampPalette(ggsci::pal_uchicago()(5))(600)) +
      labs(x = xlab, y = ylab) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13)
      ) 
    
  }


plotLambdaVSerror <- function(object,
                              xlab = "Log lambda",
                              ylab = "Mean absolute error") {
  cvm <-
    data.frame(
      lambda = object$lambda,
      df = object$glmnet.fit$df,
      cvm = object$cvm,
      cvup = object$cvup,
      cvlo = object$cvlo,
      stringsAsFactors = FALSE
    )
  
  cvm %>%
    ggplot(., aes(log(lambda), cvm)) +
    geom_vline(xintercept = log(
      c(
        lasso_regression2$lambda.min,
        lasso_regression2$lambda.1se
      )
    ),
    linetype = 2) +
    geom_errorbar(aes(ymin = cvlo, ymax = cvup), colour = "#155F83FF") +
    geom_point(size = 2, colour = "#FFA319FF") +
    scale_x_continuous(
      position = "bottom",
      sec.axis = sec_axis(
        trans = ~ .,
        breaks = log(cvm$lambda)[seq(1, 100, by = 7)],
        labels = cvm$df[seq(1, 100, by = 7)],
        name = ""
      )
    ) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 13)
          # plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"))
    )
          
          
}
