###FIRST TIME ONLY:
{
  install.packages("readxl")
  install.packages("mixtools")
  install.packages("ggpmisc")
  install.packages("ggplot2")
  install.packages("splus2R")
  install.packages("slider")
  install.packages("openxlsx")
  install.packages("data.table")
  install.packages("zoo")
  install.packages("pracma")
  install.packages("flextable")
  install.packages("officer")
  install.packages("webshot")
  install.packages("dplyr")
  install.packages("magick")
  install.packages("tidyverse")
  install.packages("styler")
  install_phantomjs()
}
###FIRST TIME ONLY ^^^


{
  library(readxl)
  library(mixtools)
  library(ggpmisc)
  library(ggplot2)
  library(splus2R)
  library(slider)
  library(openxlsx)
  library(data.table)
  library(zoo)
  library(pracma)
  library(flextable)
  library(officer)
  library(webshot)
  library(plyr)
  library(dplyr)
  library(magick)
  library(reprex)
  library(qpcR)
}



###Span good range 0.3 to 0.8, maxit good range 90 to 150, epsilon good range 1e-03 to 1e-08.
### Initial should be set to span = 0.5, Maxit=150, Epsiol=1e-06. If the code goes infinate set span = 0.7, and epsilon=1e-03.
span=0.5
maxit=150
epsilon = 1e-06

setwd("C:/Users/Owner/OneDrive/UMB/Deredge Lab/Deconvolution/DXPS_WT_and_Mutant/R420A")
data <- read_excel("R420A.xlsx")
num_peptime <- c (1, 2, 3 ,4 ,5 ,6)
charge_list <- c(2, 2, 2, 2, 4, 4)
data_list <- qpcR:::cbind.na(X41.56,X41.56.map,X183.196,X183.196.map,X278.298, X278.298.map)

{layout(mat = matrix(c(1, 2, 3, 4, 4, 4, 5, 6, 7, 8, 8, 8),
                      ncol=3, 
                      byrow=TRUE),
         heights = c(10,1, 10, 3),
         widths = c(10,10,10)
  )
  par(mar=c(1,1,1,1))
  layout.show(n=8)
  options(scipen=999)
  Peptide_Names <- c(colnames(data))
  output <- list()
}


for(k in num_peptime){
if(k>1 && k!=3){
{
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge_list[k]
  x <- data[print(Peptide_Names[k+(k-1)])]
  x <- na.omit(x[[1]])
  y.sg4 <- data[print(Peptide_Names[k+1+(k-1)])]
  y.sg4 <- na.omit(y.sg4[[1]])
  y.sg4 <- savgol(y.sg4, 3, 4, 0)
  df_smooth <- data.frame(x, y.sg4)
}
{
  mixture<-normalmixEM(x_i, sd.constr = c("a", "a"), arbmean = TRUE, arbvar = TRUE, maxit=maxit, epsilon = epsilon, ECM=TRUE)
  
  ptest <- predict(loess(y_i ~ x_i, family="gaussian", span= span), df)
  
  p_pred <- data.frame(x_i, ptest)
  
  amp1 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[1])
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1]){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nls(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
               data = df,
               start = list(mean1 = mixture$mu[1],
                            sigma1 = mixture$sigma[1],
                            C1 = amp1),
               lower =c(0,mixture$sigma[1],0),
               upper=c(Inf, mixture$sigma[1], 100),
               algorithm = "port",
    )
    
    dffit <- data.frame(x=x_i)
    dffit$y <- predict(fit, newdata=dffit)
    fit.sum <- summary(fit)
    fit.sum
    
    coef.fit <- fit.sum$coefficients
    mu.fit <- coef.fit[1]
    sigma.fit <- coef.fit[2]
    c.fit <- coef.fit[3]
  }else {
    x <- x_i
    amp1 <- approx(p_pred$x, p_pred$ptest, xout=mixture$mu[1])
    amp2 <- approx(p_pred$x, p_pred$ptest, xout=mixture$mu[2])

    
    fit <- nls(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                        C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
               data = df,
               start = list(mean1 = mixture$mu[1],
                            mean2 = mixture$mu[2],
                            sigma1 = mixture$sigma[1],
                            sigma2 = mixture$sigma[2],
                            C1 = amp1$y,
                            C2 = amp2$y),
               lower =c(0,0,mixture$sigma[1],mixture$sigma[2],0,0),
               upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], 100, 100),
               algorithm = "port",
    )
    
    dffit <- data.frame(x=x_i)
    dffit$y <- predict(fit, newdata=dffit)
    fit.sum <- summary(fit)
    fit.sum
    
    coef.fit <- fit.sum$coefficients[,1]
    mu.fit <- coef.fit[1:2]
    sigma.fit <- coef.fit[3:4]
    c.fit <- coef.fit[5:6]    
  }
}
{
  plot_x <- x_i
  plot_x <- append(x_i, x_i[1]-charge_shift, 0)
  plot_y <- append(y_i, 0, 0)
  plot_df <- data.frame(plot_x, plot_y)
  
  #original
  plot(df_smooth$y.sg4 ~ df_smooth$x, data = df_smooth, type = "l", main = print(Peptide_Names[k+(k-1)]), xlab="m/z", lwd=1.5, ylab="Intensity", font=2) 
  box(lwd=2)
  #overall fit
  lo_sum <- spline(plot_x,plot_y, method="natural")
  lines(lo_sum, col ="red", cex = 0.2, lwd=2.5) 
  #components of the fit
  for(i in 1:2){
    x <- dffit$x
    y <- (c.fit[i] *exp(-(x-mu.fit[i])**2/(2 * sigma.fit[i]**2)))
    x <- append(x, dffit$x[1]-charge_shift, 0)
    y <- append(y, 0, 0)
    lo <- spline(x, y, method="natural")
    lines(lo, col = i + 2, lwd=2.5)
  }
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  setequal(FWHM_1,FWHM_2)
  FWHM <- FWHM_1
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- plot_x
  AUC_tot_y <- plot_y
  id_2 <- order(plot_x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  centroid_idx <- which.max((x_i*y_i)/length(x_i))
  centroid <- x_i[centroid_idx]
  
  if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  } else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
  }

}
}
if(k==1){
    {
      data_x <- na.omit(data_list[k])
      data_y <- na.omit(data_list[k+1])
      x_i <- c(data_x$x)
      y_i <- c(data_y$y)
      df <- data.frame(x_i,y_i)
      charge_shift =1/charge_list[k]
      x <- data[print(Peptide_Names[k])]
      x <- na.omit(x[[1]])
      y.sg4 <- data[print(Peptide_Names[k+1])]
      y.sg4 <- na.omit(y.sg4[[1]])
      y.sg4 <- savgol(y.sg4, 3, 4, 0)
      df_smooth <- data.frame(x, y.sg4)
    }
    {
      mixture<-normalmixEM(x_i, sd.constr = c("a", "a"), arbmean = TRUE, arbvar = TRUE, maxit=maxit, epsilon = epsilon, ECM=TRUE)
      
      ptest <- predict(loess(y_i ~ x_i, family="gaussian", span= span), df)
      
      p_pred <- data.frame(x_i, ptest)
      
      amp1 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[1])
      
      amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
    }
    {
    if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1]){
      x <- x_i
      amp1 <- max(c(amp1$y, amp2$y))
      
      fit <- nls(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              sigma1 = mixture$sigma[1],
                              C1 = amp1),
                 lower =c(0,mixture$sigma[1],0),
                 upper=c(Inf, mixture$sigma[1], 100),
                 algorithm = "port",
      )
      
      dffit <- data.frame(x=x_i)
      dffit$y <- predict(fit, newdata=dffit)
      fit.sum <- summary(fit)
      fit.sum
      
      coef.fit <- fit.sum$coefficients
      mu.fit <- coef.fit[1]
      sigma.fit <- coef.fit[2]
      c.fit <- coef.fit[3]
    }else {
      x <- x_i
      amp1 <- approx(p_pred$x, p_pred$ptest, xout=mixture$mu[1])
      amp2 <- approx(p_pred$x, p_pred$ptest, xout=mixture$mu[2])
      
      
      fit <- nls(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower =c(0,0,mixture$sigma[1],mixture$sigma[2],0,0),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], 100, 100),
                 algorithm = "port",
      )
      
      dffit <- data.frame(x=x_i)
      dffit$y <- predict(fit, newdata=dffit)
      fit.sum <- summary(fit)
      fit.sum
      
      coef.fit <- fit.sum$coefficients[,1]
      mu.fit <- coef.fit[1:2]
      sigma.fit <- coef.fit[3:4]
      c.fit <- coef.fit[5:6]    
    }
  }
    {
      plot_x <- append(x_i, x_i[1]-charge_shift, 0)
      plot_y <- append(y_i, 0, 0)
      plot_df <- data.frame(plot_x, plot_y)
      
      #original
      plot(df_smooth$y.sg4 ~ df_smooth$x, data = df_smooth, type = "l", main = print(Peptide_Names[k+(k-1)]), xlab="m/z", lwd=1.5, ylab="Intensity", font=2) 
      box(lwd=2)
      #overall fit
      lo_sum <- spline(plot_x,plot_y, method="natural")
      lines(lo_sum, col ="red", cex = 0.2, lwd=2.5) 
      #components of the fit
      for(i in 1:2){
        x <- dffit$x
        y <- (c.fit[i] *exp(-(x-mu.fit[i])**2/(2 * sigma.fit[i]**2)))
        x <- append(x, dffit$x[1]-charge_shift, 0)
        y <- append(y, 0, 0)
        lo <- spline(x, y, method="natural")
        lines(lo, col = i + 2, lwd=2.5)
      }
      FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
      FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
      setequal(FWHM_1,FWHM_2)
      FWHM <- FWHM_1
      
      y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
      id <- order(dffit$x)
      AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
      AUC_1
      y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
      AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
      AUC_2
      AUC_tot_x <- plot_x
      AUC_tot_y <- plot_y
      id_2 <- order(plot_x)
      AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
      AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
      
      AUC_list <- c(AUC_1, AUC_2)
      AUC_tot <- c(AUC_tot_dat, NA)
      centroid_idx <- which.max((x_i*y_i)/length(x_i))
      centroid <- x_i[centroid_idx]
      
      if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k]," ", "Unimodal Peak")))
      } else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
        rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k]," ", "Bimodal Peak 2"))))
      }
      
    }
}
if(k==3){
    {
      data_x <- na.omit(data_list[k+(k-1)])
      data_y <- na.omit(data_list[k+1+(k-1)])
      x_i <- c(data_x$x)
      y_i <- c(data_y$y)
      df <- data.frame(x_i,y_i)
      charge_shift =1/charge_list[k]
      x <- data[print(Peptide_Names[k+(k-1)])]
      x <- na.omit(x[[1]])
      y.sg4 <- data[print(Peptide_Names[k+1+(k-1)])]
      y.sg4 <- na.omit(y.sg4[[1]])
      y.sg4 <- savgol(y.sg4, 3, 4, 0)
      df_smooth <- data.frame(x, y.sg4)
    }
    {
      mixture<-normalmixEM(x_i, sd.constr = c("a", "a"), arbmean = TRUE, arbvar = TRUE, maxit=maxit, epsilon = epsilon, ECM=TRUE)
      
      ptest <- predict(loess(y_i ~ x_i, family="gaussian", span= span), df)
      
      p_pred <- data.frame(x_i, ptest)
      
      amp1 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[1])
      
      amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
    }
    {
    if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1]){
      x <- x_i
      amp1 <- max(c(amp1$y, amp2$y))
      
      fit <- nls(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              sigma1 = mixture$sigma[1],
                              C1 = amp1),
                 lower =c(0,mixture$sigma[1],0),
                 upper=c(Inf, mixture$sigma[1], 100),
                 algorithm = "port",
      )
      
      dffit <- data.frame(x=x_i)
      dffit$y <- predict(fit, newdata=dffit)
      fit.sum <- summary(fit)
      fit.sum
      
      coef.fit <- fit.sum$coefficients
      mu.fit <- coef.fit[1]
      sigma.fit <- coef.fit[2]
      c.fit <- coef.fit[3]
    }else {
      x <- x_i
      amp1 <- approx(p_pred$x, p_pred$ptest, xout=mixture$mu[1])
      amp2 <- approx(p_pred$x, p_pred$ptest, xout=mixture$mu[2])
      
      
      fit <- nls(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower =c(0,0,mixture$sigma[1],mixture$sigma[2],0,0),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], 100, 100),
                 algorithm = "port",
      )
      
      dffit <- data.frame(x=x_i)
      dffit$y <- predict(fit, newdata=dffit)
      fit.sum <- summary(fit)
      fit.sum
      
      coef.fit <- fit.sum$coefficients[,1]
      mu.fit <- coef.fit[1:2]
      sigma.fit <- coef.fit[3:4]
      c.fit <- coef.fit[5:6]    
    }
  }
    {
      plot_x <- x_i
      plot_x <- append(x_i, x_i[1]-charge_shift, 0)
      plot_y <- append(y_i, 0, 0)
      plot_df <- data.frame(plot_x, plot_y)
      
      #original
      plot(df_smooth$y.sg4 ~ df_smooth$x, data = df_smooth, type = "l", main = print(Peptide_Names[k+(k-1)]), xlab="m/z", lwd=1.5, ylab="Intensity", font=2) 
      box(lwd=2)
      #overall fit
      lo_sum <- spline(plot_x,plot_y, method="natural")
      lines(lo_sum, col ="red", cex = 0.2, lwd=2.5) 
      #components of the fit
      for(i in 1:2){
        x <- dffit$x
        y <- (c.fit[i] *exp(-(x-mu.fit[i])**2/(2 * sigma.fit[i]**2)))
        x <- append(x, dffit$x[1]-charge_shift, 0)
        y <- append(y, 0, 0)
        lo <- spline(x, y, method="natural")
        lines(lo, col = i + 2, lwd=2.5)
      }
      FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
      FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
      setequal(FWHM_1,FWHM_2)
      FWHM <- FWHM_1
      
      y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
      id <- order(dffit$x)
      AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
      AUC_1
      y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
      AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
      AUC_2
      AUC_tot_x <- plot_x
      AUC_tot_y <- plot_y
      id_2 <- order(plot_x)
      AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
      AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
      
      AUC_list <- c(AUC_1, AUC_2)
      AUC_tot <- c(AUC_tot_dat, NA)
      centroid_idx <- which.max((x_i*y_i)/length(x_i))
      centroid <- x_i[centroid_idx]
      
      if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      } else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
        rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
      
      plot(1, type = "n", axes=F, xlab="", ylab="") # Create empty plot
        }
      
    }
  } 
}


###Finalize and Save Plot
{
  plot(1, type = "n", axes=F, xlab="", ylab="") # Create empty plot
  legend("bottom",lty = c(1,1,1), horiz = TRUE, col = c("red", "green", "blue", "black"), c("Sum", "Envelope 1", "Envelope 2", "Original"), lwd = 2)
  
  dev.copy(jpeg,'Final_Plot.jpeg', width=3840, height=2020, res=200)
  dev.size(units=c("px"))
  dev.off()
}

###Data Table and Export
{
output_total <- rbind(output[[1]],output[[2]],output[[3]],output[[4]],output[[5]],output[[6]])
  
output_total

write.csv(output_total, "Final_Data.csv")

Final_Data_frame<- output_total %>%
  as.data.frame() %>% 
  add_rownames() %>% 
  flextable()

fin_bord <- fp_border(color ="dimgray", style="solid", width=2)

Final_Data_frame<- Final_Data_frame %>%
  set_table_properties(layout = "autofit")%>%
  bg(bg="white", part="all")%>%
  hline() %>%
  set_header_labels(Final_Data_frame, rowname = "Timepoint")%>%
  bold(part="header")%>%
  bold(j=1)%>%
  vline(j=1, border = fin_bord )

autofit(Final_Data_frame)
save_as_image(Final_Data_frame, path ="Final_Data_Sample.jpeg")
}

