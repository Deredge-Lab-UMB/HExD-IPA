###FIRST TIME ONLY:
{
  #install.packages("readxl")
  #install.packages("mixtools")
  #install.packages("ggpmisc")
  #install.packages("ggplot2")
  #install.packages("splus2R")
  #install.packages("slider")
  #install.packages("openxlsx")
  #install.packages("data.table")
  #install.packages("zoo")
  #install.packages("pracma")
  #install.packages("flextable")
  #install.packages("officer")
  #install.packages("webshot")
  #install.packages("dplyr")
  #install.packages("magick")
  #install.packages("tidyverse")
  #install.packages("styler")
  #install_phantomjs()
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
  library(ggformula)
  library(ggthemes)
  library(jpeg)
  library(grid)
  library(gridExtra)
  library(readr)
  show_col_types = FALSE
}




###Span good range 0.3 to 0.8, maxit good range 90 to 150, epsilon good range 1e-03 to 1e-08.
### Initial should be set to span = 0.5, Maxit=150, Epsiol=1e-06. If the code goes infinate set span = 0.7, and epsilon=1e-03.
span=0.5
maxit=150
epsilon = 1e-06

setwd("E:/For Vincent/p38/EX1/Alpha/130-145")
data <- read_csv("alpha 130-145.csv")
show_col_types = FALSE
num_peptime <- c (1, 2, 3 ,4 ,5 ,6,7,8)
time <- c(1, 10,60, (10*60),(10*60), (2*3600), (2*3600),(3*3600))
charge <- c(3)


#XUndeut <- undeut
#XTD <- TD


data_list <- qpcR:::cbind.na(XUndeut, X10.sec,X11.sec,X10.m,X11.m,X2.hr,X2.01h, XTD)

{
  Peptide_Names <- c(colnames(data))
  output <- list()
}

###Undeut
{
  k=1
  data_x <- na.omit(data_list[k])
  data_y <- na.omit(data_list[k+1])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p1 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
          )
  #overall fit
  p1<- p1+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p1 <- p1 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p1 <- p1 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p1)
  ggsave('p1.jpeg',p1,width=7,height=7)
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[1], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k]," ", "Unimodal Peak")))
}
}

###Second Time Point
{
  k=2
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
             lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
    plot_y <-y_i
    plot_df <- data.frame(plot_x, plot_y)
    
    smooth_x <- df_smooth$x
    smooth_y <- df_smooth$y.sg4
    
    #original
    p2 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
      geom_line(lwd=0.8)+
      ggtitle(print(Peptide_Names[k+(k-1)])) +
      theme_tufte()+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
            plot.title = element_text(color="black", size=15, face="bold"),
            axis.title = element_blank(),
            axis.text.x = element_text(face="bold", color="Black",size=15),
            axis.text.y = element_text(face="bold", color="Black",size=15)
      )
    #overall fit
    p2<- p2+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
    
    #components of the fit
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p2 <- p2 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p2 <- p2 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
    }
    
    plot(p2)
    ggsave('p2.jpeg',p2,width=7,height=7)
    
    
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
        rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
      }
  }else {
  output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
  colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
  rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

###Third Time Point
{
  k=3
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p3 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k+(k-1)])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
    )
  #overall fit
  p3<- p3+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p3 <- p3 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p3 <- p3 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p3)
  ggsave('p3.jpeg',p3,width=7,height=7)
  
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

###Fourth Time Point
{
  k=4
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p4 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k+(k-1)])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
    )
  #overall fit
  p4<- p4+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p4 <- p4 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p4 <- p4 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p4)
  ggsave('p4.jpeg',p4,width=7,height=7)

  
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

###Fifth Time Point
{
  k=5
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p5 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k+(k-1)])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
    )
  #overall fit
  p5<- p5+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p5 <- p5 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p5 <- p5 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p5)
  ggsave('p5.jpeg',p5,width=7,height=7)
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

###Sixth Time Point
{
  k=6
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p6 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k+(k-1)])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
    )
  #overall fit
  p6<- p6+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p6 <- p6 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p6 <- p6 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p6)
  ggsave('p6.jpeg',p6,width=7,height=7)
  
  
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

###Seventh Time Point
{
  k=7
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p7 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k+(k-1)])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
    )
  #overall fit
  p7<- p7+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p7 <- p7 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p7 <- p7 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p7)
  ggsave('p7.jpeg',p7,width=7,height=7)
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

###Eigth Time Point
{
  k=8
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p8 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k+(k-1)])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
    )
  #overall fit
  p8<- p8+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p8 <- p8 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p8 <- p8 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p8)
  ggsave('p8.jpeg',p8,width=7,height=7)
  
  
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

###Nineth Time Point
{
  k=9
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p9 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k+(k-1)])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
    )
  #overall fit
  p9<- p9+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p9 <- p9 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p9 <- p9 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p9)
  ggsave('p9.jpeg',p9,width=7,height=7)
  
  
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

###Tenth Time Point
{
  k=10
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p10 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k+(k-1)])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
    )
  #overall fit
  p10<- p10+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p10 <- p10 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p10 <- p10 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p10)
  ggsave('p10.jpeg',p10,width=7,height=7)
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

###Eleventh Time Point
{
  k=11
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p11 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k+(k-1)])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
    )
  #overall fit
  p11<- p11+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p11 <- p11 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p11 <- p11 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p11)
  ggsave('p11.jpeg',p11,width=7,height=7)
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

###Twelfth Time Point
{
  k=12
  data_x <- na.omit(data_list[k+(k-1)])
  data_y <- na.omit(data_list[k+1+(k-1)])
  x_i <- c(data_x$x)
  y_i <- c(data_y$y)
  df <- data.frame(x_i,y_i)
  charge_shift =1/charge
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
               lower =c(0,mixture$mu1,mixture$sigma[1],mixture$sigma[2],0,0),
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
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p12 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=0.8)+
    ggtitle(print(Peptide_Names[k+(k-1)])) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1),
          plot.title = element_text(color="black", size=15, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=15),
          axis.text.y = element_text(face="bold", color="Black",size=15)
    )
  #overall fit
  p12<- p12+geom_spline(data=plot_df, aes(plot_x,plot_y), color="red", lwd=1.3)
  
  #components of the fit
  if(length(mu.fit)>1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
    lo2 <- spline(x, y2, method="natural")
    p12 <- p12 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)+ 
      geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=4, lwd=1.3)
  }else{
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p12 <- p12 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=3, lwd=1.3)
  }
  
  plot(p12)
  ggsave('p12.jpeg',p12,width=7,height=7)
  
  
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
  
  if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
    if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, AUC_Norm, AUC_tot))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
      rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
    }
  }else {
    output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], c(1), AUC_tot[1]))
    colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Normalized Area', 'Area Total')
    rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
  }    
}

 ###Finalize and Save Plot
{
  if(length(num_peptime)==1){
    plot1 <- magick::image_read('p1.jpeg') 
    row1 <- image_append(c(plot1))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white")
    plotF<- magick::image_append(c(row1,legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==2){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    row1 <- image_append(c(plot1,plot2))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "1200")
    plotF<- magick::image_append(c(row1,legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==3){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    plot3 <- magick::image_read('p3.jpeg') 
    row1 <- image_append(c(plot1,plot2))
    row2 <- image_append(c(plot3))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "1200")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==4){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    plot3 <- magick::image_read('p3.jpeg') 
    plot4 <- magick::image_read('p4.jpeg') 
    row1 <- image_append(c(plot1,plot2))
    row2 <- image_append(c(plot3,plot4))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "1200")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==5){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    plot3 <- magick::image_read('p3.jpeg') 
    plot4 <- magick::image_read('p4.jpeg') 
    plot5 <- magick::image_read('p5.jpeg') 
    row1 <- image_append(c(plot1,plot2,plot3))
    row2 <- image_append(c(plot4,plot5))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "2500")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==6){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    plot3 <- magick::image_read('p3.jpeg') 
    plot4 <- magick::image_read('p4.jpeg') 
    plot5 <- magick::image_read('p5.jpeg') 
    plot6 <- magick::image_read('p6.jpeg') 
    row1 <- image_append(c(plot1,plot2,plot3))
    row2 <- image_append(c(plot4,plot5,plot6))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "2500")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==7){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    plot3 <- magick::image_read('p3.jpeg') 
    plot4 <- magick::image_read('p4.jpeg') 
    plot5 <- magick::image_read('p5.jpeg') 
    plot6 <- magick::image_read('p6.jpeg')
    plot7 <- magick::image_read('p7.jpeg') 
    row1 <- image_append(c(plot1,plot2,plot3,plot4))
    row2 <- image_append(c(plot5,plot6,plot7))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "3500")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==8){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    plot3 <- magick::image_read('p3.jpeg') 
    plot4 <- magick::image_read('p4.jpeg') 
    plot5 <- magick::image_read('p5.jpeg') 
    plot6 <- magick::image_read('p6.jpeg')
    plot7 <- magick::image_read('p7.jpeg') 
    plot8 <- magick::image_read('p8.jpeg') 
    row1 <- image_append(c(plot1,plot2,plot3, plot4))
    row2 <- image_append(c(plot5,plot6, plot7, plot8))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "3500")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==9){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    plot3 <- magick::image_read('p3.jpeg') 
    plot4 <- magick::image_read('p4.jpeg') 
    plot5 <- magick::image_read('p5.jpeg') 
    plot6 <- magick::image_read('p6.jpeg') 
    plot7 <- magick::image_read('p7.jpeg') 
    plot8 <- magick::image_read('p8.jpeg') 
    plot9 <- magick::image_read('p9.jpeg') 
    row1 <- image_append(c(plot1,plot2,plot3))
    row2 <- image_append(c(plot4,plot5,plot6))
    row3 <- image_append(c(plot7,plot8,plot9))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "2500")
    plotF<- magick::image_append(c(row1,row2,row3, legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==10){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    plot3 <- magick::image_read('p3.jpeg') 
    plot4 <- magick::image_read('p4.jpeg') 
    plot5 <- magick::image_read('p5.jpeg') 
    plot6 <- magick::image_read('p6.jpeg') 
    plot7 <- magick::image_read('p7.jpeg') 
    plot8 <- magick::image_read('p8.jpeg') 
    plot9 <- magick::image_read('p9.jpeg') 
    plot10 <- magick::image_read('p10.jpeg') 
    row1 <- image_append(c(plot1,plot2,plot3,plot4))
    row2 <- image_append(c(plot5,plot6,plot7,plot8))
    row3 <- image_append(c(plot9,plot10))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "3500")
    plotF<- magick::image_append(c(row1,row2,row3, legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==11){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    plot3 <- magick::image_read('p3.jpeg') 
    plot4 <- magick::image_read('p4.jpeg') 
    plot5 <- magick::image_read('p5.jpeg') 
    plot6 <- magick::image_read('p6.jpeg') 
    plot7 <- magick::image_read('p7.jpeg') 
    plot8 <- magick::image_read('p8.jpeg') 
    plot9 <- magick::image_read('p9.jpeg') 
    plot10 <- magick::image_read('p10.jpeg') 
    plot11 <- magick::image_read('p11.jpeg') 
    row1 <- image_append(c(plot1,plot2,plot3,plot4))
    row2 <- image_append(c(plot5,plot6,plot7,plot8))
    row3 <- image_append(c(plot9,plot10,plot11))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "3500")
    plotF<- magick::image_append(c(row1,row2,row3, legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==12){
    plot1 <- magick::image_read('p1.jpeg') 
    plot2 <- magick::image_read('p2.jpeg') 
    plot3 <- magick::image_read('p3.jpeg') 
    plot4 <- magick::image_read('p4.jpeg') 
    plot5 <- magick::image_read('p5.jpeg') 
    plot6 <- magick::image_read('p6.jpeg') 
    plot7 <- magick::image_read('p7.jpeg') 
    plot8 <- magick::image_read('p8.jpeg') 
    plot9 <- magick::image_read('p9.jpeg') 
    plot10 <- magick::image_read('p10.jpeg') 
    plot11 <- magick::image_read('p11.jpeg') 
    plot12 <- magick::image_read('p12.jpeg')
    row1 <- image_append(c(plot1,plot2,plot3,plot4))
    row2 <- image_append(c(plot5,plot6,plot7,plot8))
    row3 <- image_append(c(plot9,plot10,plot11,plot12))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "3500")
    plotF<- magick::image_append(c(row1,row2,row3, legend), stack=TRUE)
    image_write(plotF, path = "Final.jpeg", format = "jpeg", quality = 100)
  }
  
}

###Data Table and Export
{
  output_total = data.frame()
for(z in num_peptime){
  temp_df<- (output[[z]])
  output_total <- rbind(output_total,temp_df)
}
  
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

###Normalized Area Data
{
output_tot_env1_in <- output_total[grep(pattern = "Bimodal Peak 1|Unimodal Peak", x=rownames(output_total)),]
output_tot_env2_in <-  output_total[grep(pattern = "Bimodal Peak 2", x=rownames(output_total)),]
ex_tran <- which((output_tot_env1_in$`Normalized Area`)%in%min(output_tot_env1_in$`Normalized Area`))


Uni_Chk <- output_total[grep(pattern = "Unimodal Peak", x=rownames(output_total)),]

Uni_Idx <- which((output_tot_env1_in$`Mean Value`)%in%(Uni_Chk$`Mean Value`))

Uni_Chk2 <- output_tot_env2_in[grep(pattern = "Unimodal Peak", x=rownames(output_tot_env2_in)),]

Uni_Idx2 <- which((output_tot_env2_in$`Mean Value`)%in%(Uni_Chk2$`Mean Value`))

newrow <- c(NA, NA, NA, NA, NA, NA, 0, NA)

if(length(Uni_Idx2)>0){
  Uni_Chk <- output_total[grep(pattern = "Unimodal Peak", x=rownames(output_total[1:ex_tran])),]
  
  Uni_Idx <- which((output_tot_env1$`Mean Value`)%in%(Uni_Chk$`Mean Value`))
  
  output_tot_env1_sub <- output_tot_env1_in$`Mean Value`[1:ex_tran]
  output_tot_env1 <- subset(output_tot_env1_in,output_tot_env1_in == output_tot_env1_sub)
  output_tot_env2_add <- na.omit(subset(output_tot_env1_in,output_tot_env1_in != output_tot_env1_sub))
  output_tot_env2 <- rbind(output_tot_env2_in,output_tot_env2_add)
  
for(i in Uni_Idx){
  if(i==1){
    output_tot_env2 <- rbind(newrow,output_tot_env2)
  }
}

for(i in Uni_Idx){
  if(i > 1){
    r = i
    newrow_1 <- c(NA, NA, NA, NA, NA, NA, 0, NA)
    output_tot_env2 <- rbind(output_tot_env2[1:(i-1), ],newrow_1, output_tot_env2[- (1:(i-1)), ])
  } 
}

for(i in Uni_Idx2){
  if(i==ex_tran){
    output_tot_env1 <- rbind(output_tot_env1,newrow)
  }
}

for(i in Uni_Idx2){
  if(i > ex_tran){
    r = i
    newrow_1 <- c(NA, NA, NA, NA, NA, NA, 0, NA)
    output_tot_env1 <- rbind(output_tot_env1[1:(i-1), ],newrow_1, output_tot_env1[- (1:(i-1)), ])
  } 
}
}

if(length(Uni_Idx2)==0){
  
  output_tot_env1<- output_total[grep(pattern = "Bimodal Peak 1|Unimodal Peak", x=rownames(output_total)),]
  
  output_tot_env2 <-  output_total[grep(pattern = "Bimodal Peak 2", x=rownames(output_total)),]
  
  for(i in Uni_Idx){
    if(i==1){
      output_tot_env2 <- rbind(newrow,output_tot_env2)
    }
  }
  
  for(i in Uni_Idx){
    if(i > 1){
      r = i
      newrow_1 <- c(NA, NA, NA, NA, NA, NA, 0, NA)
      output_tot_env2 <- rbind(output_tot_env2[1:(i-1), ],newrow_1, output_tot_env2[- (1:(i-1)), ])
    } 
  }
}
  
A_env1 <- (output_tot_env1$`Normalized Area`)
A_env2 <- (output_tot_env2$`Normalized Area`)
Area_df <- data.frame(time,A_env1,A_env2)

colnames(Area_df) <- c("Time (Sec)","Normalized Area of Envelope 1","Normalized Area of Envelope 2")

write.csv(Area_df, "Normalized_Data.csv")

Final_Det_frame<- Area_df %>%
  as.data.frame() %>% 
  flextable()

fin_bord <- fp_border(color ="dimgray", style="solid", width=2)

Final_Det_frame<- Final_Det_frame %>%
  set_table_properties(layout = "autofit")%>%
  bg(bg="white", part="all")%>%
  hline() %>%
  set_header_labels(Final_Det_frame, rowname = "Timepoint")%>%
  bold(part="header")%>%
  bold(j=1)%>%
  vline(j=1, border = fin_bord )

autofit(Final_Det_frame)
save_as_image(Final_Det_frame, path ="Normalized_Area_Data.jpeg")

layout(mat = matrix(c(1, 2, 3, 3, 4,4),
                    ncol=2, 
                    byrow=TRUE),
       heights = c(10,1.7,0.00000001),
       widths = c(10,10)
)
par(mar=c(4,4,1,1))
layout.show(n=4)

y_rel_max<- max(A_env1,A_env2)

plot(log(time),A_env1, col="green", type="b", xlab="Log of Time (Sec)", ylab="Normalized Areas", main="Log Linear Graph of Normalized Areas vs Time", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_rel_max+0.1)))
lines(log(time),A_env2, col="blue", type="b", lwd=2)
plot(time,A_env1, type="n", xlab="Time (Sec)", ylab="", main="Lines of Best Fit of Normalized Areas vs. Time", font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_rel_max+0.1)))
abline(lm(A_env1~time), col="green", lty=1, lwd=2)
abline(lm(A_env2~time), col="blue", lty=1, lwd=2)

A_env1_bf <-(lm(A_env1~time))
A_env2_bf <-(lm(A_env2~time))


plot(1, type = "n", axes=F, xlab="", ylab="") # Create empty plot
legend("center", lty = c(1,1), horiz = TRUE, col = c("green", "blue"), c("Envelope 1", "Envelope 2"), lwd = 2, cex=1.1)

dev.copy(jpeg,'Normalized_Area.jpeg', width=3840, height=2020, res=200)
dev.size(units=c("px"))
dev.off()

graphics.off()
}
