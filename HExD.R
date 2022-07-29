{
  library(mixtools)
  library(ggpmisc)
  library(ggplot2)
  library(splus2R)
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
  library(ggformula)
  library(ggthemes)
  library(jpeg)
  library(grid)
  library(gridExtra)
  library(minpack.lm)
  library(readr)
  show_col_types = FALSE
} #Call Packages Needed


###Span good range 0.3 to 0.8, maxit good range 90 to 150, epsilon good range 1e-03 to 1e-08.
### Initial should be set to span = 0.5, Maxit=150, Epsiol=1e-06. If the code goes infinate set span = 0.7, and epsilon=1e-03.
{
span=0.5
maxit=150
epsilon = 1e-06
}


#Data Setup
{
setwd("C:/Users/Owner/OneDrive/UMB/Deredge Lab/For Vincent/FliD/EX1/361-375")
data <- read_csv("361 to 375.csv")
show_col_types = FALSE
num_peptime <- c (1, 2, 3 ,4 ,5 ,6,7)
charge_list <- c(2,5)
Pep_Name <- c("FliD_361-375")
}


#Plot Setup
{
lwd_mz <- 0.8 #line width of mass spec data graph
lwd_env <-1.3 #line width of all envelopes
axis_font <- 15 #axis font size
border_width <- 1.1 #width of black border
}


#Select Analysis Mode
{
time <- c() #time series data, if none leave as c()
mutant <- c() #time series data with, if none leave as c()
temp <- c() #temperature change data, if none leave as c()
conc <- c(1,2,3,4,5,6) #Concentration change data, if none leave as c()
}


#Analysis Graph Modes
{
Log=TRUE #used to graph the X-axis as a log of the value = TRUE if log, = FALSE if linear
Linear=FALSE #in reference to a linear graphing of x-values (not taking the log of it). = TRUE if linear, = FALSE if log
Line = TRUE #Plot a line graph
Bar = FALSE #plot a bar graph
}


###Code Cleaning (Once after IPA, then comment out again):
#XUndeut <- undeut
#XTD <- TD


data_list <- qpcR:::cbind.na(XTD, XUndeut, X10.sec, X1.min, X10.min, X1.hr, X2.hr) #Load in IPA data

#Initial Values Call, No need to change ever
{
  Peptide_Names <- c(colnames(data))
  output <- list()
  c.fit_tot <- list()
  x_max <- max(na.omit(data))
  x_min <- min(na.omit(data$TD))
}

#Important notes:
#1a) A good number of iterations is between 15 to 60 Iterations for a simple system
#1b) If the system shows unimodality and a bad fit/if the system goes infinate within NLSLM fitting then aim for 3 iterations
#2) The first thing to change with any errors is the Span. A good test is to go up and down within the listed range
#3) If the span doesn't work, change the epsilon
#4) If nothing works, check the IPA to ensure the proper peaks were selected
#5) If the timepoint flags a false bimodal at a later timepoint, change Num_envelops to =1 and color_1 = 4
#6) if the colors of the two envelops are switched, change color_1 to =4 and color_3 to =3 and rerun just the plot section
#7) If there is a Chi squared error it means that the fit is too far from the original data and should be rerun
#8) When running the code, be sure to have the plots tab open to visually catch any fitting issues
#9) If a block in a timepoint section fails, address the error if needed and then rerun startng from the data assignment block of that section
#10) Each block under each timepiont section is the same, so the naming is only applied to the first 4 blocks for reference.

###First Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE


{
  k=1
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
} #Data Assignment and Initial Calls
{
  mixture<-normalmixEM(x_i, sd.constr = c("a", "a"), arbmean = TRUE, arbvar = TRUE, maxit=maxit, epsilon = epsilon, ECM=TRUE)
  
  ptest <- predict(loess(y_i ~ x_i, family="gaussian", span= span), df)
  
  p_pred <- data.frame(x_i, ptest)
  
  amp1 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[1])
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
} #Mixed Effects Modeling and Amplitude Prediction via regression (finds peak at non-discrete datapoint)
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }else{
      c1_high <- 100
      c2_low <- 0
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
} #NLSLM (Non-Linear Least Squares Levenberg-Marquardt) fitting. Fits the two bimodal envelops
{
  plot_x <- x_i
  plot_y <-y_i
  plot_df <- data.frame(plot_x, plot_y)
  
  smooth_x <- df_smooth$x
  smooth_y <- df_smooth$y.sg4
  
  #original
  p1 <- ggplot(data=df_smooth, aes(smooth_x,smooth_y), colour="black") +
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p1<- p1+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p1 <- p1 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
    ADI_1 <- abs(0)
    ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p1 <- p1 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p1 <- p1 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- 0
      ADI_2 <- 0
      ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p1 <- p1 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
      if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- 0
      ADI_df <- c(ADI_1, 0)
      }
    }
  }
  plot(p1)
  ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p1,width=7,height=7)
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
  if(Num_Envelops==1 && Color_1==3){
    if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
      if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }
  }
  
  if(Num_Envelops==1 && Color_1==4){
    if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
      if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
    }
  }
  
  if(Num_Envelops==2){
    if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
      if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
        rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
      }
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }
  }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
} #Plot creation, Data table creation, and Chi Square Testing

###Second Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=2
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p2<- p2+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p2 <- p2 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
    ADI_1 <- abs(centroid-output[[1]]$Centroid)
    ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p2 <- p2 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p2 <- p2 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
      ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p2 <- p2 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p2))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p2,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
  if(Num_Envelops==1 && Color_1==3){
    if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
      if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }
  }
  
  if(Num_Envelops==1 && Color_1==4){
    if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
      if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
    }
  }
  
  if(Num_Envelops==2){
    if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
      if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
        rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
      }
    }else {
      output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
      colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
      rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
    }
  }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}

###Third Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=3
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p3<- p3+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p3 <- p3 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p3 <- p3 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p3 <- p3 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
        ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p3 <- p3 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(centroid-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p3))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p3,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}

###Fourth Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=4
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p4<- p4+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p4 <- p4 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p4 <- p4 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p4 <- p4 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
        ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p4 <- p4 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(centroid-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p4))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p4,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}

###Fifth Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=5
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p5<- p5+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p5 <- p5 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p5 <- p5 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p5 <- p5 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
        ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p5 <- p5 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(centroid-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p5))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p5,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}

###Sixth Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=6
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p6<- p6+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p6 <- p6 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p6 <- p6 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p6 <- p6 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
        ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p6 <- p6 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(centroid-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p6))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p6,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}

###Seventh Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=7
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p7<- p7+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p7 <- p7 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p7 <- p7 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p7 <- p7 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
        ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p7 <- p7 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(centroid-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p7))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p7,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}

###Eigth Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=8
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p8<- p8+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p8 <- p8 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p8 <- p8 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p8 <- p8 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
        ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p8 <- p8 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(centroid-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p8))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p8,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}

###Nineth Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=9
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p9<- p9+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p9 <- p9 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p9 <- p9 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p9 <- p9 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
        ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p9 <- p9 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(centroid-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p9))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p9,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}

###Tenth Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=10
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p10<- p10+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p10 <- p10 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p10 <- p10 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p10 <- p10 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
        ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p10 <- p10 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(centroid-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p10))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p10,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}

###Eleventh Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=11
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p11<- p11+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p11 <- p11 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p11 <- p11 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p11 <- p11 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
        ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p11 <- p11 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(centroid-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p11))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p11,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}

###Twelfth Time Point
Color_1=3
Color_2=4
Num_Envelops=2
Replicate=FALSE

{
  k=12
  if (length(charge_list)>1){
    charge <- charge_list[k]
  }else{
    charge <- charge_list[1]
  }
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
  
  amp1_test <- amp1
  
  amp2 <- approx(p_pred$x_i, p_pred$ptest, xout=mixture$mu[2])
  
  amp2_test <-amp2
}
{
  if(mixture$mu[1] == mixture$mu[2] | (mixture$mu[2] - 1) < mixture$mu[1] | Num_Envelops == 1){
    x <- x_i
    amp1 <- max(c(amp1$y, amp2$y))
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
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
    if(length(time)!=0 && k!=1 | length(conc)!=0 && k!=1){
      if(Replicate=FALSE){
        c1_high <- output[[k-1]]$Height[1]
        c2_low <- output[[k-1]]$Height[2]
      }
      if(Replicate=TRUE){
        c1_high <- 100
        c2_low <- 0
      }
      
    }
    
    fit <- nlsLM(y_i ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                          C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
                 data = df,
                 start = list(mean1 = mixture$mu[1],
                              mean2 = mixture$mu[2],
                              sigma1 = mixture$sigma[1],
                              sigma2 = mixture$sigma[2],
                              C1 = amp1$y,
                              C2 = amp2$y),
                 lower=c(0,mixture$mu[1],mixture$sigma[1],mixture$sigma[2],0,c2_low),
                 upper=c(Inf, Inf, mixture$sigma[1], mixture$sigma[2], c1_high, 100),
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
    geom_line(lwd=lwd_mz)+
    xlim(x_min,x_max)+
    ggtitle(c(print(paste0(Peptide_Names[k+(k-1)]," ", paste0(Pep_Name))))) +
    theme_tufte()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
          plot.title = element_text(color="black", size=axis_font, face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold", color="Black",size=axis_font),
          axis.text.y = element_text(face="bold", color="Black",size=axis_font)
    )
  #overall fit
  p12<- p12+geom_line(data=plot_df, aes(plot_x,plot_y), color="red", lwd=lwd_env)
  
  #components of the fit
  if(Num_Envelops==1){
    x <- dffit$x
    y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
    lo <- spline(x, y, method="natural")
    p12 <- p12 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
    centroid <- ((sum(x_i*y_i))/sum(y_i))
    plot_sum_y <- y
    plot_sum_df <- data.frame(x,plot_sum_y)
    if(length(time)!=0 | length(conc)!=0){
      ADI_1 <- abs(centroid-output[[1]]$Centroid)
      ADI_df <- c(ADI_1, 0)
    }
    FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
  }else{
    if(length(mu.fit)>1){
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      y2 <- (c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2)))
      lo2 <- spline(x, y2, method="natural")
      p12 <- p12 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)+ 
        geom_line(data=data.frame(lo2), aes(lo2$x,lo2$y),color=Color_2, lwd=lwd_env)
      plot_sum_y <- y+y2
      plot_sum_df <- data.frame(x,plot_sum_y)
      p12 <- p12 + geom_line(data=plot_sum_df, aes(x,plot_sum_y), color="#E69F00", lwd=lwd_env)
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(mu.fit[1]-output[[1]]$Centroid)
        ADI_2 <- abs(mu.fit[2]-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, ADI_2)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*(sqrt((sigma.fit[1]^2)+(sigma.fit[2]^2)))
    }else{
      x <- dffit$x
      y <- (c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2)))
      lo <- spline(x, y, method="natural")
      p12 <- p12 + geom_line(data=data.frame(lo), aes(lo$x,lo$y),color=Color_1, lwd=lwd_env)
      plot_sum_y <- y
      plot_sum_df <- data.frame(x,plot_sum_y)
      centroid <- ((sum(x_i*y_i))/sum(y_i))
      if(length(time)!=0 | length(conc)!=0){
        ADI_1 <- abs(centroid-output[[1]]$Centroid)
        ADI_df <- c(ADI_1, 0)
      }
      FWHM_Sum = 2*sqrt(log(2)*2)*sigma.fit[1]
    }
  }
  suppressWarnings(plot(p12))
  suppressWarnings(ggsave(paste0('p',k,"_",Pep_Name,".jpeg"),p12,width=7,height=7))
  
  FWHM_1 = 2*sqrt(log(2)*2)*sigma.fit[1]
  FWHM_2 = 2*sqrt(log(2)*2)*sigma.fit[2]
  FWHM <- FWHM_1
  FWHM_Tot <- c(FWHM_Sum,NA)
  
  y_AUC_1 <- c.fit[1] *exp(-(x-mu.fit[1])**2/(2 * sigma.fit[1]**2))
  id <- order(dffit$x)
  AUC_1 <- sum(diff(dffit$x[id])*rollmean(y_AUC_1[id],2))
  AUC_1
  y_AUC_2 <- c.fit[2] *exp(-(x-mu.fit[2])**2/(2 * sigma.fit[2]**2))
  AUC_2 <- sum(diff(dffit$x[id])*rollmean(y_AUC_2[id],2))
  AUC_2
  AUC_tot_x <- dffit$x
  AUC_tot_y <- plot_sum_y
  id_2 <- order(dffit$x)
  AUC_tot_dat <- sum(diff(AUC_tot_x[id_2])*rollmean(AUC_tot_y[id_2],2))
  AUC_peaks_x <- plot_x
  AUC_peaks_y <- plot_y
  id_3 <- order(plot_x)
  AUC_peaks_dat <- sum(diff(AUC_peaks_x[id_3])*rollmean(AUC_peaks_y[id_3],2))
  AUC_Norm <- c((AUC_1/AUC_tot_dat),(AUC_2/AUC_tot_dat))
  
  AUC_list <- c(AUC_1, AUC_2)
  AUC_tot <- c(AUC_tot_dat, NA)
  
  Chi_sqr <- (sum(plot_sum_y - plot_y)^2)/(sum(plot_y))
  Chi_crit <- qchisq(p=0.95, df=(length(plot_y)-1))
  if(Chi_sqr > Chi_crit){
    stop("Chi Squared Value Past Critical Point")
  }
  
  centroid <- ((sum(x_i*y_i))/sum(y_i))
  if(amp1_test$y<amp2_test$y){
    c.fit_tot[[k]] <- c(0,c.fit)
  }else{
    c.fit_tot[[k]] <-c.fit
  }
  if(length(time)!=0 | length(conc)!=0){
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, ADI_df, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1),ADI_df[1], AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area','Relative DA', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }else{
    if(Num_Envelops==1 && Color_1==3){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area','Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area',, 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
    
    if(Num_Envelops==1 && Color_1==4){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1),AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak Env2")))
      }
    }
    
    if(Num_Envelops==2){
      if(is.na(AUC_1) == FALSE && is.na(AUC_2)==FALSE){
        if(AUC_Norm[1] >= 0.9 | AUC_Norm[2] >= 0.9 | AUC_Norm[1]<= 0.1 | AUC_Norm[2] <=0.1){
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1], Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area', 'Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
        }else {
          output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit, centroid, c.fit, sigma.fit, FWHM, AUC_list, Chi_sqr, AUC_Norm, AUC_tot, FWHM_Tot))
          colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total', 'Total Width')
          rownames(output[[k]]) <- c((print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 1"))), (print(paste0(Peptide_Names[k+(k-1)]," ", "Bimodal Peak 2"))))
        }
      }else {
        output[[k]] <- assign(Peptide_Names[k+(k-1)], data.frame(mu.fit[1], centroid[1], c.fit[1], sigma.fit[1], FWHM[1], AUC_tot[1],Chi_sqr, c(1), AUC_tot[1], FWHM_Tot[1]))
        colnames(output[[k]]) <- c('Mean Value','Centroid','Height','SD','FWHM', 'Area','Chi Squared', 'Normalized Area', 'Area Total','Total Width')
        rownames(output[[k]]) <- c(print(paste0(Peptide_Names[k+(k-1)]," ", "Unimodal Peak")))
      }
    }
  }
}



###Finalize and Save Plot, Creates the Grid Plot Graph
{
  if(length(num_peptime)==1){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot1 <- image_border(plot1, "white","1000")
    row1 <- image_append(c(plot1))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "300")
    plotF<- magick::image_append(c(row1,legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==2){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))  
    row1 <- image_append(c(plot1,plot2))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "300")
    plotF<- magick::image_append(c(row1,legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==3){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))
    plot3 <- magick::image_read(paste0('p3',"_",Pep_Name,'.jpeg')) 
    row1 <- image_append(c(plot1,plot2))
    row2 <- image_append(c(plot3))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "300")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==4){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))
    plot3 <- magick::image_read(paste0('p3',"_",Pep_Name,'.jpeg')) 
    plot4 <- magick::image_read(paste0('p4',"_",Pep_Name,'.jpeg'))
    row1 <- image_append(c(plot1,plot2))
    row2 <- image_append(c(plot3,plot4))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "300")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==5){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))
    plot3 <- magick::image_read(paste0('p3',"_",Pep_Name,'.jpeg')) 
    plot4 <- magick::image_read(paste0('p4',"_",Pep_Name,'.jpeg'))
    plot5 <- magick::image_read(paste0('p5',"_",Pep_Name,'.jpeg'))
    row1 <- image_append(c(plot1,plot2,plot3))
    row2 <- image_append(c(plot4,plot5))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "1500")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==6){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))
    plot3 <- magick::image_read(paste0('p3',"_",Pep_Name,'.jpeg')) 
    plot4 <- magick::image_read(paste0('p4',"_",Pep_Name,'.jpeg'))
    plot5 <- magick::image_read(paste0('p5',"_",Pep_Name,'.jpeg'))
    plot6 <- magick::image_read(paste0('p6',"_",Pep_Name,'.jpeg'))
    row1 <- image_append(c(plot1,plot2,plot3))
    row2 <- image_append(c(plot4,plot5,plot6))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "1500")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==7){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))
    plot3 <- magick::image_read(paste0('p3',"_",Pep_Name,'.jpeg')) 
    plot4 <- magick::image_read(paste0('p4',"_",Pep_Name,'.jpeg'))
    plot5 <- magick::image_read(paste0('p5',"_",Pep_Name,'.jpeg'))
    plot6 <- magick::image_read(paste0('p6',"_",Pep_Name,'.jpeg'))
    plot7 <- magick::image_read(paste0('p7',"_",Pep_Name,'.jpeg'))
    row1 <- image_append(c(plot1,plot2,plot3,plot4))
    row2 <- image_append(c(plot5,plot6,plot7))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "2500")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==8){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))
    plot3 <- magick::image_read(paste0('p3',"_",Pep_Name,'.jpeg')) 
    plot4 <- magick::image_read(paste0('p4',"_",Pep_Name,'.jpeg'))
    plot5 <- magick::image_read(paste0('p5',"_",Pep_Name,'.jpeg'))
    plot6 <- magick::image_read(paste0('p6',"_",Pep_Name,'.jpeg'))
    plot7 <- magick::image_read(paste0('p7',"_",Pep_Name,'.jpeg'))
    plot8 <- magick::image_read(paste0('p8',"_",Pep_Name,'.jpeg'))
    row1 <- image_append(c(plot1,plot2,plot3, plot4))
    row2 <- image_append(c(plot5,plot6, plot7, plot8))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "2500")
    plotF<- magick::image_append(c(row1,row2,legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==9){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))
    plot3 <- magick::image_read(paste0('p3',"_",Pep_Name,'.jpeg')) 
    plot4 <- magick::image_read(paste0('p4',"_",Pep_Name,'.jpeg'))
    plot5 <- magick::image_read(paste0('p5',"_",Pep_Name,'.jpeg'))
    plot6 <- magick::image_read(paste0('p6',"_",Pep_Name,'.jpeg'))
    plot7 <- magick::image_read(paste0('p7',"_",Pep_Name,'.jpeg'))
    plot8 <- magick::image_read(paste0('p8',"_",Pep_Name,'.jpeg')) 
    plot9 <- magick::image_read(paste0('p9',"_",Pep_Name,'.jpeg')) 
    row1 <- image_append(c(plot1,plot2,plot3))
    row2 <- image_append(c(plot4,plot5,plot6))
    row3 <- image_append(c(plot7,plot8,plot9))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "1500")
    plotF<- magick::image_append(c(row1,row2,row3, legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==10){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))
    plot3 <- magick::image_read(paste0('p3',"_",Pep_Name,'.jpeg')) 
    plot4 <- magick::image_read(paste0('p4',"_",Pep_Name,'.jpeg'))
    plot5 <- magick::image_read(paste0('p5',"_",Pep_Name,'.jpeg'))
    plot6 <- magick::image_read(paste0('p6',"_",Pep_Name,'.jpeg'))
    plot7 <- magick::image_read(paste0('p7',"_",Pep_Name,'.jpeg'))
    plot8 <- magick::image_read(paste0('p8',"_",Pep_Name,'.jpeg')) 
    plot9 <- magick::image_read(paste0('p9',"_",Pep_Name,'.jpeg'))
    plot10 <- magick::image_read(paste0('p10',"_",Pep_Name,'.jpeg'))
    row1 <- image_append(c(plot1,plot2,plot3,plot4))
    row2 <- image_append(c(plot5,plot6,plot7,plot8))
    row3 <- image_append(c(plot9,plot10))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "2500")
    plotF<- magick::image_append(c(row1,row2,row3, legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==11){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))
    plot3 <- magick::image_read(paste0('p3',"_",Pep_Name,'.jpeg')) 
    plot4 <- magick::image_read(paste0('p4',"_",Pep_Name,'.jpeg'))
    plot5 <- magick::image_read(paste0('p5',"_",Pep_Name,'.jpeg'))
    plot6 <- magick::image_read(paste0('p6',"_",Pep_Name,'.jpeg'))
    plot7 <- magick::image_read(paste0('p7',"_",Pep_Name,'.jpeg'))
    plot8 <- magick::image_read(paste0('p8',"_",Pep_Name,'.jpeg')) 
    plot9 <- magick::image_read(paste0('p9',"_",Pep_Name,'.jpeg'))
    plot10 <- magick::image_read(paste0('p10',"_",Pep_Name,'.jpeg'))
    plot11 <- magick::image_read(paste0('p11',"_",Pep_Name,'.jpeg')) 
    row1 <- image_append(c(plot1,plot2,plot3,plot4))
    row2 <- image_append(c(plot5,plot6,plot7,plot8))
    row3 <- image_append(c(plot9,plot10,plot11))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "2500")
    plotF<- magick::image_append(c(row1,row2,row3, legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  if(length(num_peptime)==12){
    plot1 <- magick::image_read(paste0('p1',"_",Pep_Name,'.jpeg')) 
    plot2 <- magick::image_read(paste0('p2',"_",Pep_Name,'.jpeg'))
    plot3 <- magick::image_read(paste0('p3',"_",Pep_Name,'.jpeg')) 
    plot4 <- magick::image_read(paste0('p4',"_",Pep_Name,'.jpeg'))
    plot5 <- magick::image_read(paste0('p5',"_",Pep_Name,'.jpeg'))
    plot6 <- magick::image_read(paste0('p6',"_",Pep_Name,'.jpeg'))
    plot7 <- magick::image_read(paste0('p7',"_",Pep_Name,'.jpeg'))
    plot8 <- magick::image_read(paste0('p8',"_",Pep_Name,'.jpeg')) 
    plot9 <- magick::image_read(paste0('p9',"_",Pep_Name,'.jpeg'))
    plot10 <- magick::image_read(paste0('p10',"_",Pep_Name,'.jpeg'))
    plot11 <- magick::image_read(paste0('p11',"_",Pep_Name,'.jpeg')) 
    plot12 <- magick::image_read(paste0('p12',"_",Pep_Name,'.jpeg')) 
    row1 <- image_append(c(plot1,plot2,plot3,plot4))
    row2 <- image_append(c(plot5,plot6,plot7,plot8))
    row3 <- image_append(c(plot9,plot10,plot11,plot12))
    legend <- magick::image_read("Legend_1.jpeg")
    legend <- image_scale(legend,"2000")
    legend <- image_scale(legend,"x150")
    legend <- image_border(legend, "white", "2500")
    plotF<- magick::image_append(c(row1,row2,row3, legend), stack=TRUE)
    image_write(plotF, path = paste0("Final_Plot", "_", Pep_Name, ".jpeg"), format = "jpeg", quality = 100)
  }
  
}


###Data Table and Export, Creates and exports the final data tabel
{
  output_total = data.frame()
  for(z in num_peptime){
    temp_df<- (output[[z]])
    output_total <- rbind(output_total,temp_df)
  }
  
  c.fit_total = data.frame()
  for(z in num_peptime[2:length(num_peptime)]){
    temp_df<- (c.fit_tot[[z]])
    c.fit_total <-  qpcR:::rbind.na(c.fit_total,temp_df)
  }    
  
  output_total
  c.fit_total[is.na(c.fit_total)] <- 0
  
  write.csv(output_total, paste0(Pep_Name, "_", "Final_Data.csv"))
  
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
  save_as_image(Final_Data_frame, path =paste0(Pep_Name, "_","Final_Data_Table.jpeg"))
}


###Normalized Area Data
{
  output_tot_env1<- output_total[grep(pattern = "Bimodal Peak 1|Unimodal Peak", x=rownames(output_total)),]
  output_tot_env2 <-  output_total[grep(pattern = "Bimodal Peak 2|Unimodal Peak Env2", x=rownames(output_total)),]
  
  output_tot_env1_rmv <- output_tot_env1[grep(pattern = "Unimodal Peak Env2", x=rownames(output_tot_env1)),]
  output_tot_env1 <- subset(output_tot_env1, !(row.names(output_tot_env1) %in% row.names(output_tot_env1_rmv)))
  
  Uni_Chk <- output_tot_env1[grep(pattern = "Unimodal Peak", x=rownames(output_tot_env1)),]
  
  Uni_Idx <- which((output_tot_env1$`Mean Value`)%in%(Uni_Chk$`Mean Value`))
  
  newrow <- c(NA, NA, NA, NA, NA, NA, NA, 0, NA, NA, NA)
  
  for(i in Uni_Idx){
    if(i==1){
      output_tot_env2 <- rbind(newrow,output_tot_env2)
    }
  }
  
  for(i in Uni_Idx){
    if(i > 1){
      r = i
      newrow_1 <- c(NA, NA, NA, NA, NA, NA, NA, 0, NA, NA, NA)
      output_tot_env2 <- rbind(output_tot_env2[1:(i-1), ],newrow_1, output_tot_env2[- (1:(i-1)), ])
    } 
  }
  
  Uni_Chk2 <- output_tot_env2[grep(pattern = "Unimodal Peak Env2", x=rownames(output_tot_env2)),]
  
  Uni_Idx2 <- which((output_tot_env2$`Mean Value`)%in%(Uni_Chk2$`Mean Value`))
  
  for(j in Uni_Idx2){
    if(j > 1){
      r = j
      newrow_1 <- c(NA, NA, NA, NA, NA, NA, NA, 0, NA, NA, NA)
      output_tot_env1 <- rbind(output_tot_env1[1:(j-1), ],newrow_1, output_tot_env1[- (1:(j-1)), ])
    } 
  }
  
  A_env1 <- (output_tot_env1$`Normalized Area`)
  A_env1 <- A_env1[2:length(A_env1)]
  A_env2 <- (output_tot_env2$`Normalized Area`)
  A_env2 <- A_env2[2:length(A_env2)]
  if(length(time)!=0){
  Area_df <- data.frame(time,A_env1,A_env2)
  Area_df[is.na(Area_df)] <- 0
  colnames(Area_df) <- c("Time (Sec)","Normalized Area of Envelope 1","Normalized Area of Envelope 2")
  }
  if(length(conc)!=0){
    Area_df <- data.frame(conc,A_env1,A_env2)
    Area_df[is.na(Area_df)] <- 0
    colnames(Area_df) <- c("Concentration","Normalized Area of Envelope 1","Normalized Area of Envelope 2")
  }
  if(length(temp)!=0){
    Area_df <- data.frame(temp,A_env1,A_env2)
    Area_df[is.na(Area_df)] <- 0
    colnames(Area_df) <- c("Temperature (K)","Normalized Area of Envelope 1","Normalized Area of Envelope 2")
  }
  if(length(mutant)!=0){
    Area_df <- data.frame(mutant,A_env1,A_env2)
    Area_df[is.na(Area_df)] <- 0
    colnames(Area_df) <- c("Mutant Number","Normalized Area of Envelope 1","Normalized Area of Envelope 2")
  }
  write.csv(Area_df, paste0(Pep_Name,"_","Normalized_Area_Data.csv"))
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
  save_as_image(Final_Det_frame, path =paste0(Pep_Name,"_","Normalized_Area_Data.jpeg"))
  
  y_rel_max<- max(A_env1,A_env2)
  
  
  if(length(time)!=0){
  DA_env1 <- output_tot_env1$`Relative DA`
  DA_env2 <- output_tot_env2$`Relative DA`
  DA_env1 <- DA_env1[2:length(DA_env1)]
  DA_env2 <- DA_env2[2:length(DA_env2)]
  DA_env2[1] <- 0
  DA_env1[1] <- 0
  DA_df <- data.frame(time,DA_env1,DA_env2)
  DA_df$DA_env1 <- na.locf(DA_df$DA_env1)
  DA_df$DA_env2 <- na.locf(DA_df$DA_env2)
  colnames(DA_df) <- c("Time (Sec)","Absolute Deterium Incorperation of Envelope 1","Absolute Deterium Incorperation of Envelope 2")
  write.csv(DA_df, paste0(Pep_Name,"_","Absolute_DA_Data.csv") )
  Final_DA_frame<- DA_df %>%
    as.data.frame() %>% 
    flextable()
  
  fin_bord <- fp_border(color ="dimgray", style="solid", width=2)
  
  Final_DA_frame<- Final_DA_frame %>%
    set_table_properties(layout = "autofit")%>%
    bg(bg="white", part="all")%>%
    hline() %>%
    set_header_labels(Final_DA_frame, rowname = "Timepoint")%>%
    bold(part="header")%>%
    bold(j=1)%>%
    vline(j=1, border = fin_bord )
  
  autofit(Final_DA_frame)
  save_as_image(Final_DA_frame, path =paste0(Pep_Name,"_","Absolute_DA_Data.jpeg"))
  
  y_abi_max <- max(DA_df[,c(2,3)])
  DA_1 <- DA_df[,2]
  DA_2 <- DA_df[,3]
  }
  
  if(length(conc)!=0){
    DA_env1 <- output_tot_env1$`Relative DA`
    DA_env2 <- output_tot_env2$`Relative DA`
    DA_env1 <- DA_env1[2:length(DA_env1)]
    DA_env2 <- DA_env2[2:length(DA_env2)]
    DA_env2[1] <- 0
    DA_env1[1] <- 0
    DA_df <- data.frame(conc,DA_env1,DA_env2)
    DA_df$DA_env1 <- na.locf(DA_df$DA_env1)
    DA_df$DA_env2 <- na.locf(DA_df$DA_env2)
    colnames(DA_df) <- c("Concentration","Absolute Deterium Incorperation of Envelope 1","Absolute Deterium Incorperation of Envelope 2")
    write.csv(DA_df, paste0(Pep_Name,"_","Absolute_DA_Data.csv") )
    Final_DA_frame<- DA_df %>%
      as.data.frame() %>% 
      flextable()
    
    fin_bord <- fp_border(color ="dimgray", style="solid", width=2)
    
    Final_DA_frame<- Final_DA_frame %>%
      set_table_properties(layout = "autofit")%>%
      bg(bg="white", part="all")%>%
      hline() %>%
      set_header_labels(Final_DA_frame, rowname = "Timepoint")%>%
      bold(part="header")%>%
      bold(j=1)%>%
      vline(j=1, border = fin_bord )
    
    autofit(Final_DA_frame)
    save_as_image(Final_DA_frame, path =paste0(Pep_Name,"_","Absolute_DA_Data.jpeg"))
    
    y_abi_max <- max(DA_df[,c(2,3)])
    DA_1 <- DA_df[,2]
    DA_2 <- DA_df[,3]
    
    
  }
  
  if(length(time)!=0){
    layout(mat = matrix(c(1, 2, 3, 3, 4,4),
                        ncol=2, 
                        byrow=TRUE),
           heights = c(10,1.7,0.00000001),
           widths = c(10,10)
    )
    par(mar=c(4,4,1,1))
    layout.show(n=4)
    if(Log==TRUE){
    plot(log(time),A_env1, col="green", type="b", xlab="Log of Time (Sec)", ylab="Normalized Areas", main="Log Linear Graph of Normalized Areas vs Time", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_rel_max+0.1)))
    lines(log(time),A_env2, col="blue", type="b", lwd=2)
    plot(log(time),DA_1, col="green", type="b", xlab="Log of Time (Sec)", ylab="Absolute Deterium Incorperation", main="Log Linear Graph of Absolute Deterium vs Time", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_abi_max+0.1)))
    lines(log(time),DA_2, col="blue", type="b", lwd=2)
    }
    if(Linear==TRUE){
      plot((time),A_env1, col="green", type="b", xlab="Log of Time (Sec)", ylab="Normalized Areas", main="Log Linear Graph of Normalized Areas vs Time", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_rel_max+0.1)))
      lines((time),A_env2, col="blue", type="b", lwd=2)
      plot((time),DA_1, col="green", type="b", xlab="Log of Time (Sec)", ylab="Absolute Deterium Incorperation", main="Log Linear Graph of Absolute Deterium vs Time", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_abi_max+0.1)))
      lines((time),DA_2, col="blue", type="b", lwd=2)
    }
    
    plot(1, type = "n", axes=F, xlab="", ylab="") # Create empty plot
    legend("center", lty = c(1,1), horiz = TRUE, col = c("green", "blue"), c("Envelope 1", "Envelope 2"), lwd = 2, cex=1.1)
    }
  
  if(length(conc)!=0){
    if(Line==TRUE){
        layout(mat = matrix(c(1, 2, 3, 3, 4,4),
                            ncol=2, 
                            byrow=TRUE),
               heights = c(10,1.7,0.00000001),
               widths = c(10,10)
        )
        par(mar=c(4,4,1,1))
        layout.show(n=4)
      if(Log==TRUE){
        plot(log(conc),A_env1, col="green", type="b", xlab="Concentration", ylab="Normalized Areas", main="Graph of Normalized Areas vs Log of Concentration", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_rel_max+0.1)))
        lines(log(conc),A_env2, col="blue", type="b", lwd=2)
        plot(log(conc),DA_1, col="green", type="b", xlab="Concentration", ylab="Absolute Deterium Incorperation", main="Log Linear Graph of Absolute Deterium vs Log of Concentration", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_abi_max+0.1)))
        lines(log(conc),DA_2, col="blue", type="b", lwd=2)
      }
      if(Linear==TRUE){
        plot((conc),A_env1, col="green", type="b", xlab="Concentration", ylab="Normalized Areas", main="Graph of Normalized Areas vs Concentration", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_rel_max+0.1)))
        lines((conc),A_env2, col="blue", type="b", lwd=2)
        plot((conc),DA_1, col="green", type="b", xlab="Concentration", ylab="Absolute Deterium Incorperation", main="Graph of Absolute Deterium vs Concentration", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_abi_max+0.1)))
        lines((conc),DA_2, col="blue", type="b", lwd=2)
      }
      plot(1, type = "n", axes=F, xlab="", ylab="") # Create empty plot
      legend("center", lty = c(1,1), horiz = TRUE, col = c("green", "blue"), c("Envelope 1", "Envelope 2"), lwd = 2, cex=1.1)
    }
    if(Bar==TRUE){
      if(Log==TRUE){
        plot_df <- data.frame(log(conc),A_env1,A_env2)
        colnames(plot_df) <- c("Concentration", "Normalized Area of Envelope 1", "Normalized Area of Enevelope 2")
        plot_df <- melt(plot_df, id="Concentration")
        pc <- ggplot()+
        geom_bar(data=plot_df,aes(x=Concentration,y=value, fill = variable), position = "dodge", stat = "identity")+
          geom_hline(yintercept=0, lwd=1.3)+
          scale_x_continuous(breaks=log(conc), labels=c(round(log(conc),1)))+
          ggtitle(c(print(paste0("Normalized Area vs. Log of Concentration"," ", paste0(Pep_Name))))) +
          theme_tufte()+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
                plot.title = element_text(color="black", size=axis_font, face="bold"),
                axis.title = element_blank(),
                axis.text.x = element_text(face="bold", color="Black",size=axis_font),
                axis.text.y = element_text(face="bold", color="Black",size=axis_font),
                )+
          guides(color = guide_legend(override.aes = list(size = 0.5)),
                 shape = guide_legend(override.aes = list(size = 0.5)),
                 legend.title = element_text(size = 3), 
                 legend.text = element_text(size = 3))
      }
      if(Linear==TRUE){
        plot_df <- data.frame(conc,A_env1,A_env2)
        colnames(plot_df) <- c("Concentration", "Normalized Area of Envelope 1", "Normalized Area of Enevelope 2")
        plot_df <- melt(plot_df, id="Concentration")
        pc <-ggplot()+
          geom_bar(data=plot_df,aes(x=Concentration,y=value, fill = variable), position = "dodge", stat = "identity")+
          geom_hline(yintercept=0, lwd=1.3)+
          scale_x_continuous(breaks=conc, labels=c(round(conc,1)))+
          ggtitle(c(print(paste0("Normalized Area vs. Concentration"," ", paste0(Pep_Name))))) +
          theme_tufte()+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
                plot.title = element_text(color="black", size=axis_font, face="bold"),
                axis.title = element_blank(),
                axis.text.x = element_text(face="bold", color="Black",size=axis_font),
                axis.text.y = element_text(face="bold", color="Black",size=axis_font),
          )+
          guides(color = guide_legend(override.aes = list(size = 0.5)),
                 shape = guide_legend(override.aes = list(size = 0.5)),
                 legend.title = element_text(size = 3), 
                 legend.text = element_text(size = 3))
      }
    plot(pc)
    }
  }
  
  if(length(temp)!=0){
    if(Line==TRUE){
      if(Log==TRUE){
        plot(log(temp),A_env1, col="green", type="b", xlab="Temperature", ylab="Normalized Areas", main="Graph of Normalized Areas vs Log of Temperature", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_rel_max+0.1)))
        lines(log(temp),A_env2, col="blue", type="b", lwd=2)
      }
      if(Linear==TRUE){
        plot((temp),A_env1, col="green", type="b", xlab="Temperature", ylab="Normalized Areas", main="Graph of Normalized Areas vs Temperature", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_rel_max+0.1)))
        lines((temp),A_env2, col="blue", type="b", lwd=2)
      }
    }
    if(Bar==TRUE){
      if(Log==TRUE){
        plot_df <- data.frame(log(temp),A_env1,A_env2)
        colnames(plot_df) <- c("Temperature", "Normalized Area of Envelope 1", "Normalized Area of Enevelope 2")
        plot_df <- melt(plot_df, id="Temperature")
        pt <- ggplot()+
          geom_bar(data=plot_df,aes(x=Temperature,y=value, fill = variable), position = "dodge", stat = "identity")+
          geom_hline(yintercept=0, lwd=1.3)+
          scale_x_continuous(breaks=log(temp), labels=c(round(log(temp),1)))+
          ggtitle(c(print(paste0("Normalized Area vs. Log of Temperature"," ", paste0(Pep_Name))))) +
          theme_tufte()+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
                plot.title = element_text(color="black", size=axis_font, face="bold"),
                axis.title = element_blank(),
                axis.text.x = element_text(face="bold", color="Black",size=axis_font),
                axis.text.y = element_text(face="bold", color="Black",size=axis_font),
          )+
          guides(color = guide_legend(override.aes = list(size = 0.5)),
                 shape = guide_legend(override.aes = list(size = 0.5)),
                 legend.title = element_text(size = 3), 
                 legend.text = element_text(size = 3))
      }
      if(Linear==TRUE){
        plot_df <- data.frame(temp,A_env1,A_env2)
        colnames(plot_df) <- c("Temperature", "Normalized Area of Envelope 1", "Normalized Area of Enevelope 2")
        plot_df <- melt(plot_df, id="Temperature")
        pt <- ggplot()+
          geom_bar(data=plot_df,aes(x=Temperature,y=value, fill = variable), position = "dodge", stat = "identity")+
          geom_hline(yintercept=0, lwd=1.3)+
          scale_x_continuous(breaks=temp, labels=c(round(temp,1)))+
          ggtitle(c(print(paste0("Normalized Area vs. Temperature"," ", paste0(Pep_Name))))) +
          theme_tufte()+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
                plot.title = element_text(color="black", size=axis_font, face="bold"),
                axis.title = element_blank(),
                axis.text.x = element_text(face="bold", color="Black",size=axis_font),
                axis.text.y = element_text(face="bold", color="Black",size=axis_font),
          )+
          guides(color = guide_legend(override.aes = list(size = 0.5)),
                 shape = guide_legend(override.aes = list(size = 0.5)),
                 legend.title = element_text(size = 3), 
                 legend.text = element_text(size = 3))
      }
      plot(pt)
      }
    

    }
  
  if(length(mutant)!=0){
    if(Line==TRUE){
        plot((mutant),A_env1, col="green", type="b", xlab="Mutant", ylab="Normalized Areas", main="Graph of Normalized Areas vs Mutant", lwd=2, font=2, cex.lab=1.2, cex.axis=1.2, ylim=c(0,(y_rel_max+0.1)))
        lines((mutant),A_env2, col="blue", type="b", lwd=2)
      }
    if(Bar==TRUE){
        plot_df <- data.frame(mutant,A_env1,A_env2)
        colnames(plot_df) <- c("Mutant", "Normalized Area of Envelope 1", "Normalized Area of Enevelope 2")
        plot_df <- melt(plot_df, id="Mutant")
        pm <- ggplot()+
          geom_bar(data=plot_df,aes(x=Mutant,y=value, fill = variable), position = "dodge", stat = "identity")+
          geom_hline(yintercept=0, lwd=1.3)+
          scale_x_continuous(breaks=mutant, labels=c(round(mutant,1)))+
          ggtitle(c(print(paste0("Normalized Area vs. Mutant Number"," ", paste0(Pep_Name))))) +
          theme_tufte()+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=border_width),
                plot.title = element_text(color="black", size=axis_font, face="bold"),
                axis.title = element_blank(),
                axis.text.x = element_text(face="bold", color="Black",size=axis_font),
                axis.text.y = element_text(face="bold", color="Black",size=axis_font),
          )+
          guides(color = guide_legend(override.aes = list(size = 0.5)),
                 shape = guide_legend(override.aes = list(size = 0.5)),
                 legend.title = element_text(size = 3), 
                 legend.text = element_text(size = 3))
        
        plot(pm)
    }
    }
  

} #Creates the data tables and plots for the normalized area (for any analysis method) and the Deuteration (if applicable)
{
  dev.copy(jpeg,paste0(Pep_Name, "_","Normalized_Area_Graph.jpeg"), width=3840, height=2020, res=200)
  dev.size(units=c("px"))
  dev.off()
  
  graphics.off()
  
  } #Save Plot, run this after a plot appears in the plots tab
