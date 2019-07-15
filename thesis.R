library(dplyr)
library(ggplot2)
library(stringr)
library(magrittr)
library(lubridate)
library(bdscale)
library(fda)
library(fda.usc)
library(TTR)
library(pastecs)
require(reshape2)
multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}

####---basic processing---####
daily <- read.table("daily.txt", header = T, sep = "\t", stringsAsFactors = F, fill = T)
daily <- daily[, -3]
colnames(daily) <- c("code", "name", "yyyymmdd", "ret")
daily$code <- gsub("\\s+", "", daily$code)
daily$name <- gsub("\\s+", "", daily$name)
daily$yyyymmdd <- as.numeric(daily$yyyymmdd)
daily$ret <- as.numeric(daily$ret)
daily$date <- as.Date(paste0(substring(daily$yyyymmdd, 1, 4), "/", substring(daily$yyyymmdd, 5, 6), "/", substring(daily$yyyymmdd, 7, 8)))
daily$weekday <- wday(daily$date, label = T) 

daily %<>%
  mutate(yyyy = paste0(substring(daily$yyyymmdd, 1, 4)),
         yearday = yday(daily$date),
         nweek = format(date,'%U'))

#descriptive statitics
library(e1071)
library(normtest)
library(tseries)
indList <- paste0(unique(daily$code), " ", unique(daily$name))
codeList <- unique(daily$code)
nameList <- unique(daily$name)
dateList <- as.data.frame(unique(daily$yyyymmdd))
colnames(dateList) <- c("date")

select.ind <- function(i){
  df <- daily %>% filter(code %in% codeList[i])
  
  return(df)
}

mean <- NULL
sd <- NULL
cv <- NULL
skew <- NULL
kurt <- NULL
jb <- NULL
pvalue <- NULL
AC.1 <-NULL
AC.5 <-NULL
AC.20 <-NULL
AC.60 <-NULL
AC.120 <-NULL
AC.252 <-NULL
for(i in 1:length(indList)){
  df <- select.ind(i)
  
  mean[i] <- mean(df$ret)
  sd[i] <- sd(df$ret)
  cv[i] <- mean[i]/sd[i]
  skew[i] <- skewness(df$ret)
  kurt[i] <- kurtosis(df$ret)
  jb.test <- jb.norm.test(df$ret)
  jb[i] <- jb.test$statistic
  pvalue[i] <- jb.test$p.value
  acfList<- acf(df$ret, plot = FALSE, lag.max = 252)$acf
  AC.1[i] <- acfList[2]
  AC.5[i] <- acfList[6]
  AC.20[i] <- acfList[21]
  AC.60[i] <- acfList[61]
  AC.120[i] <- acfList[121]
  AC.252[i] <- acfList[253]
}
d.stat <- data.frame(Variable = indList, Mean = mean, Sd = sd, CV = cv, Skew. = skew, Kurt. = kurt, JB = jb, pvalue = pvalue,
                    AC.1 = AC.1, AC.5 = AC.5, AC.20 = AC.20, AC.60 = AC.60, AC.120 = AC.120, AC.252 = AC.252 ,row.names = NULL)

####---term processing & B-spine mean term structure---####
date.2009 <- NULL
date.2010 <- NULL
date.2011 <- NULL
date.2012 <- NULL
date.2013 <- NULL
date.2014 <- NULL
date.2015 <- NULL
date.2016 <- NULL
date.2017 <- NULL
date.2018 <- NULL

for(i in 1:12){
  selected <- dateList %>% 
    filter(date >= 20090000 + 100 * i)
  date.2009[i] <- min(selected)
  selected <- dateList %>% 
    filter(date >= 20100000 + 100 * i)
  date.2010[i] <- min(selected)
  selected <- dateList %>% 
    filter(date >= 20110000 + 100 * i)
  date.2011[i] <- min(selected)
  selected <- dateList %>% 
    filter(date >= 20120000 + 100 * i)
  date.2012[i] <- min(selected)
  selected <- dateList %>% 
    filter(date >= 20130000 + 100 * i)
  date.2013[i] <- min(selected)
  selected <- dateList %>% 
    filter(date >= 20140000 + 100 * i)
  date.2014[i] <- min(selected)
  selected <- dateList %>% 
    filter(date >= 20150000 + 100 * i)
  date.2015[i] <- min(selected)
  selected <- dateList %>% 
    filter(date >= 20160000 + 100 * i)
  date.2016[i] <- min(selected)
  selected <- dateList %>% 
    filter(date >= 20170000 + 100 * i)
  date.2017[i] <- min(selected)
  selected <- dateList %>% 
    filter(date >= 20180000 + 100 * i)
  date.2018[i] <- min(selected)
}
firstDate <- c(date.2009, date.2010, date.2011, date.2012, date.2013, date.2014, date.2015, date.2016, date.2017, date.2018)

for(i in 1:12){
  selected <- dateList %>% 
    filter(date >= 20090000 + 100 * i & date <= 20090000 + 100 * (i + 1))  
  date.2009[i] <- nrow(selected)
  selected <- dateList %>% 
    filter(date >= 20100000 + 100 * i & date <= 20100000 + 100 * (i + 1))  
  date.2010[i] <- nrow(selected)
  selected <- dateList %>% 
    filter(date >= 20110000 + 100 * i & date <= 20110000 + 100 * (i + 1))  
  date.2011[i] <- nrow(selected)
  selected <- dateList %>% 
    filter(date >= 20120000 + 100 * i & date <= 20120000 + 100 * (i + 1))  
  date.2012[i] <- nrow(selected)
  selected <- dateList %>% 
    filter(date >= 20130000 + 100 * i & date <= 20130000 + 100 * (i + 1))  
  date.2013[i] <- nrow(selected)
  selected <- dateList %>% 
    filter(date >= 20140000 + 100 * i & date <= 20140000 + 100 * (i + 1))  
  date.2014[i] <- nrow(selected)
  selected <- dateList %>% 
    filter(date >= 20150000 + 100 * i & date <= 20150000 + 100 * (i + 1))  
  date.2015[i] <- nrow(selected)
  selected <- dateList %>% 
    filter(date >= 20160000 + 100 * i & date <= 20160000 + 100 * (i + 1))  
  date.2016[i] <- nrow(selected)
  selected <- dateList %>% 
    filter(date >= 20170000 + 100 * i & date <= 20170000 + 100 * (i + 1))  
  date.2017[i] <- nrow(selected)
  selected <- dateList %>% 
    filter(date >= 20180000 + 100 * i & date <= 20180000 + 100 * (i + 1))  
  date.2018[i] <- nrow(selected)
}
monthLength <- c(date.2009, date.2010, date.2011, date.2012, date.2013, date.2014, date.2015, date.2016, date.2017, date.2018)


daily %<>% 
  group_by(code) %>% 
  mutate(ret.1 = SMA(ret, 5) * 250, ret.2 = SMA(ret, 10) * 250, ret.3 = SMA(ret, 20) * 250, 
         ret.4 = SMA(ret, 30) * 250, ret.5 = SMA(ret, 60) * 250, ret.6 = SMA(ret, 90) * 250,
         ret.7 = SMA(ret, 120) * 250, ret.8 = SMA(ret, 180) * 250, ret.9 = SMA(ret, 250) * 250, 
         ret.10 = SMA(ret, 120) * 250, ret.11 = SMA(ret, 180) * 250, ret.12 = SMA(ret, 250) * 250)

fang <- function(y, m){
  n <- max(which(firstDate <= (y * 10000 + (m + 1) * 100)))

  results <- daily %>% 
    filter(yyyymmdd <= (y * 10000 + (m + 1) * 100)) %>%
    group_by(code) %>%
    mutate(ret.1 = SMA(ret, sum(monthLength[(n-1):(n-1)])) * 250, ret.2 = SMA(ret, sum(monthLength[(n-2):(n-1)])) * 250,
           ret.3 = SMA(ret, sum(monthLength[(n-3):(n-1)])) * 250, ret.4 = SMA(ret, sum(monthLength[(n-4):(n-1)])) * 250, 
           ret.5 = SMA(ret, sum(monthLength[(n-5):(n-1)])) * 250, ret.6 = SMA(ret, sum(monthLength[(n-6):(n-1)])) * 250,
           ret.7 = SMA(ret, sum(monthLength[(n-7):(n-1)])) * 250, ret.8 = SMA(ret, sum(monthLength[(n-8):(n-1)])) * 250, 
           ret.9 = SMA(ret, sum(monthLength[(n-9):(n-1)])) * 250, ret.10 = SMA(ret, sum(monthLength[(n-10):(n-1)])) * 250, 
           ret.11 = SMA(ret, sum(monthLength[(n-11):(n-1)])) * 250, ret.12 = SMA(ret, sum(monthLength[(n-12):(n-1)])) * 250) %>%
    filter(yyyymmdd == firstDate[n]) %>%
    dplyr::select(code, name, ret.1, ret.2, ret.3, ret.4, ret.5, ret.6, ret.7, ret.8, ret.9, ret.10, ret.11, ret.12)
  
  return(results)
}  

result <- matrix(nrow = 108, ncol = 12)

fang.plot.raw <- function(i){
  for(y in 1:9){
    for(m in 1:12){
      hi <- fang((y + 2009), m) %>%
        filter(code == codeList[i])
      result[12*(y-1)+m, 1:12] <- unlist(hi[, 3:14])
    }
  }
  ggplot(data = melt(result), aes(x = Var2, y = value, group=factor(Var1))) +
    geom_line(aes(color=factor(Var1))) +
    ggtitle(indList[i]) +
    labs(x = "Frequncy", y = "Return(%)") +
    scale_x_continuous(breaks = seq(1, 12, 1)) +
    theme_light() +
    theme(legend.position="none",
          plot.title = element_text(size=12, hjust = 0.5, face = "bold")) 
}
fang.plot.0 <- function(i){
  for(y in 1:9){
    for(m in 1:12){
      hi <- fang((y + 2009), m) %>%
        filter(code == codeList[i])
      result[12*(y-1)+m, 1:12] <- unlist(hi[, 3:14])
    }
  }
  ret.fdata <- fdata(result)
  ret.fd <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 12)
  p <- plot(ret.fd, xlab = "Frequency (month)", ylab = "Return (%)", main = indList[i])
  
  return(p)
}
fang.plot.1 <- function(i){
  for(y in 1:9){
    for(m in 1:12){
      hi <- fang((y + 2009), m) %>%
        filter(code == codeList[i])
      result[12*(y-1)+m, 1:12] <- unlist(hi[, 3:14])
    }
  }
  ret.fdata <- fdata(result)
  ret.fd <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 12, nderiv = 1)
  p <- plot(ret.fd, xlab = "Frequency (month)", ylab = "Return (%)", main = indList[i])
  
  return(p)
}
fang.plot.2 <- function(i){
  for(y in 1:9){
    for(m in 1:12){
      hi <- fang((y + 2009), m) %>%
        filter(code == codeList[i])
      result[12*(y-1)+m, 1:12] <- unlist(hi[, 3:14])
    }
  }
  ret.fdata <- fdata(result)
  ret.fd <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 12, nderiv = 2)
  p <- plot(ret.fd, xlab = "Frequency (month)", ylab = "Return (%)", main = indList[i])
  
  return(p)
}
fang.plot.mean <- function(i){
  for(y in 1:9){
    for(m in 1:12){
      hi <- fang((y + 2009), m) %>%
        filter(code == codeList[i])
      result[12*(y-1)+m, 1:12] <- unlist(hi[, 3:14])
    }
  }
  ret.fdata <- fdata(result)
  ret.fd <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 15)
  meanFd <- mean.fd(ret.fd)
  plot(meanFd, lwd = 3, xlab = "Frequency (month)", ylab = "Return (%)", main = indList[i])
  abline(h = 0, lty = 2, lwd = 1, col="red")
}
fang.plot.mean.1 <- function(i){
  for(y in 1:9){
    for(m in 1:12){
      hi <- fang((y + 2009), m) %>%
        filter(code == codeList[i])
      result[12*(y-1)+m, 1:12] <- unlist(hi[, 3:14])
    }
  }
  ret.fdata <- fdata(result)
  ret.fd <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 15, nderiv = 1)
  meanFd <- mean.fd(ret.fd)
  plot(meanFd, lwd = 3, xlab = "Frequency (month)", ylab = "Slope", main = indList[i])
  abline(h = 0, lty = 2, lwd = 1, col="red")
}
fang.plot.mean.2 <- function(i){
  for(y in 1:9){
    for(m in 1:12){
      hi <- fang((y + 2009), m) %>%
        filter(code == codeList[i])
      result[12*(y-1)+m, 1:12] <- unlist(hi[, 3:14])
    }
  }
  ret.fdata <- fdata(result)
  ret.fd <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 15, nderiv = 2)
  meanFd <- mean.fd(ret.fd)
  plot(meanFd, lwd = 3, xlab = "Frequency (month)", ylab = "Curvature", main = indList[i])
  abline(h = 0, lty = 2, lwd = 1, col="red")
}
fang.plot.mean.mix <- function(i){
  for(y in 1:9){
    for(m in 1:12){
      hi <- fang((y + 2009), m) %>%
        filter(code == codeList[i])
      result[12*(y-1)+m, 1:12] <- unlist(hi[, 3:14])
    }
  }
  ret.fdata <- fdata(result)
  ret.fd.1 <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 15, nderiv = 1)
  ret.fd.2 <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 15, nderiv = 2)
  meanFd.1 <- mean.fd(ret.fd.1)
  meanFd.2 <- mean.fd(ret.fd.2)
  plot(meanFd.2, lwd = 3, xlab = "Frequency (month)", ylab = "", main = indList[i])
  lines(meanFd.1, lwd = 3, col = "blue")
  abline(h = 0, lty = 2, lwd = 1, col="red")
}


for(i in 1:length(indList)){
  assign(paste0("p", i), fang.plot.raw(i)) 
}
gridExtra::grid.arrange(p1, p2, p3, p4, p5,
                        p6, p7, p8, p9, p10,
                        p11, p12, p13, p14, p15,
                        p16, p17, p18, p19, p20,
                        p21, p22, p23, p24, p25,
                        p26, p27, p28, p29, p30,
                        nrow = 10)

for(i in 1:length(indList)){
  png(filename = paste0(indList[i], "0.png"), width=650, height=500)
  fang.plot.0(i)
  dev.off()
}
for(i in 1:length(indList)){
  png(filename = paste0(indList[i], "1.png"), width=650, height=500)
  fang.plot.1(i)
  dev.off()
}
for(i in 1:length(indList)){
  png(filename = paste0(indList[i], "2.png"), width=650, height=500)
  fang.plot.2(i)
  dev.off()
}
for(i in 1:length(indList)){
  png(filename = paste0(indList[i], "raw.png"), width=650, height=500)
  print(fang.plot.raw(i))
  dev.off()
}
for(i in 1:length(indList)){
  png(filename = paste0(indList[i], "mean.png"), width=650, height=500)
  fang.plot.mean(i)
  dev.off()
}
for(i in 1:length(indList)){
  png(filename = paste0(indList[i], "mean1.png"), width=650, height=500)
  fang.plot.mean.1(i)
  dev.off()
}
for(i in 1:length(indList)){
  png(filename = paste0(indList[i], "mean2.png"), width=650, height=500)
  fang.plot.mean.2(i)
  dev.off()
}
for(i in 1:length(indList)){
  png(filename = paste0(indList[i], "mix.png"), width=650, height=500)
  fang.plot.mean.mix(i)
  dev.off()
}

png("fang_plot_0.png", width=1464, height=2380)
par(mfrow = c(10, 3))
par(mar = c(1,1,1,1))
for(i in 1:length(indList)){
  fang.plot.0(i)
}
dev.off()

png("fang_plot_1.png", width=1464, height=2380)
par(mfrow = c(10, 3))
par(mar = c(1,1,1,1))
for(i in 1:length(indList)){
  fang.plot.1(i) 
}
dev.off()

png("fang_plot_2.png", width=1464, height=2380)
par(mfrow = c(10, 3))
par(mar = c(1,1,1,1))
for(i in 1:length(indList)){
  fang.plot.2(i) 
}
dev.off()

png("fang_plot_m1.png", width=1464, height=2380)
par(mfrow = c(10, 3))
par(mar = c(1,1,1,1))
for(i in 1:length(indList)){
  fang.plot.mean.1(i) 
}
dev.off()

png("fang_plot_m2.png", width=1464, height=2380)
par(mfrow = c(10, 3))
par(mar = c(1,1,1,1))
for(i in 1:length(indList)){
  fang.plot.mean.2(i) 
}
dev.off()

####---using functional model---####
reg <- as.data.frame(matrix(nrow = 917, ncol = 2))
colnames(reg) <- c("q", "r2")
fm.results <- as.data.frame(matrix(nrow = 108, ncol = 8))
colnames(fm.results) <- c("q", "c0", "p.c0", "c1", "p.c1", "c2", "p.c2", "r2")
optimal <- NULL

fang.fm <- function(i){
  for(y in 1:9){
    for(m in 1:12){
      hi <- fang((y + 2009), m) %>%
        filter(code == codeList[i])
      result[12*(y-1)+m, 1:12] <- unlist(hi[, 3:14])
    }
  }
  for(j in 1:108){
    for(i in seq(1/12, 1, 0.001)){
      #find optimal q
      q <- i
      t <- seq(1/12, 1 , 1/12)
      f1.t <- exp(-t/q)
      f2.t <- (t/q)*exp(1-t/q)
      
      func <- lm(result[j, ] ~ f1.t + f2.t)
      reg$q[round(i*1000)-82] <- q
      reg$r2[round(i*1000)-82] <- summary(func)$r.squared
      
      reg %<>% arrange(desc(r2))
      optimal <- reg[1, ]
      
      #record it
      q <- optimal$q
      t <- seq(1/12, 1 , 1/12)
      f1.t <- exp(-t/q)
      f2.t <- (t/q)*exp(1-t/q)
      
      func <- lm(result[1, ] ~ f1.t + f2.t)
      fm.results$q[j] <- q
      fm.results$c0[j] <- summary(func)$coefficients[1, 1]
      fm.results$p.c0[j] <- summary(func)$coefficients[1, 4]
      fm.results$c1[j] <- summary(func)$coefficients[2, 1]
      fm.results$p.c1[j] <- summary(func)$coefficients[2, 4]
      fm.results$c2[j] <- summary(func)$coefficients[3, 1]
      fm.results$p.c2[j] <- summary(func)$coefficients[3, 4]
      fm.results$r2[j] <- summary(func)$r.squared
    }
  }
  return(fm.results)
}

#save statistics results
for(i in 2:length(indList)){
  fm.results <- data.frame(fang.fm(i))
  write.csv(fm.results, file = paste0(indList[i], ".csv"), row.names = F)
}

#load statistics results
for(i in 1:length(indList)){
  assign(paste0(codeList[i]), read.csv(paste0(indList[i], ".csv"))) 
}

####---distribution---####
distribution.q <- function(df){
  list1 <- c(mean(df$q), sd(df$q), quantile(df$q, 0.25), quantile(df$q, 0.5), quantile(df$q, 0.75))
  return(list1)
}
distribution.r2 <- function(df){
  list1 <- c(mean(df$r2), sd(df$r2), quantile(df$r2, 0.25), quantile(df$r2, 0.5), quantile(df$r2, 0.75))
  return(list1)
}
distribution.c0 <- function(df){
  list1 <- c(mean(df$c0), sd(df$c0), quantile(df$c0, 0.25), quantile(df$c0, 0.5), quantile(df$c0, 0.75), sum(df$p.c0 <= 0.05))
  return(list1)
}
distribution.c1 <- function(df){
  list1 <- c(mean(df$c1), sd(df$c1), quantile(df$c1, 0.25), quantile(df$c1, 0.5), quantile(df$c1, 0.75), sum(df$p.c1 <= 0.05))
  return(list1)
}
distribution.c2 <- function(df){
  list1 <- c(mean(df$c2), sd(df$c2), quantile(df$c2, 0.25), quantile(df$c2, 0.5), quantile(df$c2, 0.75), sum(df$p.c2 <= 0.05))
  return(list1)
}

df1 <- NULL
for(i in 1:length(codeList)){
  list1 <- distribution.q(get(codeList[i]))
  df1 <- rbind(df1, list1)
}
df.q <- df1
colnames(df.q) <- c("Mean", "Std.", "Q1", "Median", "Q3")
row.names(df.q) <- indList

#duration
fang.duration <- function(df){
  q <- quantile(df$q, 0.5)
  f0.t <- -t
  f1.t <- exp(-t/q)*-t
  f2.t <- (t/q)*exp(1-t/q)*-t
  duration <- cbind(f0.t, f1.t, f2.t)
  
  return(duration)
}
duration.ind <- as.data.frame(matrix(nrow = 12, ncol = 90))
colnames(duration.ind) <- rep(1:3, 30)
for(i in 1:30){
  duration.ind[, (3*i-2):(3*i)] <- fang.duration(get(codeList[i]))
}
write.csv(duration.ind, "duration.csv")

#q
df1 <- NULL
for(i in 1:length(codeList)){
  list1 <- get(codeList[i])$q
  df1 <- cbind(df1, list1)
}
plot.q <- as.data.frame(df1)
colnames(plot.q) <- unique(daily$name)

png("distribution of q.png", width=850, height=600)
ggplot(data = melt(plot.q), aes(x = reorder(variable, value, FUN = median, desc = T), y = value)) +
  geom_boxplot(aes(fill=variable)) +
  theme(legend.position="none",
        plot.title = element_text(size=20, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

#r2
df1 <- NULL
for(i in 1:length(codeList)){
  list1 <- get(codeList[i])$r2
  df1 <- cbind(df1, list1)
}
plot.r2 <- as.data.frame(df1)
colnames(plot.r2) <- unique(daily$name)

png("distribution of q.png", width=850, height=600)
ggplot(data = melt(plot.r2), aes(x = reorder(variable, value, FUN = median, desc = T), y = value)) +
  geom_boxplot(aes(fill=variable)) +
  theme(legend.position="none",
        plot.title = element_text(size=20, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

#c0
df1 <- NULL
for(i in 1:length(codeList)){
  list1 <- get(codeList[i])$c0
  df1 <- cbind(df1, list1)
}
plot.c0 <- as.data.frame(df1)
colnames(plot.c0) <- unique(daily$name)

require(reshape2)

png("distribution of c0.png", width=850, height=600)
ggplot(data = melt(plot.c0), aes(x = reorder(variable, value, FUN = median, desc = T), y = value)) +
  geom_boxplot(aes(fill=variable)) +
  theme(legend.position="none",
        plot.title = element_text(size=20, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

#c1
df1 <- NULL
for(i in 1:length(codeList)){
  list1 <- get(codeList[i])$c1
  df1 <- cbind(df1, list1)
}
plot.c1 <- as.data.frame(df1)
colnames(plot.c1) <- unique(daily$name)

require(reshape2)

png("distribution of c1.png", width=850, height=600)
ggplot(data = melt(plot.c1), aes(x = reorder(variable, value, FUN = median, desc = T), y = value)) +
  geom_boxplot(aes(fill=variable)) +
  theme(legend.position="none",
        plot.title = element_text(size=20, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

#c2
df1 <- NULL
for(i in 1:length(codeList)){
  list1 <- get(codeList[i])$c2
  df1 <- cbind(df1, list1)
}
plot.c2 <- as.data.frame(df1)
colnames(plot.c2) <- unique(daily$name)

require(reshape2)

png("distribution of c2.png", width=850, height=600)
ggplot(data = melt(plot.c2), aes(x = reorder(variable, value, FUN = median, desc = T), y = value)) +
  geom_boxplot(aes(fill=variable)) +
  theme(legend.position="none",
        plot.title = element_text(size=20, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

