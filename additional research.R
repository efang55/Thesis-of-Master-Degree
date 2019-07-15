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

####---flm---####
fang.ver2 <- function(y, m){
  n <- max(which(firstDate <= (y * 10000 + (m + 1) * 100)))
  
  results <- daily %>% 
    filter(yyyymmdd <= (y * 10000 + (m + 1) * 100)) %>%
    group_by(code) %>%
    mutate(ret.1 = SMA(ret, sum(monthLength[(n-1):(n-1)])) * 250, ret.2 = SMA(ret, sum(monthLength[(n-2):(n-1)])) * 250,
           ret.3 = SMA(ret, sum(monthLength[(n-3):(n-1)])) * 250, ret.4 = SMA(ret, sum(monthLength[(n-4):(n-1)])) * 250, 
           ret.5 = SMA(ret, sum(monthLength[(n-5):(n-1)])) * 250, ret.6 = SMA(ret, sum(monthLength[(n-6):(n-1)])) * 250,
           ret.7 = SMA(ret, sum(monthLength[(n-7):(n-1)])) * 250, ret.8 = SMA(ret, sum(monthLength[(n-8):(n-1)])) * 250, 
           ret.9 = SMA(ret, sum(monthLength[(n-9):(n-1)])) * 250, ret.10 = SMA(ret, sum(monthLength[(n-10):(n-1)])) * 250, 
           ret.11 = SMA(ret, sum(monthLength[(n-11):(n-1)])) * 250, ret.12 = SMA(ret, sum(monthLength[(n-12):(n-1)])) * 250,
           ret.13 = SMA(ret, sum(monthLength[(n):(n)])) * 250) %>%
    filter(yyyymmdd == firstDate[n]) %>%
    dplyr::select(code, name, ret.1, ret.2, ret.3, ret.4, ret.5, ret.6, ret.7, ret.8, ret.9, ret.10, ret.11, ret.12, ret.13)
  
  return(results)
}  
result <- matrix(nrow = 108, ncol = 13)
fang.flm <- function(i){
  for(y in 1:9){
    for(m in 1:12){
      hi <- fang.ver2((y + 2009), m) %>%
        filter(code == codeList[i])
      result[12*(y-1)+m, 1:13] <- unlist(hi[, 3:15])
    }
  }
  ret.fdata <- fdata(result[, -13])
  fd.raw <- as.data.frame(result)
  tt <- ret.fdata[["argvals"]]
  nbasis.x <- 11
  nbasis.b <- 11
  basis1 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.x, norder = 5)
  basis2 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.b, norder = 5)
  basis.x <- list("x" = basis1)
  basis.b <- list("x" = basis2)
  ldata <- list("df" = fd.raw, "x" = ret.fdata)
  
  res <- fregre.lm(V13 ~ x, ldata, basis.x=basis.x, basis.b=basis.b)
  
  g <- summary(res)
  
  return(g)
}
flm.result <- as.data.frame(matrix(nrow = 60, ncol = 13))
for(i in 1:length(codeList)){
  flm.list <- fang.flm(i)
  flm.result[(2*i-1), 1:12] <- round(flm.list$coefficients[, 1], 3)
  flm.result[(2*i), 1:12] <- round(flm.list$coefficients[, 4], 3)
  flm.result[(2*i-1), 13] <- round(flm.list$adj.r.squared, 3)
}


#winsorized
daily %<>%
  filter(ret <= quantile(ret, 0.75) + 1.5*IQR(ret) & ret >= quantile(ret, 0.25) - 1.5*IQR(ret)) %>%
  arrange(code, yyyymmdd)

#raw to fd
par(mfrow = c(10, 3))
par(mar = c(1,1,1,1))
for(i in 1:length(indList)){
  df <- daily %>% 
    filter(code == codeList[i])
  mon <- df %>% 
    filter(weekday == "Mon") %>%
    dplyr::select(yyyy, nweek, ret)
  tue <- df %>% 
    filter(weekday == "Tues") %>%
    dplyr::select(ret, nweek, yyyy)
  wed <- df %>% 
    filter(weekday == "Wed") %>%
    dplyr::select(ret, nweek, yyyy)
  thu <- df %>% 
    filter(weekday == "Thurs") %>%
    dplyr::select(ret, nweek, yyyy)
  fri <- df %>% 
    filter(weekday == "Fri") %>%
    dplyr::select(ret, nweek, yyyy)
  
  fd.raw <- left_join(mon, tue, by = c("nweek", "yyyy")) %>%
    left_join(., wed, by = c("nweek", "yyyy")) %>%
    left_join(., thu, by = c("nweek", "yyyy")) %>%
    left_join(., fri, by = c("nweek", "yyyy"))
  colnames(fd.raw) <- c("yyyy", "nweek", "mon", "tue", "wed", "thu", "fri")
  fd.raw <- as.matrix(fd.raw[, 3:7])
  
  ret.fdata <- fdata(fd.raw)
  ret.fd <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 15)
  plot(ret.fd, xlab = "Weekday", ylab = "Return (%)", main = indList[i])
}

#d1
for(i in 1:length(indList)){
  df <- daily %>% 
    filter(code == codeList[i] & yyyy >= 2016)
  mon <- df %>% 
    filter(weekday == "Mon") %>%
    dplyr::select(yyyy, nweek, ret)
  tue <- df %>% 
    filter(weekday == "Tues") %>%
    dplyr::select(ret, nweek, yyyy)
  wed <- df %>% 
    filter(weekday == "Wed") %>%
    dplyr::select(ret, nweek, yyyy)
  thu <- df %>% 
    filter(weekday == "Thurs") %>%
    dplyr::select(ret, nweek, yyyy)
  fri <- df %>% 
    filter(weekday == "Fri") %>%
    dplyr::select(ret, nweek, yyyy)
  
  fd.raw <- left_join(mon, tue, by = c("nweek", "yyyy")) %>%
    left_join(., wed, by = c("nweek", "yyyy")) %>%
    left_join(., thu, by = c("nweek", "yyyy")) %>%
    left_join(., fri, by = c("nweek", "yyyy"))
  colnames(fd.raw) <- c("yyyy", "nweek", "mon", "tue", "wed", "thu", "fri")
  fd.raw <- as.matrix(fd.raw[, 3:7])
  
  ret.fdata <- fdata(fd.raw)
  ret.fd <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 15, nderiv = 1)
  plot(ret.fd, xlab = "Weekday", ylab = "d1(Return (%))", main = indList[i])
}

#d2
for(i in 1:length(indList)){
  df <- daily %>% 
    filter(code == codeList[i] & yyyy >= 2016)
  mon <- df %>% 
    filter(weekday == "Mon") %>%
    dplyr::select(yyyy, nweek, ret)
  tue <- df %>% 
    filter(weekday == "Tues") %>%
    dplyr::select(ret, nweek, yyyy)
  wed <- df %>% 
    filter(weekday == "Wed") %>%
    dplyr::select(ret, nweek, yyyy)
  thu <- df %>% 
    filter(weekday == "Thurs") %>%
    dplyr::select(ret, nweek, yyyy)
  fri <- df %>% 
    filter(weekday == "Fri") %>%
    dplyr::select(ret, nweek, yyyy)
  
  fd.raw <- left_join(mon, tue, by = c("nweek", "yyyy")) %>%
    left_join(., wed, by = c("nweek", "yyyy")) %>%
    left_join(., thu, by = c("nweek", "yyyy")) %>%
    left_join(., fri, by = c("nweek", "yyyy"))
  colnames(fd.raw) <- c("yyyy", "nweek", "mon", "tue", "wed", "thu", "fri")
  fd.raw <- as.matrix(fd.raw[, 3:7])
  
  ret.fdata <- fdata(fd.raw)
  ret.fd <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 15, nderiv = 2)
  plot(ret.fd, xlab = "Weekday", ylab = "d2(Return (%))", main = indList[i])
}

#mean
for(i in 1:length(indList)){
  df <- daily %>% 
    filter(code == codeList[i])
  mon <- df %>% 
    filter(weekday == "Mon") %>%
    dplyr::select(yyyy, nweek, ret)
  tue <- df %>% 
    filter(weekday == "Tues") %>%
    dplyr::select(ret, nweek, yyyy)
  wed <- df %>% 
    filter(weekday == "Wed") %>%
    dplyr::select(ret, nweek, yyyy)
  thu <- df %>% 
    filter(weekday == "Thurs") %>%
    dplyr::select(ret, nweek, yyyy)
  fri <- df %>% 
    filter(weekday == "Fri") %>%
    dplyr::select(ret, nweek, yyyy)
  
  fd.raw <- left_join(mon, tue, by = c("nweek", "yyyy")) %>%
    left_join(., wed, by = c("nweek", "yyyy")) %>%
    left_join(., thu, by = c("nweek", "yyyy")) %>%
    left_join(., fri, by = c("nweek", "yyyy"))
  colnames(fd.raw) <- c("yyyy", "nweek", "mon", "tue", "wed", "thu", "fri")
  fd.raw <- as.matrix(fd.raw[, 3:7]) %>% 
    na.omit()
  
  ret.fdata <- fdata(fd.raw)
  ret.fd <- fdata2fd(ret.fdata, type.basis="bspline", nbasis = 5)
  meanFd <- mean.fd(ret.fd)
  plot(meanFd, lwd = 3, xlab = "Weekday", ylab = "Return (%)", main = indList[i])
  abline(h = 0, lty = 2, lwd = 1, col="red")
}
dev.off()

#lm.test
lm.result <- matrix(nrow = 60, ncol = 6)
for(i in 1:length(indList)){
  df <- daily %>% 
    filter(code == codeList[i] & weekday %not in% "Sat")
  df$mon <- ifelse(df$weekday %in% "Mon", 1, 0)
  df$tue <- ifelse(df$weekday %in% "Tues", 1, 0)
  df$wed <- ifelse(df$weekday %in% "Wed", 1, 0)
  df$thu <- ifelse(df$weekday %in% "Thurs", 1, 0)
  df.raw <- df %>%
    dplyr::select(ret, mon, tue, wed, thu) %>% 
    na.omit()
  
  reg <- lm(ret ~ mon + tue + wed + thu, df.raw)
  lm.result[(2*i-1):(2*i), 1] <- rep(indList[i], 2)
  lm.result[(2*i-1), 2:6] <- summary(reg)$coefficients[c(2:5, 1), 1]
  lm.result[(2*i), 2:6] <- summary(reg)$coefficients[c(2:5, 1), 4]
}
lm.result <- as.data.frame(lm.result)
colnames(lm.result) <- c("Industry", "Mon", "Tue", "Wed", "Thu", "Fri")

#flm-Friday
yearList <- unique(daily$yyyy)
flm.result <- matrix(nrow = 60, ncol = 41)
for(i in 1:length(indList)){
  for(j in 1:length(yearList)){
    df <- daily %>% 
      filter(code == codeList[i] & yyyy == yearList[j])
    mon <- df %>% 
      filter(weekday == "Mon") %>%
      dplyr::select(yyyy, nweek, ret)
    tue <- df %>% 
      filter(weekday == "Tues") %>%
      dplyr::select(ret, nweek, yyyy)
    wed <- df %>% 
      filter(weekday == "Wed") %>%
      dplyr::select(ret, nweek, yyyy)
    thu <- df %>% 
      filter(weekday == "Thurs") %>%
      dplyr::select(ret, nweek, yyyy)
    fri <- df %>% 
      filter(weekday == "Fri") %>%
      dplyr::select(ret, nweek, yyyy)
    
    fd.raw <- left_join(mon, tue, by = c("nweek", "yyyy")) %>%
      left_join(., wed, by = c("nweek", "yyyy")) %>%
      left_join(., thu, by = c("nweek", "yyyy")) %>%
      left_join(., fri, by = c("nweek", "yyyy"))
    colnames(fd.raw) <- c("yyyy", "nweek", "mon", "tue", "wed", "thu", "fri")
    fd.raw <- as.matrix(fd.raw[, 3:7]) %>% na.omit()
    fd.raw.flm <- fd.raw[, -5]
    ret.fdata <- fdata(fd.raw.flm)
    fd.raw <- as.data.frame(fd.raw) %>% na.omit()
    
    tt <- ret.fdata[["argvals"]]
    nbasis.x <- 4
    nbasis.b <- 4
    basis1 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.x, norder = 4)
    basis2 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.b, norder = 4)
    basis.x <- list("x" = basis1)
    basis.b <- list("x" = basis2)
    ldata <- list("df" = fd.raw, "x" = ret.fdata)
    
    res <- fregre.lm(fri ~ x, ldata, basis.x=basis.x, basis.b=basis.b)
    
    plot(res$beta.l$x, main = paste0(indList[i], "(", yearList[j], ")"))
    flm.result[(i*2-1), (4*j-3):(4*j)] <- summary(res)$coefficients[2:5, 1]
    flm.result[(i*2), (4*j-3):(4*j)] <- summary(res)$coefficients[2:5, 4]
  }
  flm.result[(2*i-1):(2*i), 41] <- rep(indList[i], 2)
}
colnames(flm.result) <- c(rep(c("Mon", "Tue", "Wed", "Thu"), 10), "Industry")

#flm-Monday
yearList <- unique(daily$yyyy)
flm.result <- matrix(nrow = 60, ncol = 51)
for(i in 1:length(indList)){
  for(j in 1:length(yearList)){
    df <- daily %>% 
      filter(code == codeList[i] & yyyy == yearList[j])
    mon <- df %>% 
      filter(weekday == "Mon") %>%
      dplyr::select(yyyy, nweek, ret)
    tue <- df %>% 
      filter(weekday == "Tues") %>%
      dplyr::select(ret, nweek, yyyy)
    wed <- df %>% 
      filter(weekday == "Wed") %>%
      dplyr::select(ret, nweek, yyyy)
    thu <- df %>% 
      filter(weekday == "Thurs") %>%
      dplyr::select(ret, nweek, yyyy)
    fri <- df %>% 
      filter(weekday == "Fri") %>%
      dplyr::select(ret, nweek, yyyy)
    
    fd.raw <- left_join(mon, tue, by = c("nweek", "yyyy")) %>%
      left_join(., wed, by = c("nweek", "yyyy")) %>%
      left_join(., thu, by = c("nweek", "yyyy")) %>%
      left_join(., fri, by = c("nweek", "yyyy"))
    colnames(fd.raw) <- c("yyyy", "nweek", "mon", "tue", "wed", "thu", "fri")
    fd.raw %<>% mutate(lag.mon = lag(mon))
    fd.raw <- as.matrix(fd.raw[, -c(1:2)]) %>% na.omit()
    fd.raw.flm <- fd.raw[, -6]
    ret.fdata <- fdata(fd.raw.flm)
    fd.raw <- as.data.frame(fd.raw) %>% na.omit()
    
    tt <- ret.fdata[["argvals"]]
    nbasis.x <- 5
    nbasis.b <- 5
    basis1 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.x, norder = 5)
    basis2 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.b, norder = 5)
    basis.x <- list("x" = basis1)
    basis.b <- list("x" = basis2)
    ldata <- list("df" = fd.raw, "x" = ret.fdata)
    
    res <- fregre.lm(lag.mon ~ x, ldata, basis.x=basis.x, basis.b=basis.b)
    
    plot(res$beta.l$x, main = paste0(indList[i], "(", yearList[j], ")"))
    flm.result[(i*2-1), (5*j-4):(5*j)] <- summary(res)$coefficients[2:6, 1]
    flm.result[(i*2), (5*j-4):(5*j)] <- summary(res)$coefficients[2:6, 4]
  }
  flm.result[(2*i-1):(2*i), 51] <- rep(indList[i], 2)
}
colnames(flm.result) <- c(rep(c("Mon", "Tue", "Wed", "Thu", "Fri"), 10), "Industry")

fd.raw <- na.omit(fd.raw)
fd.raw.flm <- fd.raw[, -5]
ret.fdata <- fdata(fd.raw.flm)
fd.raw <- as.data.frame(fd.raw) %>% na.omit()
tt <- ret.fdata[["argvals"]]
nbasis.x <- 4
nbasis.b <- 4
basis1 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.x, norder = 4)
basis2 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.b, norder = 4)
basis.x <- list("x" = basis1)
basis.b <- list("x" = basis2)
ldata <- list("df" = fd.raw, "x" = ret.fdata)

res <- fregre.lm(fri ~ x, ldata, basis.x=basis.x, basis.b=basis.b)
plot(res$coefficients)
res <- fregre.basis(ret.fdata, fd.raw$fri, basis.x=basis1, basis.b=basis2)
summary(res)

m1100 <- daily %>% filter(code %in% "M1100" & weekday %not in% "Sat")
m1100$mon <- ifelse(m1100$weekday %in% "Mon", 1, 0)
m1100$tue <- ifelse(m1100$weekday %in% "Tues", 1, 0)
m1100$wed <- ifelse(m1100$weekday %in% "Wed", 1, 0)
m1100$thu <- ifelse(m1100$weekday %in% "Thurs", 1, 0)
m1100.raw <- m1100 %>%
  dplyr::select(mon, tue, wed, thu)
m1100.raw <- as.matrix(m1100.raw)
m1100.fdata <- fdata(m1100.raw)
m1100.raw <- m1100 %>%
  dplyr::select(ret, mon, tue, wed, thu)

tt <- m1100.fdata[["argvals"]]
nbasis.x <- 4
nbasis.b <- 4
norder.x <- 3
norder.b <- 3
basis1 <- create.bspline.basis(rangeval = range(tt), norder = norder.x, nbasis = nbasis.x)
basis2 <- create.bspline.basis(rangeval = range(tt), norder = norder.b, nbasis = nbasis.b)
basis.x <- list("x" = basis1)
basis.b <- list("x" = basis2)
ldata <- list("df" = m1100.raw, "x" = m1100.fdata)

reg <- lm(ret ~ mon + tue + wed + thu, m1100.raw)
summary(reg)
res <- fregre.lm(ret ~ mon, ldata, basis.x=basis.x, basis.b=basis.b)
res <- fregre.basis(m1100.fdata, m1100.raw$ret, basis.x=basis1, basis.b=basis2)
summary(res)
plot(res)
lm(ret ~ fri, m1100.raw)

####---monthly---####
monthly <- read.table("monthly.txt", header = T, sep = "\t", stringsAsFactors = F, fill = T)
monthly <- monthly[, -3]
colnames(monthly) <- c("code", "name", "yyyymm", "ret")
monthly$code <- gsub("\\s+", "", monthly$code)
monthly$name <- gsub("\\s+", "", monthly$name)
monthly$yyyymm <- as.numeric(monthly$yyyymm)
monthly$ret <- as.numeric(monthly$ret)

monthly %<>%
  mutate(yyyy = paste0(substring(monthly$yyyymm, 1, 4)),
         mm = paste0(substring(monthly$yyyymm, 5, 6))) %>%
  arrange(code, yyyymm)

monthly %>% filter(code %in% "M1100") %>%
  ggplot(., aes(x = mm, y = ret, group = yyyy, color = yyyy)) +
  geom_line(size = 1) +
  labs(title = "M1100",
       x = "Month", y = "Return (%)") +
  theme_bw() 

#fpca
fd.raw %<>% na.omit()
ret.fdata %<>% na.omit()
pc.svd <- fdata2pc(ret.fdata, ncomp = 4) 
summary(pc.svd, fd.raw, biplot=T)
plot(pc.svd$rotation)
