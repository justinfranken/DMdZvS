###############################################################################
### ----- Funktionen --- ######################################################
# ----- Daten ----- ###########################################################
daten = function(start = 1985, end = 2007){
  temp = tempfile()
  download.file("https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2021-10.csv",temp)
  Makrodata = read.csv(temp)
  unlink(temp)
  temp <- tempfile()
  download.file("http://qed.econ.queensu.ca/jae/2014-v29.1/jungbacker-koopman-van_der_wel/jkv-data.zip",temp)
  data <- read.table(unz(temp, "UnsmFB_70-09.txt"))
  unlink(temp)
  jungbacker.data = ts(data, start = c(1970, 1), frequency = 12)
  colnames(jungbacker.data) = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
  CU = Makrodata[c(-1,-755),"CUMFNS"]
  CU.full = ts(CU, start = 1959, end = c(2021,9), frequency = 12)
  FFR = Makrodata[c(-1,-755),"FEDFUNDS"]
  FFR.full = ts(FFR, start = 1959, end = c(2021,9), frequency = 12)
  Infl = Makrodata[c(-1,-755),"CPIAUCSL"]
  INFL = numeric(753)
  for (i in 1:753) {
    INFL[i+12] = (Infl[i+12]/Infl[i]*100)-100
  }
  INFL.full = ts(INFL, start = 1959, end = c(2021,9), frequency = 12)
  CU = window(CU.full, start = start, end = end)
  FFR = window(FFR.full, start = start, end = end)
  INFL = window(INFL.full, start = start, end = end)
  yields = window(jungbacker.data, start = start, end = end)
  output = list("CU.full" = CU.full, "FFR.full" = FFR.full, "INFL.full" = INFL.full, "CU" = CU, "FFR" = FFR, 
                "INFL"= INFL, "jungbacker.data" = jungbacker.data, "yields" = yields)
  return(output)
}
data = daten()
yields = data$yields
CU = data$CU
FFR = data$FFR
INFL = data$INFL
rm(daten, data)
# ----- Modellierung --- ######################################################
faktoranalyse = function(x, r = 3){ # x: Zinsdaten , r : Anzahl an Faktoren
  yield.means.jb = colMeans(x)
  centered_jungbacker = scale(x, scale = FALSE)
  yield.cov = cov(centered_jungbacker)
  yields.eigen = eigen(yield.cov)
  yields.eigenvec = yields.eigen$vectors[,1:r]
  yields.eigenval = yields.eigen$values[1:r]
  factors = matrix(0,nrow = dim(centered_jungbacker)[1], ncol = r)
  for (i in 1:dim(centered_jungbacker)[1]){
    factors[i,] = lm(centered_jungbacker[i,] ~ yields.eigenvec[,1:r] - 1)$coefficients
  }
  yields.f = matrix(nrow = dim(x)[1], ncol = 17)
  for (i in 1:dim(x)[1]) {
    yields.f[i,] = (factors[i,] %*% t(yields.eigenvec[,1:r])) + yield.means.jb
  }
  output = list("zinsen.f" = yields.f, "faktoren" = factors, "eigenvec" = yields.eigenvec, 
                "mu" = yield.means.jb, "eigenval" = yields.eigenval)
  return(output)
}
diebold.li = function(x, lambda = 0.0609){
  NS <- matrix(nrow=17, ncol=3)
  laufzeit = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
  for (i in 1:17){
    NS[i, 1] = 1
    NS[i, 2] = (1-exp(-lambda*laufzeit[i]))/(lambda * laufzeit[i])
    NS[i, 3] = (1-exp(-lambda*laufzeit[i]))/(lambda * laufzeit[i]) - exp(- lambda * laufzeit[i])}
  beta.dl <- matrix(nrow = dim(x)[1], ncol = 3)
  for (i in 1:dim(x)[1]){
    beta.dl[i,] = lm(yields[i,] ~ NS - 1)$coefficients}
  fine.grid = 3:120
  NS2 <- matrix(nrow=118, ncol=3)
  for (i in 1:118){
    NS2[i, 1] = 1
    NS2[i, 2] = (1-exp(-lambda*fine.grid[i]))/(lambda * fine.grid[i])
    NS2[i, 3] = (1-exp(-lambda*fine.grid[i]))/(lambda * fine.grid[i]) - exp(- lambda * fine.grid[i])}
  yields.dl = beta.dl %*% t(NS2)
  output = list("zinsen.dl" = yields.dl, "beta.dl" = beta.dl, "NS.dl" = NS, "fine.grid" = NS2)
  return(output)
}
svensson = function(x, tau1 = 0.0609, tau2 = 0.25){
  Sline = matrix(nrow = 17, ncol = 4)
  laufzeit = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
  for (i in 1:17){
    Sline[i, 1] <- 1
    Sline[i, 2] <- (1-exp(-tau1*laufzeit[i]))/(tau1 * laufzeit[i])
    Sline[i, 3] <- (1-exp(-tau1*laufzeit[i]))/(tau1 * laufzeit[i]) - exp(- tau1 * laufzeit[i])
    Sline[i, 4] <- (1-exp(-tau2*laufzeit[i]))/(tau2 * laufzeit[i]) - exp(- tau2 * laufzeit[i])
  }
  beta.s = matrix(nrow=dim(x)[1], ncol = 4)
  for (i in 1:dim(x)[1]){
    beta.s[i,] = lm(yields[i,] ~ Sline - 1)$coefficients
  }
  Sline2 = matrix(nrow = 118, ncol = 4)
  fine.grid = 3:120
  for (i in 1:118){
    Sline2[i, 1] <- 1
    Sline2[i, 2] <- (1-exp(-tau1*fine.grid[i]))/(tau1 * fine.grid[i])
    Sline2[i, 3] <- (1-exp(-tau1*fine.grid[i]))/(tau1 * fine.grid[i]) - exp(- tau1 * fine.grid[i])
    Sline2[i, 4] <- (1-exp(-tau2*fine.grid[i]))/(tau2 * fine.grid[i]) - exp(- tau2 * fine.grid[i])
  }
  yields.s = beta.s %*% t(Sline2)
  output = list("zinsen.s" = yields.s, "beta.s" = beta.s, "NS.s" = Sline, "fine.grid" = Sline2)
  return(output)
}
# ----- Vorhersage --- ########################################################
prognose.faktoranalyse = function(x, r = 3, h = 1, VAR = FALSE, makro = FALSE, makro.var = cbind(CU,FFR,INFL)){
  if(r == 1)stop("r muss größer 1 sein!")
  factors = faktoranalyse(x,r)$faktoren
  mu = faktoranalyse(x,r)$mu
  eigenvec = faktoranalyse(x,r)$eigenvec
  if(makro == FALSE){
    mean.factors = colMeans(factors)
    for (i in 1:dim(factors)[1]){
      factors[i,] = factors[i,] - mean.factors
    }
    if(VAR == FALSE){
      theta.f = numeric(r)
      for (i in 1:r) {
        theta.f[i] = lm(embed(factors,2)[,i] ~ -1 + embed(factors, 2)[,i+r])$coefficients
      }
      pred.fac = theta.f * factors[dim(factors)[1],]
      if(h > 1){
        pred.factor = matrix(nrow = h, ncol = r)
        pred.factor[1,] = pred.fac
        for (i in 1:(h-1)){
          pred.factor[i+1,] = theta.f * pred.factor[i,]
        }
        for (i in 1:h) {
          pred.factor[i,] = pred.factor[i,] + mean.factors
        }
      }
      else{pred.factor = pred.fac + mean.factors}
      pred.yield.f = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h){
          pred.yield.f[i,] = mu + pred.factor[i,] %*% t(eigenvec)
        }
      }
      else{pred.yield.f = mu + pred.factor %*% t(eigenvec)}
    }
    else{
      Theta.f = lm(embed(factors,2)[,1:r] ~ embed(factors,2)[,-(1:r)] - 1)$coefficients
      dimnames(Theta.f) = list(1:dim(Theta.f)[1], 1:dim(Theta.f)[2])
      pred.fac = Theta.f %*% factors[dim(factors)[1],]
      if(h > 1){
        pred.factor = matrix(nrow = h, ncol = r)
        pred.factor[1,] = pred.fac
        for (i in 1:(h-1)) {
          pred.factor[i+1,] = Theta.f %*% pred.factor[i,]
        }
        for (i in 1:h) {
         pred.factor[i,] = pred.factor[i,] + mean.factors 
        }
      }
      else{pred.factor = pred.fac + mean.factors}
      pred.yield.f = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h){
          pred.yield.f[i,] = mu + pred.factor[i,] %*% t(eigenvec)
        }
      }
      else{pred.yield.f = mu + t(pred.factor) %*% t(eigenvec)}
    }
  }
  if(makro == TRUE){
    factors.makro = cbind(factors, makro.var)
    mean.factors.makro = colMeans(factors.makro)
    for (i in 1:dim(factors.makro)[1]){
      factors.makro[i,] = factors.makro[i,] - mean.factors.makro
    }
    if(VAR == FALSE){
      theta.f.makro = matrix(nrow = r+1, ncol = r)
      for (i in 1:r){
        theta.f.makro[,i] = lm(embed(factors.makro,2)[,i] ~ -1 + embed(factors.makro,2)[,c(i+r+3, dim(embed(factors.makro,2))[2]-2,dim(embed(factors.makro,2))[2]-1,dim(embed(factors.makro,2))[2])])$coefficients
      }
      pred.fac = theta.f.makro[1,] * factors.makro[dim(factors.makro)[1],1:r]
      if(h > 1){
        pred.factor = matrix(nrow = h, ncol = r)
        pred.factor[1,] = pred.fac
        for (i in 1:(h-1)){
          pred.factor[i+1,] = theta.f.makro[1,] * pred.factor[i,]
        }
        for (i in 1:h) {
          pred.factor[i,] = pred.factor[i,] + mean.factors.makro[1:r]
        }
      }
      else{pred.factor = pred.fac + mean.factors.makro[1:r]}
      pred.yield.f = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h){
          pred.yield.f[i,] = mu + pred.factor[i,] %*% t(eigenvec)
        }
      }
      else{pred.yield.f = mu + pred.factor %*% t(eigenvec)}
    }
    else{
      Theta.f.makro = lm(embed(factors.makro,2)[,1:(r+3)] ~ -1 + embed(factors.makro,2)[,(dim(embed(factors.makro,2))[2]-r-2):dim(embed(factors.makro,2))[2]])$coefficients
      dimnames(Theta.f.makro) = list(1:dim(Theta.f.makro)[1], 1:dim(Theta.f.makro)[2])
      pred.fac = Theta.f.makro[1:r,1:r] %*% factors.makro[dim(factors.makro)[1],1:r]
      if(h > 1){
        pred.factor = matrix(nrow = h, ncol = r)
        pred.factor[1,] = pred.fac
        for (i in 1:(h-1)) {
          pred.factor[i+1,] = Theta.f.makro[1:r,1:r] %*% pred.factor[i,]
        }
        for (i in 1:h) {
          pred.factor[i,] = pred.factor[i,] + mean.factors.makro[1:r]
        }
      }
      else{pred.factor = pred.fac + mean.factors.makro[1:r]}
      pred.yield.f = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h){
          pred.yield.f[i,] = mu + pred.factor[i,] %*% t(eigenvec)
        }
      }
      else{pred.yield.f = mu + t(pred.factor) %*% t(eigenvec)}
    }
  }
  if(VAR == FALSE & makro == FALSE){theta.f = theta.f}
  if(VAR == TRUE & makro == FALSE){theta.f = Theta.f}
  if(VAR == FALSE & makro == TRUE){theta.f = theta.f.makro}
  if(VAR == TRUE & makro == TRUE){theta.f = Theta.f.makro}
  output = list("prog.zins" = pred.yield.f, "prog.faktor" = pred.factor, "theta.f" = theta.f)
  return(output)
}
prognose.dieboldli = function(x, h = 1, lambda = 0.0609, VAR = FALSE, makro = FALSE, makro.var = cbind(CU,FFR,INFL)){
  beta.dl = diebold.li(x, lambda = lambda)$beta.dl
  NS = diebold.li(x, lambda = lambda)$NS.dl
  if(makro == FALSE){
    mean.beta.dl = colMeans(beta.dl)
    for (i in 1:dim(beta.dl)[1]){
      beta.dl[i,] = beta.dl[i,] - mean.beta.dl
    }
    if(VAR == FALSE){
      theta.dl = numeric(3)
      for (i in 1:3) {
        theta.dl[i] = lm(embed(beta.dl,2)[,i] ~ embed(beta.dl, 2)[,i+3] - 1)$coefficients
      }     
      pred.dl = theta.dl * beta.dl[dim(beta.dl)[1],]
      if(h > 1){
        pred.beta.dl = matrix(nrow = h, ncol = 3)
        pred.beta.dl[1,] = pred.dl
        for (i in 1:(h-1)) {
          pred.beta.dl[i+1,] = theta.dl * pred.beta.dl[i,]
        }
        for (i in 1:h) {
          pred.beta.dl[i,] = pred.beta.dl[i,] + mean.beta.dl
        }
      }
      else{pred.beta.dl = pred.dl + mean.beta.dl}
      if(h > 1){
        pred.yield.dl = matrix(nrow = h, ncol = 17)
        for (i in 1:h) {
         pred.yield.dl[i,] = pred.beta.dl[i,] %*% t(NS) 
        } 
      }
      else{pred.yield.dl = pred.beta.dl %*% t(NS)}
    }
    else{
      Theta.dl = lm(embed(beta.dl,2)[,1:3] ~ embed(beta.dl,2)[,-(1:3)] - 1)$coefficients
      dimnames(Theta.dl) = list(1:dim(Theta.dl)[1], 1:dim(Theta.dl)[2])
      pred.dl = Theta.dl %*% beta.dl[dim(beta.dl)[1],]
      if(h > 1){
        pred.beta.dl = matrix(nrow = h, ncol = 3)
        pred.beta.dl[1,] = pred.dl
        for (i in 1:(h-1)) {
          pred.beta.dl[i+1,] = Theta.dl %*% pred.beta.dl[i,]
        }
        for (i in 1:h) {
          pred.beta.dl[i,] = pred.beta.dl[i,] + mean.beta.dl
        }
      }
      else{pred.beta.dl = pred.dl + mean.beta.dl}
      pred.yield.dl = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h) {
          pred.yield.dl[i,] = pred.beta.dl[i,] %*% t(NS)
        }
      }
      else{pred.yield.dl = t(pred.beta.dl) %*% t(NS)}
    }
  }
  if(makro == TRUE){
    beta.dl.makro = cbind(beta.dl, makro.var)
    mean.beta.dl.makro = colMeans(beta.dl.makro)
    for (i in 1:dim(beta.dl.makro)[1]) {
      beta.dl.makro[i,] = beta.dl.makro[i,] - mean.beta.dl.makro
    }
    if(VAR == FALSE){
      theta.dl.makro = matrix(nrow = 4, ncol = 3)
      for (i in 1:3) {
        theta.dl.makro[,i] = lm(embed(beta.dl.makro,2)[,i] ~ -1 + embed(beta.dl.makro,2)[,c((i+6),10,11,12)])$coefficients
      }
      pred.dl = theta.dl.makro[1,] * beta.dl.makro[dim(beta.dl.makro)[1],1:3]
      if(h > 1){
        pred.beta.dl = matrix(nrow = h, ncol = 3)
        pred.beta.dl[1,] = pred.dl
        for (i in 1:(h-1)) {
          pred.beta.dl[i+1,] = theta.dl.makro[1,] * pred.beta.dl[i,]
        }
        for (i in 1:h) {
          pred.beta.dl[i,] = pred.beta.dl[i,] + mean.beta.dl.makro[1:3]
        }
      }
      else{pred.beta.dl = pred.dl + mean.beta.dl.makro[1:3]}
      pred.yield.dl = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h) {
          pred.yield.dl[i,] = pred.beta.dl[i,] %*% t(NS) 
        }
      }
      else{pred.yield.dl = pred.beta.dl %*% t(NS)}
    }
    else{
      Theta.dl.makro = lm(embed(beta.dl.makro,2)[,1:6] ~ -1 + embed(beta.dl.makro,2)[,7:12])$coefficients
      dimnames(Theta.dl.makro) = list(1:dim(Theta.dl.makro)[1], 1:dim(Theta.dl.makro)[2])
      pred.dl = Theta.dl.makro[1:3,1:3] %*% beta.dl.makro[dim(beta.dl.makro)[1],1:3]
      if(h > 1){
        pred.beta.dl = matrix(nrow = h, ncol = 3)
        pred.beta.dl[1,] = pred.dl
        for (i in 1:(h-1)) {
          pred.beta.dl[i+1,] = Theta.dl.makro[1:3,1:3] %*% pred.beta.dl[i,]
        }
        for (i in 1:h) {
          pred.beta.dl[i,] = pred.beta.dl[i,] + mean.beta.dl.makro[1:3]
        }
      }
      else{pred.beta.dl = pred.dl + mean.beta.dl.makro[1:3]}
      pred.yield.dl = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h) {
          pred.yield.dl[i,] = pred.beta.dl[i,] %*% t(NS)
        }
      }
      else{pred.yield.dl = t(pred.beta.dl) %*% t(NS)}
    }
  }
  if(VAR == FALSE & makro == FALSE){theta.dl = theta.dl}
  if(VAR == TRUE & makro == FALSE){theta.dl = Theta.dl}
  if(VAR == FALSE & makro == TRUE){theta.dl = theta.dl.makro}
  if(VAR == TRUE & makro == TRUE){theta.dl = Theta.dl.makro}
  output = list("prog.zins" = pred.yield.dl, "prog.faktor" = pred.beta.dl, "theta.dl" = theta.dl)
  return(output)
}
prognose.svensson = function(x, h = 1, tau1 = 0.0609, tau2 = 0.25, VAR = FALSE, makro = FALSE, makro.var = cbind(CU,FFR,INFL)){
  beta.s = svensson(x, tau1 = tau1, tau2 = tau2)$beta.s
  NS = svensson(x, tau1 = tau1, tau2 = tau2)$NS.s
  if(makro == FALSE){
    mean.beta.s = colMeans(beta.s)
    for (i in 1:dim(beta.s)[1]){
      beta.s[i,] = beta.s[i,] - mean.beta.s
    }
    if(VAR == FALSE){
      theta.s = numeric(4)
      for (i in 1:4) {
        theta.s[i] = lm(embed(beta.s,2)[,i] ~ embed(beta.s, 2)[,i+4] - 1)$coefficients
      }     
      pred.s = theta.s * beta.s[dim(beta.s)[1],]
      if(h > 1){
        pred.beta.s = matrix(nrow = h, ncol = 4)
        pred.beta.s[1,] = pred.s
        for (i in 1:(h-1)) {
          pred.beta.s[i+1,] = theta.s * pred.beta.s[i,]
        }
        for (i in 1:h) {
          pred.beta.s[i,] = pred.beta.s[i,] + mean.beta.s
        }
      }
      else{pred.beta.s = pred.s + mean.beta.s}
      if(h > 1){
        pred.yield.s = matrix(nrow = h, ncol = 17)
        for (i in 1:h) {
          pred.yield.s[i,] = pred.beta.s[i,] %*% t(NS) 
        } 
      }
      else{pred.yield.s = pred.beta.s %*% t(NS)}
    }
    else{
      Theta.s = lm(embed(beta.s,2)[,1:4] ~ embed(beta.s,2)[,-(1:4)] - 1)$coefficients
      dimnames(Theta.s) = list(1:dim(Theta.s)[1], 1:dim(Theta.s)[2])
      pred.s = Theta.s %*% beta.s[dim(beta.s)[1],]
      if(h > 1){
        pred.beta.s = matrix(nrow = h, ncol = 4)
        pred.beta.s[1,] = pred.s
        for (i in 1:(h-1)) {
          pred.beta.s[i+1,] = Theta.s %*% pred.beta.s[i,]
        }
        for (i in 1:h) {
          pred.beta.s[i,] = pred.beta.s[i,] + mean.beta.s
        }
      }
      else{pred.beta.s = pred.s + mean.beta.s}
      pred.yield.s = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h) {
          pred.yield.s[i,] = pred.beta.s[i,] %*% t(NS)
        }
      }
      else{pred.yield.s = t(pred.beta.s) %*% t(NS)}
    }
  }
  if(makro == TRUE){
    beta.s.makro = cbind(beta.s, makro.var)
    mean.beta.s.makro = colMeans(beta.s.makro)
    for (i in 1:dim(beta.s.makro)[1]) {
      beta.s.makro[i,] = beta.s.makro[i,] - mean.beta.s.makro
    }
    if(VAR == FALSE){
      theta.s.makro = matrix(nrow = 4, ncol = 4)
      for (i in 1:4) {
        theta.s.makro[,i] = lm(embed(beta.s.makro,2)[,i] ~ -1 + embed(beta.s.makro,2)[,c((i+7),12,13,14)])$coefficients
      }
      pred.s = theta.s.makro[1,] * beta.s.makro[dim(beta.s.makro)[1],1:4]
      if(h > 1){
        pred.beta.s = matrix(nrow = h, ncol = 4)
        pred.beta.s[1,] = pred.s
        for (i in 1:(h-1)) {
          pred.beta.s[i+1,] = theta.s.makro[1,] * pred.beta.s[i,]
        }
        for (i in 1:h) {
          pred.beta.s[i,] = pred.beta.s[i,] + mean.beta.s.makro[1:4]
        }
      }
      else{pred.beta.s = pred.s + mean.beta.s.makro[1:4]}
      pred.yield.s = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h) {
          pred.yield.s[i,] = pred.beta.s[i,] %*% t(NS) 
        }
      }
      else{pred.yield.s = pred.beta.s %*% t(NS)}
    }
    else{
      Theta.s.makro = lm(embed(beta.s.makro,2)[,1:7] ~ embed(beta.s.makro,2)[,8:14] - 1)$coefficients
      dimnames(Theta.s.makro) = list(1:dim(Theta.s.makro)[1], 1:dim(Theta.s.makro)[2])
      pred.s = Theta.s.makro[1:4,1:4] %*% beta.s.makro[dim(beta.s.makro)[1],1:4]
      if(h > 1){
        pred.beta.s = matrix(nrow = h, ncol = 4)
        pred.beta.s[1,] = pred.s
        for (i in 1:(h-1)) {
          pred.beta.s[i+1,] = Theta.s.makro[1:4,1:4] %*% pred.beta.s[i,]
        }
        for (i in 1:h) {
          pred.beta.s[i,] = pred.beta.s[i,] + mean.beta.s.makro[1:4]
        }
      }
      else{pred.beta.s = pred.s + mean.beta.s.makro[1:4]}
      pred.yield.s = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h) {
          pred.yield.s[i,] = pred.beta.s[i,] %*% t(NS)
        }
      }
      else{pred.yield.s = t(pred.beta.s) %*% t(NS)}
    }
  }
  if(VAR == FALSE & makro == FALSE){theta.s = theta.s}
  if(VAR == TRUE & makro == FALSE){theta.s = Theta.s}
  if(VAR == FALSE & makro == TRUE){theta.s = theta.s.makro}
  if(VAR == TRUE & makro == TRUE){theta.s = Theta.s.makro}
  output = list("prog.zins" = pred.yield.s, "prog.faktor" = pred.beta.s, "theta.s" = theta.s)
  return(output)
}
# ----- RMSE --- ##############################################################
# Faktoranalyse
rmse.1.f = function(x, r, VAR = FALSE, makro = FALSE, start.training = 100){ 
  wiederholungen = dim(x)[1]-start.training
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    training.makro.var = cbind(CU[1:(start.training+i-1)],FFR[1:(start.training+i-1)],INFL[1:(start.training+i-1)])
    y_hat = prognose.faktoranalyse(trainings.data, VAR = VAR, makro = makro, makro.var = training.makro.var)$prog.zins
    y = yields[start.training+i,]
    pred.error[i] = mean(y-y_hat)^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
rmse.6.f = function(x, r, VAR = FALSE, makro = FALSE, start.training = 100){ 
  wiederholungen = dim(x)[1]-start.training-6
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    training.makro.var = cbind(CU[1:(start.training+i-1)],FFR[1:(start.training+i-1)],INFL[1:(start.training+i-1)])
    y_hat = prognose.faktoranalyse(trainings.data, VAR = VAR, makro = makro, h = 6, makro.var = training.makro.var)$prog.zins
    y = yields[start.training+i+5,]
    pred.error[i] = mean(y-y_hat[6,])^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
rmse.12.f = function(x, r, VAR = FALSE, makro = FALSE, start.training = 100){ 
  wiederholungen = dim(x)[1]-start.training-12
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    training.makro.var = cbind(CU[1:(start.training+i-1)],FFR[1:(start.training+i-1)],INFL[1:(start.training+i-1)])
    y_hat = prognose.faktoranalyse(trainings.data, VAR = VAR, makro = makro, h = 12, makro.var = training.makro.var)$prog.zins
    y = yields[start.training+i+11,]
    pred.error[i] = mean(y-y_hat[12,])^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
# Diebold und Li
rmse.1.dl = function(x, lambda = 0.0609, VAR = FALSE, makro = FALSE, start.training = 100){
  wiederholungen = dim(x)[1]-start.training
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    training.makro.var = cbind(CU[1:(start.training+i-1)],FFR[1:(start.training+i-1)],INFL[1:(start.training+i-1)])
    y_hat = prognose.dieboldli(trainings.data, lambda = lambda, VAR = VAR, makro = makro, makro.var = training.makro.var)$prog.zins
    y = yields[start.training+i,]
    pred.error[i] = mean(y-y_hat)^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
rmse.6.dl = function(x, lambda = 0.0609, VAR = FALSE, makro = FALSE, start.training = 100){ 
  wiederholungen = dim(x)[1]-start.training-6
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    training.makro.var = cbind(CU[1:(start.training+i-1)],FFR[1:(start.training+i-1)],INFL[1:(start.training+i-1)])
    y_hat = prognose.dieboldli(trainings.data, lambda = lambda, VAR = VAR, makro = makro, h = 6, makro.var = training.makro.var)$prog.zins
    y = yields[start.training+i+5,]
    pred.error[i] = mean(y-y_hat[6,])^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
rmse.12.dl = function(x, lambda = 0.0609, VAR = FALSE, makro = FALSE, start.training = 100){ 
  wiederholungen = dim(x)[1]-start.training-12
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    training.makro.var = cbind(CU[1:(start.training+i-1)],FFR[1:(start.training+i-1)],INFL[1:(start.training+i-1)])
    y_hat = prognose.dieboldli(trainings.data, lambda = lambda, VAR = VAR, makro = makro, h = 12, makro.var = training.makro.var)$prog.zins
    y = yields[start.training+i+11,]
    pred.error[i] = mean(y-y_hat[12,])^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
# svensson
rmse.1.s = function(x, tau1 = 0.0609, tau2 = 0.25, VAR = FALSE, makro = FALSE, start.training = 100){
  wiederholungen = dim(x)[1]-start.training
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    training.makro.var = cbind(CU[1:(start.training+i-1)],FFR[1:(start.training+i-1)],INFL[1:(start.training+i-1)])
    y_hat = prognose.svensson(trainings.data, tau1 = tau1, tau2 = tau2, VAR = VAR, makro = makro, makro.var = training.makro.var)$prog.zins
    y = yields[start.training+i,]
    pred.error[i] = mean(y-y_hat)^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
rmse.6.s = function(x, tau1 = 0.0609, tau2 = 0.25, VAR = FALSE, makro = FALSE, start.training = 100){ 
  wiederholungen = dim(x)[1]-start.training-6
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    training.makro.var = cbind(CU[1:(start.training+i-1)],FFR[1:(start.training+i-1)],INFL[1:(start.training+i-1)])
    y_hat = prognose.svensson(trainings.data, tau1 = tau1, tau2 = tau2, VAR = VAR, makro = makro, h = 6, makro.var = training.makro.var)$prog.zins
    y = yields[start.training+i+5,]
    pred.error[i] = mean(y-y_hat[6,])^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
rmse.12.s = function(x, tau1 = 0.0609, tau2 = 0.25, VAR = FALSE, makro = FALSE, start.training = 100){ 
  wiederholungen = dim(x)[1]-start.training-12
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    training.makro.var = cbind(CU[1:(start.training+i-1)],FFR[1:(start.training+i-1)],INFL[1:(start.training+i-1)])
    y_hat = prognose.svensson(trainings.data, tau1 = tau1, tau2 = tau2, VAR = VAR, makro = makro, h = 12, makro.var = training.makro.var)$prog.zins
    y = yields[start.training+i+11,]
    pred.error[i] = mean(y-y_hat[12,])^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
### ----- Tabellen --- ########################################################
### Tabelle 1/3: Vorhersage der Modelle ####
# f
# VAR = FALSE, makro = FALSE
rmse.1.f. = rmse.1.f(yields, VAR = FALSE, makro = FALSE) 
# rmse: 0.241589, mse: 0.05836523, max.mse: 0.4861628, min.mse: 5.619721e-06, sd.mse: 0.08455746
rmse.6.f. = rmse.6.f(yields, VAR = FALSE, makro = FALSE) 
#rmse: 0.7562787, mse: 0.5719575, max.mse: 4.148918, min.mse: 4.674727e-08, sd.mse: 0.7489498
rmse.12.f. = rmse.12.f(yields, VAR = FALSE, makro = FALSE) 
#rmse: 1.188187, mse: 1.411788, max.mse: 10.41933, min.mse: 0.0007847895, sd.mse: 1.824749
# VAR = TRUE, makro = FALSE
rmse.1.f.var = rmse.1.f(yields, VAR = TRUE, makro = FALSE) 
#rmse: 0.2481085, mse: 0.06155781, max.mse: 0.5336497, min.mse: 7.113894e-06, sd.mse: 0.0883579
rmse.6.f.var = rmse.6.f(yields, VAR = TRUE, makro = FALSE) 
#rmse: 0.788221, mse: 0.6212923, max.mse: 4.974624, min.mse: 6.431074e-06, sd.mse: 0.8477025
rmse.12.f.var = rmse.12.f(yields, VAR = TRUE, makro = FALSE) 
#rmse: 1.223953, mse: 1.49806, max.mse: 10.64188, min.mse: 0.0002900088, sd.mse: 2.017627
# VAR = FALSE, makro = TRUE
rmse.1.f.makro = rmse.1.f(yields, VAR = FALSE, makro = TRUE) 
#rmse: 0.2380791, mse: 0.05668167, max.mse: 0.4297689, min.mse: 3.106965e-07, sd.mse: 0.08022008
rmse.6.f.makro = rmse.6.f(yields, VAR = FALSE, makro = TRUE) 
#rmse: 0.7035496, mse: 0.494982, max.mse: 3.061104, min.mse: 9.227663e-06, sd.mse: 0.6102703
rmse.12.f.makro = rmse.12.f(yields, VAR = FALSE, makro = TRUE) 
#rmse: 1.04821, mse: 1.098745, max.mse: 10.19903, min.mse: 8.17595e-05, sd.mse: 1.544282
# VAR = TRUE, makro = TRUE
rmse.1.f.makro.var = rmse.1.f(yields, VAR = TRUE, makro = TRUE) 
#rmse: 0.2940126, mse: 0.08644342, max.mse: 0.7076973, min.mse: 5.975465e-06, sd.mse: 0.1248889
rmse.6.f.makro.var = rmse.6.f(yields, VAR = TRUE, makro = TRUE) 
#rmse: 1.019944, mse: 1.040285, max.mse: 6.578723, min.mse: 4.620514e-06, sd.mse: 1.27031
rmse.12.f.makro.var = rmse.12.f(yields, VAR = TRUE, makro = TRUE) 
#rmse: 1.537018, mse: 2.362424, max.mse: 11.25366, min.mse: 0.0006439672, sd.mse: 2.825939
# dl
# VAR = FALSE, makro = FALSE
rmse.1.dl. = rmse.1.dl(yields, VAR = FALSE, makro = FALSE) 
#rmse: 0.2537664, mse: 0.06439736, max.mse: 0.5513198, min.mse: 1.490709e-05, sd.mse: 0.09224411
rmse.6.dl. = rmse.6.dl(yields, VAR = FALSE, makro = FALSE) 
#rmse: 0.8250929, mse: 0.6807782, max.mse: 4.361179, min.mse: 6.666446e-05, sd.mse: 0.8511514
rmse.12.dl. = rmse.12.dl(yields, VAR = FALSE, makro = FALSE) 
#rmse: 1.264988, mse: 1.600194, max.mse: 12.28705, min.mse: 9.240488e-07, sd.mse: 2.091235
# VAR = TRUE, makro = FALSE
rmse.1.dl.var = rmse.1.dl(yields, VAR = TRUE, makro = FALSE) 
#rmse: 0.2535021, mse: 0.06426329, max.mse: 0.5663474, min.mse: 4.216018e-07, sd.mse: 0.09702723
rmse.6.dl.var = rmse.6.dl(yields, VAR = TRUE, makro = FALSE) 
#rmse: 0.7997477, mse: 0.6395965, max.mse: 4.833398, min.mse: 3.696638e-06, sd.mse: 0.8213084
rmse.12.dl.var = rmse.12.dl(yields, VAR = TRUE, makro = FALSE) 
#rmse: 1.303503, mse: 1.699121, max.mse: 10.73264, min.mse: 0.000189318, sd.mse: 2.007788
# VAR = FALSE, makro = TRUE
rmse.1.dl.makro = rmse.1.dl(yields, VAR = FALSE, makro = TRUE) 
#rmse: 0.2679634, mse: 0.07180439, max.mse: 0.5949264, min.mse: 3.458666e-06, sd.mse: 0.1045168
rmse.6.dl.makro = rmse.6.dl(yields, VAR = FALSE, makro = TRUE) 
#rmse: 0.898805, mse: 0.8078505, max.mse: 5.473214, min.mse: 1.611047e-05, sd.mse: 1.049348
rmse.12.dl.makro = rmse.12.dl(yields, VAR = FALSE, makro = TRUE) 
#rmse: 1.341472, mse: 1.799548, max.mse: 10.12184, min.mse: 6.932526e-05, sd.mse: 2.313549
# VAR = TRUE, makro = TRUE
rmse.1.dl.makro.var = rmse.1.dl(yields, VAR = TRUE, makro = TRUE) 
#rmse: 0.9652989, mse: 0.9318019, max.mse: 11.57259, min.mse: 2.424098e-05, sd.mse: 1.349574
rmse.6.dl.makro.var = rmse.6.dl(yields, VAR = TRUE, makro = TRUE) 
#rmse: 3.225115, mse: 10.40136, max.mse: 156.0724, min.mse: 0.001813201, sd.mse: 16.45599
rmse.12.dl.makro.var = rmse.12.dl(yields, VAR = TRUE, makro = TRUE) 
#rmse: 3.927738, mse: 15.42713, max.mse: 260.5435, min.mse: 3.203196e-05, sd.mse: 29.69846
# s
# VAR = FALSE, makro = FALSE
rmse.1.s. = rmse.1.s(yields, VAR = FALSE, makro = FALSE) 
#rmse: 0.2593807, mse: 0.06727833, max.mse: 0.6134914, min.mse: 7.794579e-06, sd.mse: 0.09792129
rmse.6.s. = rmse.6.s(yields, VAR = FALSE, makro = FALSE) 
#rmse: 0.8506164, mse: 0.7235483, max.mse: 5.084359, min.mse: 0.0003254435, sd.mse: 0.946452
rmse.12.s. = rmse.12.s(yields, VAR = FALSE, makro = FALSE) 
#rmse: 1.27991, mse: 1.638169, max.mse: 10.16419, min.mse: 3.168329e-05, sd.mse: 2.170542
# VAR = TRUE, makro = FALSE
rmse.1.s.var = rmse.1.s(yields, VAR = TRUE, makro = FALSE) 
#rmse: 0.3169027, mse: 0.1004273, max.mse: 0.7909522, min.mse: 2.111915e-05, sd.mse: 0.141352
rmse.6.s.var = rmse.6.s(yields, VAR = TRUE, makro = FALSE) 
#rmse: 1.054405, mse: 1.111771, max.mse: 6.124955, min.mse: 1.134298e-07, sd.mse: 1.291428
rmse.12.s.var = rmse.12.s(yields, VAR = TRUE, makro = FALSE) 
#rmse: 1.49889, mse: 2.246672, max.mse: 11.86595, min.mse: 0.000149315, sd.mse: 2.547265
# VAR = FALSE, makro = TRUE
rmse.1.s.makro = rmse.1.s(yields, VAR = FALSE, makro = TRUE) 
#rmse: 0.2763983, mse: 0.07639601, max.mse: 0.6742248, min.mse: 2.845866e-05, sd.mse: 0.1134993
rmse.6.s.makro = rmse.6.s(yields, VAR = FALSE, makro = TRUE) 
#rmse: 0.9236925, mse: 0.8532079, max.mse: 6.033014, min.mse: 9.947908e-07, sd.mse: 1.149848
rmse.12.s.makro = rmse.12.s(yields, VAR = FALSE, makro = TRUE) 
#rmse: 1.352935, mse: 1.830432 max.mse: 9.594866, min.mse: 1.051637e-05, sd.mse: 2.384919
# VAR = TRUE, makro = TRUE
rmse.1.s.makro.var = rmse.1.s(yields, VAR = TRUE, makro = TRUE) 
#rmse: 1.277118, mse: 1.631031, max.mse: 12.65537, min.mse: 0.001452254, sd.mse: 1.771475
rmse.6.s.makro.var = rmse.6.s(yields, VAR = TRUE, makro = TRUE) 
#rmse: 3.525747, mse: 12.43089, max.mse: 258.3894, min.mse: 0.003813188, sd.mse: 22.96441
rmse.12.s.makro.var = rmse.12.s(yields, VAR = TRUE, makro = TRUE) 
#rmse: 4.144647, mse: 17.1781, max.mse: 560.4187, min.mse: 6.569575e-05, sd.mse: 48.68193
### Tabelle 2: Statistiken der Zins und Makro Daten ####
yield3 = yields[,1]
summary(yield3)
round(sd(yield3),3)
yield6 = yields[,2]
yield9 = yields[,3]
yield12 = yields[,4]
yield15 = yields[,5]
yield18 = yields[,6]
yield21 = yields[,7]
yield24 = yields[,8]
yield30 = yields[,9]
yield36 = yields[,10]
yield48 = yields[,11]
yield60 = yields[,12]
yield72 = yields[,13]
yield84 = yields[,14]
yield96 = yields[,15]
yield108 = yields[,16]
yield120 = yields[,17]
summary(yield108)
### ----- Abbildungen --- #####################################################
### Abbildung 3.2: Nelson-Siegel Faktorladungen ####
lambda = 0.0609
NSline = matrix(nrow=120, ncol=3)
for (i in 1:120){ # 120 ist größte Beobachtete Laufzeit
  NSline[i, 1] <- 1
  NSline[i, 2] <- (1-exp(-lambda*i))/(lambda * i)
  NSline[i, 3] <- (1-exp(-lambda*i))/(lambda * i) - exp(- lambda * i)
}
pdf('NelsonSiegel.pdf', width=12, height=6, pointsize = 10)
plot(1:120, NSline[,1], cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.8,
     type = "l", ylim = c(0,1.05), lty=1, lwd = 2.3, xlab = "Laufzeit in Monaten (m)", ylab = "Faktorladung",
     main = "Abbildung 3.2: Nelson-Siegel Faktorladungen")
lines(1:120, NSline[,2], lty = 5, lwd = 2.3)
lines(1:120, NSline[,3], lty = 3, lwd = 2.3)
legend(x = "right", legend = c('Niveau', 'Steigung', 'Krümmung'), lty = c(1,5,3), 
       lwd = c(1.5,1.5,1.5), cex = 1.5)
dev.off()
### Abbildung 5.1: US-amerikanische Staatsanleihen von 1970-2009  ####
library(plot3Drgl)
Abbildung5.1 = function(data, observationgrid = NULL){
  if(is.null(observationgrid)) observationgrid = as.numeric(colnames(data))
  if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers")
  if(length(observationgrid) != ncol(data) ) stop('please specify a valid observationgrid')
  TT <- dim(data)[1]
  time <- replicate(length(observationgrid), c(time(data)))
  surf3D(x = time,
         y = t(replicate(TT, observationgrid)),
         z = data,
         bty="b2",
         theta=-35,
         phi=11,
         colvar = data,
         col = NULL,
         shade = 0.25,
         colkey = FALSE,
         border = "black",
         zlim = range(data) * c(0.9,1.1),
         xlab = 'Zeit',
         ylab ='Laufzeit in Monaten',
         zlab ='Zins (Prozent)',
         ticktype = "detailed",
         main = "Abbildung 5.1: US-amerikanische Staatsanleihen von 1970-2009",
         cex.main=1.8,
         cex.lab=1,
         cex.axis=1,
         expand = c(0.4,0.6,1)
  )
}
jungbacker.data = daten()$jungbacker.data
pdf('Jungbacker-3DData.pdf', width=12, height=6)
Abbildung5.1(jungbacker.data)
dev.off()
### Abbildung 5.1.1: Screeplot ####
eigenwerte = faktoranalyse(yields, r = 17)$eigenval
pdf('Screeplot.pdf', width=6, height=4.2, pointsize = 10)
plot(x = 1:17, y = eigenwerte, ylab = "", xlab = "Eigenwerte der zentrierten Kovarianzmatrix", 
     type = "b", main = "Abbildung 5.1.1: Screeplot", cex.main = 2, cex.lab = 1.5, cex.axis = 1.3)
dev.off()
# erklärte Varianz
eigenwerte3 = (eigenwerte[1]+eigenwerte[2]+eigenwerte[3]) / sum(eigenwerte)
eigenwerte2 = (eigenwerte[1]+eigenwerte[2]) / sum(eigenwerte)
eigenwerte1 = eigenwerte[1]/sum(eigenwerte)
### Abbildung 5.2: Faktorladungen ####
par(mfrow = c(1,3))
# f
loading.f = faktoranalyse(yields, r = 3)$eigenvec
loading.f[,1:3] = -loading.f[,1:3]
plot(x = 1:17, xaxt = "n",y = loading.f[,1], ylim = c(-0.5,0.6), type = "b", lty=1, lwd = 2.5,
     xlab = "Eigenvektoreintrag", ylab = "Faktorladung", cex.axis = 1.5, cex.lab = 1.7)
axis(1, at=1:17, labels = 1:17, cex.axis = 1.7)
lines(loading.f[,2], lty = 10, lwd = 2.3, type = "b")
lines(loading.f[,3], lty = 2, lwd = 2.3, type = "b")
legend(x = "topright", legend = c('Ladung 1', 'Ladung 2', 'Ladung 3'), lty = c(1,10,3), 
       lwd = c(2,1.7,1.7), cex = 1.3)
# dl
loading.dl = diebold.li(yields)$fine.grid
plot(loading.dl[,1], type = "l", ylim = c(0,1.05), lty=1, lwd = 2.5, xlab = "Laufzeit in Monaten (m)", 
     ylab = "Faktorladung", main = "Abbildung 5.2: Faktorladungen",
     cex.main = 2.3, cex.lab = 1.7, cex.axis = 1.7)
lines(loading.dl[,2], lty = 5, lwd = 2.3)
lines(loading.dl[,3], lty = 3, lwd = 2.3)
legend(x = "right", legend = c('Niveau', 'Steigung', 'Krümmung'), lty = c(1,5,3), 
       lwd = c(2,1.7,1.7), cex = 1.3)
# s
loading.s = svensson(yields)$fine.grid
plot(loading.s[,1], type = "l", ylim = c(0,1.05), lty=1, lwd = 2.5, xlab = "Laufzeit in Monaten (m)", 
     ylab = "Faktorladung", cex.main = 2, cex.lab = 1.7, cex.axis = 1.7)
lines(loading.s[,2], lty = 5, lwd = 2.3)
lines(loading.s[,3], lty = 3, lwd = 2.3)
lines(loading.s[,4], lty = 4, lwd = 2.3)
legend(x = "right", legend = c('Niveau', 'Steigung', 'Krümmung 1', 'Krümmung 2'), lty = c(1,5,3,4), 
       lwd = c(2,1.7,1.7,1.7), cex = 1.3)
par(mfrow = c(1,1))
### Abbildung 5.2.1: Modellvergleich ####
# f
f.monoton = faktoranalyse(yields)$zinsen.f[120,]
f.komplex = faktoranalyse(yields)$zinsen.f[256,]
# dl
dl.monoton = diebold.li(yields)$zinsen.dl[120,]
dl.komplex = diebold.li(yields)$zinsen.dl[256,]
# s
s.monoton = svensson(yields)$zinsen.s[120,]
s.komplex = svensson(yields)$zinsen.s[256,]
par(mfrow = c(3,2))
plot(x = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120), yields[256,], pch = 8, xlab = "Laufzeit in Monaten", 
     ylab = "Zins (Prozent)", main = "Zinsstruktur der Faktoranalyse vom April 2006", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 1.5)
lines(x  = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120), f.komplex, type = "b", lwd = 2)
plot(x = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120), yields[120,], pch = 8, ylim = c(5.6, 8)
     ,ylab = "Zins (Prozent)", xlab = "Laufzeit in Monaten", main = "Zinsstruktur der Faktoranalyse vom Dezember 1994",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(x  = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120), f.monoton, type = "b", lwd = 2)
plot(x = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120), yields[256,], pch = 8, xlab = "Laufzeit in Monaten",
     ylab = "Zins (Prozent)", main = "Zinsstruktur nach Diebold und Li vom April 2006",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(x = 3:120, y = dl.komplex, lwd = 2)
plot(x = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120), yields[120,], pch = 8, ylim = c(5.6, 8),
     ylab = "Zins (Prozent)", xlab = "Laufzeit in Monaten", main = "Zinsstruktur nach Diebold und Li vom Dezember 1994",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(x = 3:120, y = dl.monoton, lwd = 2)
plot(x = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120), yields[256,], pch = 8, xlab = "Laufzeit in Monaten",
     ylab = "Zins (Prozent)", main = "Zinsstruktur nach Svensson vom April 2006",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(x = 3:120, y = s.komplex, lwd = 2)
plot(x = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120), yields[120,], pch = 8, ylim = c(5.6, 8),
     xlab = "Laufzeit in Monaten", ylab = "Zins (Prozent)", main = "Zinsstruktur nach Svensson vom Dezember 1994",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(x = 3:120, y = s.monoton, lwd = 2)
legend(x = "bottomright", legend = c("echter Zins", "Zinsstrukturenkurve"), pch = c(8,NA), lty = c(NA, 1), 
       lwd = c(1,2), cex = 1.7)
par(mfrow = c(1,1))










