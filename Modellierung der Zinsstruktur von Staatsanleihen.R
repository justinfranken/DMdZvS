###### Abbildungen ##################################################
###### Abbildung 5.1 ########
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
         xlab = 'Time',
         ylab ='',
         zlab ='Yield (Percent)',
         ticktype = "detailed",
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
##### Abbildung 3.2 #########
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
##### eventuell noch Svensson
lambda1 = 0.0609
lambda2 = 0.5
Sline = matrix(nrow=120, ncol=4)
for (i in 1:120){ # 120 ist größte Beobachtete Laufzeit
  Sline[i, 1] <- 1
  Sline[i, 2] <- (1-exp(-lambda1*i))/(lambda1 * i)
  Sline[i, 3] <- (1-exp(-lambda1*i))/(lambda1 * i) - exp(- lambda1 * i)
  Sline[i, 4] <- (1-exp(-lambda2*i))/(lambda2 * i) - exp(- lambda2 * i)
}
plot(1:120, Sline[,1], cex.lab = 1.3, cex.axis = 1.3,
     type = "l", ylim = c(0,1.05), lty=1, lwd = 2, xlab = "Laufzeit in Monaten", ylab = "Faktorladung")
lines(1:120, Sline[,2], lty = 5, lwd = 2)
lines(1:120, Sline[,3], lty = 3, lwd = 2)
lines(1:120, Sline[,4], lty = 2, lwd = 2)
###### Daten ########################################################
## Makro Daten
temp = tempfile()
download.file("https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2021-10.csv",temp)
Makrodata = read.csv(temp)
unlink(temp)
## Zins Daten
temp <- tempfile()
download.file("http://qed.econ.queensu.ca/jae/2014-v29.1/jungbacker-koopman-van_der_wel/jkv-data.zip",temp)
data <- read.table(unz(temp, "UnsmFB_70-09.txt"))
unlink(temp)
jungbacker.data = ts(data, start = c(1970, 1), frequency = 12)
colnames(jungbacker.data) = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
rm(temp)
rm(data)
## CU,FFR, INFL
CU = Makrodata[c(-1,-755),"CUMFNS"]
CU = ts(CU, start = 1959, end = c(2021,9), frequency = 12)
FFR = Makrodata[c(-1,-755),"FEDFUNDS"]
FFR = ts(FFR, start = 1959, end = c(2021,9), frequency = 12)
Infl = Makrodata[c(-1,-755),"CPIAUCSL"]
INFL = numeric(753)
for (i in 1:753) {
  INFL[i+12] = (Infl[i+12]/Infl[i]*100)-100
}
INFL = ts(INFL, start = 1959, end = c(2021,9), frequency = 12)
## Daten ab 1985 bis 2007
CU = window(CU, start = 1985, end = 2007)
FFR = window(FFR, start = 1985, end = 2007)
INFL = window(INFL, start = 1985, end = 2007)
yields = window(jungbacker.data, start = 1985, end = 2007)
rm(Makrodata, Infl,i)
###### Zusammenfassungsstatistik ####################################
summary(yields)
summary(CU)
summary(FFR)
summary(INFL)
sd.yields = numeric(17)
for (i in 1:17) {
  sd.yields[i] = sd(yields[,i])
}
sd.yields
sd(CU)
sd(FFR)
sd(INFL)
###### Faktoranalyse ################################################
yield.means.jb = colMeans(yields) # means
centered_jungbacker = scale(yields, scale = FALSE) # Daten um mean bereinigt
yield.cov = cov(centered_jungbacker) # Kovarianz berechnet
yields.eigen = eigen(yield.cov) # Eigenvektor/values berechnet
yields.eigenvec = yields.eigen$vectors
yields.eigenval = yields.eigen$values
t(yields.eigenvec[,1]) %*% yields.eigenvec[,1] # Eigenvektoren auf 1 normiert
t(yields.eigenvec[,1]) %*% yields.eigenvec[,2] # und unabhängig voneinander
plot(yields.eigenval, type = "l") # Screeplot
plot(-yields.eigenvec[,1], ylim = c(-0.4,0.65), type = "l", lwd = 2)
lines(-yields.eigenvec[,2], lty = 3) # Ladungsplot (warum L2 und L3 umgekehrt?)
lines(-yields.eigenvec[,3], lty = 2) # warum L1 bei ca. -0.2 ?
r = 3 # Anzahl Faktoren
factors = matrix(0,nrow = dim(centered_jungbacker)[1], ncol = r)
for (i in 1:dim(centered_jungbacker)[1]){
  factors[i,] = lm(centered_jungbacker[i,] ~ yields.eigenvec[,1:r] - 1)$coefficients
} # Regression der Faktoren
yields.f = factors %*% t(yields.eigenvec[,1:r]) + yield.means.jb # yields
Faktoranalyse = function(x, r = 3){ # x: Zinsdaten , r : Anzahl an Faktoren
  yield.means.jb = colMeans(x)
  centered_jungbacker = scale(x, scale = FALSE)
  yield.cov = cov(centered_jungbacker)
  yields.eigen = eigen(yield.cov)
  yields.eigenvec = yields.eigen$vectors[,1:r]
  yields.eigenval = yields.eigen$values[1:r]
  factors = matrix(0,nrow = dim(centered_jungbacker)[1], ncol = r)
  for (i in 1:dim(centered_jungbacker)[1]){
    factors[i,] = lm(centered_jungbacker[i,] ~ yields.eigenvec - 1)$coefficients
  }
  yields.f = factors %*% t(yields.eigenvec) + yield.means.jb
  output = list("Faktoren" = factors, "Zinsen.f" = yields.f, "mu" = yield.means.jb, 
                "eigenvec" = yields.eigenvec, "eigenval" = yields.eigenval)
  return(output)
}
t = Faktoranalyse(yields, r = 17)
t$eigenval
sum.all = sum(t$eigenval)
VAR1 = t$eigenval[1]/sum.all
VAR1and2 = (t$eigenval[1]+t$eigenval[2])/sum.all
VAR1and2and3 = (t$eigenval[1]+t$eigenval[2]+t$eigenval[3])/sum.all
###### Nelson Siegel nach Diebold Li ################################
NS <- matrix(nrow=17, ncol=3)
laufzeit = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
lambda = 0.0609
for (i in 1:17){
  NS[i, 1] = 1
  NS[i, 2] = (1-exp(-lambda*laufzeit[i]))/(lambda * laufzeit[i])
  NS[i, 3] = (1-exp(-lambda*laufzeit[i]))/(lambda * laufzeit[i]) - exp(- lambda * laufzeit[i])}
beta.dl <- matrix(nrow=265, ncol = 3)
for (i in 1:265){
  beta.dl[i,] = lm(yields[i,] ~ NS - 1)$coefficients}
yields.dl = beta.dl %*% t(NS)
yields.dl = ts(yields.dl, start = 1970, frequency = 12)
colnames(yields.dl) = colnames(yields)
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
  yields.dl = beta.dl %*% t(NS)
  output = list("Zinsen.dl" = yields.dl, "beta.dl" = beta.dl, "NS.dl" = NS)
  return(output)
}
###### Nelson Siegel nach Svensson ##################################
lambda1 = 0.0609
lambda2 = 0.5
Sline = matrix(nrow=17, ncol=4)
laufzeit = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
for (i in 1:17){ # 120 ist größte Beobachtete Laufzeit
  Sline[i, 1] <- 1
  Sline[i, 2] <- (1-exp(-tau1*laufzeit[i]))/(tau1 * laufzeit[i])
  Sline[i, 3] <- (1-exp(-tau1*laufzeit[i]))/(tau1 * laufzeit[i]) - exp(- tau1 * laufzeit[i])
  Sline[i, 4] <- (1-exp(-tau2*laufzeit[i]))/(tau2 * laufzeit[i]) - exp(- tau2 * laufzeit[i])
}
beta.s = matrix(nrow=265, ncol = 4)
for (i in 1:265){
  beta.s[i,] = lm(yields[i,] ~ Sline - 1)$coefficients}
yields.s = beta.s %*% t(Sline)
yields.s = ts(yields.s, start = 1970, frequency = 12)
colnames(yields.s) = colnames(yields)
svensson = function(x, tau1 = 0.0609, tau2 = 0.5){
  Sline = matrix(nrow = 17, ncol = 4)
  laufzeit = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
  for (i in 1:17){
    Sline[i, 1] <- 1
    Sline[i, 2] <- (1-exp(-tau1*laufzeit[i]))/(tau1 * laufzeit[i])
    Sline[i, 3] <- (1-exp(-tau1*laufzeit[i]))/(tau1 * laufzeit[i]) - exp(- tau1 * laufzeit[i])
    Sline[i, 4] <- (1-exp(-tau2*laufzeit[i]))/(tau2 * laufzeit[i]) - exp(- tau2 * laufzeit[i])
  }
  beta.s = matrix(nrow=dim(x)[1], ncol = 4)
  for (i in 1:265){
    beta.s[i,] = lm(yields[i,] ~ Sline - 1)$coefficients}
  yields.s = beta.s %*% t(Sline)
  output = list("zinsen.s" = yields.s, "beta.s" = beta.s, "NS.s" = Sline)
  return(output)
}
###### Autoregressives Modell ohne Makrovariablen ###################
### Faktoranalyse
mean.factors = colMeans(factors)
for (i in 1:dim(factors)[1]){
  factors[i,] = factors[i,] - mean.factors
}
theta.f = numeric(r)
for (i in 1:r) {
  theta.f[i] = lm(embed(factors,2)[,i] ~ embed(factors, 2)[,i+r])$coefficients[2]
}
pred.factor.ar = theta.f * factors[dim(factors)[1],]
pred.yield.f.ar = yield.means.jb + pred.factor.ar %*% t(yields.eigenvec[,1:r])
### Diebold und Li
### Svensson
###### Vektorautoregressives Modell ohne Makrovariablen #############
### Faktoranalyse
mean.factors = colMeans(factors)
for (i in 1:dim(factors)[1]){
  factors[i,] = factors[i,] - mean.factors
}
Theta.f = lm(embed(factors,2)[,1:r] ~ embed(factors,2)[,-(1:r)])$coefficients[-1,]
pred.factor.var = Theta.f %*% factors[dim(factors)[1],]
pred.yield.f.var = yield.means.jb + t(pred.factor.var) %*% t(yields.eigenvec[,1:r])
### Diebold und Li
### Svensson
###### Autoregressives Modell mit Makrovariablen ####################
### Faktoranalyse
factors.makro = cbind(factors, CU, FFR, INFL)
mean.factors.makro = colMeans(factors.makro)
for (i in 1:dim(factors.makro)[1]){
  factors.makro[i,] = factors.makro[i,] - mean.factors.makro
}
theta.f.makro = numeric(r)
for (i in 1:r){
  theta.f.makro[i] = lm(embed(factors.makro,2)[,i] ~ -1 + embed(factors.makro,2)[,c(i+r+3,r+7,r+8,r+9)])$coefficients[1]
}
pred.factor.ar.makro = theta.f.makro * factors[dim(factors)[1],]
pred.yield.f.ar.makro = yield.means.jb + pred.factor.ar.makro %*% t(yields.eigenvec[,1:r])
### Diebold und Li
### Svensson
###### Vektorautoregressives Modell mit Makrovariablen ##############
### Faktoranalyse
factors.makro = cbind(factors, CU, FFR, INFL)
mean.factors.makro = colMeans(factors.makro)
for (i in 1:dim(factors.makro)[1]) {
  factors.makro[i,] = factors.makro[i,] - mean.factors.makro
}
Theta.f.makro = lm(embed(factors.makro,2)[,1:(r+3)] ~ -1 + embed(factors.makro,2)[,(r+4):(r+9)])$coefficients
pred.factor.var.makro = Theta.f.makro[1:r,1:r] %*% factors[dim(factors)[1],]
pred.yield.f.var.makro = yield.means.jb + t(pred.factor.var.makro) %*% t(yields.eigenvec[,1:r])
### Diebold und Li
### Svensson
###### Vorhersage - Funktionen ############################
### Faktoranalyse
# Es müssen 3 makrovariablen rein, wenn makro = TRUE!
prognose.faktoranalyse = function(x, r = 3, h = 1, VAR = FALSE, makro = FALSE, makro.var = cbind(CU,FFR,INFL)){
  if(r == 1)stop("r muss größer 1 sein!")
  factors = Faktoranalyse(x,r)$Faktoren
  mu = Faktoranalyse(x,r)$mu
  eigenvec = Faktoranalyse(x,r)$eigenvec
  if(makro == FALSE){
    mean.factors = colMeans(factors)
    for (i in 1:dim(factors)[1]){
      factors[i,] = factors[i,] - mean.factors
    }
    if(VAR == FALSE){
      theta.f = numeric(r)
      for (i in 1:r) {
        theta.f[i] = lm(embed(factors,2)[,i] ~ embed(factors, 2)[,i+r])$coefficients[2]
      }
      pred.fac = theta.f * factors[dim(factors)[1],]
      if(h > 1){
        pred.factor = matrix(nrow = h, ncol = r)
        pred.factor[1,] = pred.fac
        for (i in 1:(h-1)){
          pred.factor[i+1,] = theta.f * pred.factor[i,]
        }
      }
      else{pred.factor = pred.fac}
      pred.yield.f = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h){
          pred.yield.f[i,] = mu + pred.factor[i,] %*% t(eigenvec)
        }
      }
      else{pred.yield.f = mu + pred.factor %*% t(eigenvec)}
    }
    else{
      Theta.f = lm(embed(factors,2)[,1:r] ~ embed(factors,2)[,-(1:r)])$coefficients[-1,]
      dimnames(Theta.f) = list(1:dim(Theta.f)[1], 1:dim(Theta.f)[2])
      pred.fac = Theta.f %*% factors[dim(factors)[1],]
      if(h > 1){
        pred.factor = matrix(nrow = h, ncol = r)
        pred.factor[1,] = pred.fac
        for (i in 1:(h-1)) {
          pred.factor[i+1,] = Theta.f %*% pred.factor[i,]
        }
      }
      else{pred.factor = pred.fac}
      pred.yield.f = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h){
          pred.yield.f[i,] = mu + pred.factor[i,] %*% t(eigenvec)
        }
      }
      else{pred.yield.f = mu + pred.factor %*% t(eigenvec)}
    }
  }
  if(makro == TRUE){
    factors.makro = cbind(factors, makro.var)
    mean.factors.makro = colMeans(factors.makro)
    for (i in 1:dim(factors.makro)[1]){
      factors.makro[i,] = factors.makro[i,] - mean.factors.makro
    }
    if(VAR == FALSE){
      theta.f.makro = numeric(r)
      for (i in 1:r){
        theta.f.makro[i] = lm(embed(factors.makro,2)[,i] ~ -1 + embed(factors.makro,2)[,c(i+r+3, dim(embed(factors.makro,2))[2]-2,dim(embed(factors.makro,2))[2]-1,dim(embed(factors.makro,2))[2])])$coefficients[1]
      }
      pred.fac = theta.f.makro * factors[dim(factors)[1],]
      if(h > 1){
        pred.factor = matrix(nrow = h, ncol = r)
        pred.factor[1,] = pred.fac
        for (i in 1:(h-1)){
          pred.factor[i+1,] = theta.f.makro * pred.factor[i,]
        }
      }
      else{pred.factor = pred.fac}
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
      pred.fac = Theta.f.makro[1:r,1:r] %*% factors[dim(factors)[1],]
      if(h > 1){
        pred.factor = matrix(nrow = h, ncol = r)
        pred.factor[1,] = pred.fac
        for (i in 1:(h-1)) {
          pred.factor[i+1,] = Theta.f.makro[1:r,1:r] %*% pred.factor[i,]
        }
      }
      else{pred.factor = pred.fac}
      pred.yield.f = matrix(nrow = h, ncol = 17)
      if(h > 1){
        for (i in 1:h){
          pred.yield.f[i,] = mu + pred.factor[i,] %*% t(eigenvec)
        }
      }
      else{pred.yield.f = mu + pred.factor %*% t(eigenvec)}
    }
  }
  if(VAR == FALSE & makro == FALSE){theta.f = theta.f}
  if(VAR == TRUE & makro == FALSE){theta.f = Theta.f}
  if(VAR == FALSE & makro == TRUE){theta.f = theta.f.makro}
  if(VAR == TRUE & makro == TRUE){theta.f = Theta.f.makro}
  output = list("prog.zins" = pred.yield.f, "prog.faktor" = pred.factor, "theta.f" = theta.f)
  return(output)
}
###### RMSE ###############################################
### Faktoranalyse
# 
rmse.1.f = function(x,r,VAR = FALSE, makro = FALSE, start.training = 100){ 
  wiederholungen = dim(x)[1]-start.training
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    y_hat = prognose.faktoranalyse(trainings.data, VAR = VAR, makro = makro)$prog.zins
    y = yields[start.training+i,]
    pred.error[i] = mean(y-y_hat)^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
rmse.6.f = function(x,r,VAR = FALSE, makro = FALSE, start.training = 100){ 
  wiederholungen = dim(x)[1]-start.training-6
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    y_hat = prognose.faktoranalyse(trainings.data, VAR = VAR, makro = makro, h = 6)$prog.zins
    y = yields[start.training+i+5,]
    pred.error[i] = mean(y-y_hat[6,])^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
rmse.12.f = function(x,r,VAR = FALSE, makro = FALSE, start.training = 100){ 
  wiederholungen = dim(x)[1]-start.training-12
  pred.error = numeric(wiederholungen)
  for (i in 1:wiederholungen) {
    trainings.data = yields[1:(start.training+i-1),]
    y_hat = prognose.faktoranalyse(trainings.data, VAR = VAR, makro = makro, h = 12)$prog.zins
    y = yields[start.training+i+11,]
    pred.error[i] = mean(y-y_hat[12,])^{2}
  }
  mse = mean(pred.error)
  rmse = sqrt(mean(pred.error)) 
  output = list("rmse" = rmse, "mse" = mse, "pred.error" = pred.error)
  return(output)
}
### Diebold und Li
### Svensson













##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
###############################################################################
###############################################################################
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
###############################################################################
###############################################################################
# ----- Funktionen ----- ######################################################
###############################################################################
### --- Modellierung --- ##############
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
  yields.f = factors %*% t(yields.eigenvec[,1:r]) + yield.means.jb
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
svensson = function(x, tau1 = 0.0609, tau2 = 0.5){
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
### --- Vorhersage --- ################
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
      theta.f.makro = matrix(nrow = 1+r, ncol = 3)
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
          pred.beta.dl[i,] = pred.beta.dl[i,] + mean.beta.dl.makro[1:r]
        }
      }
      else{pred.beta.dl = pred.dl + mean.beta.dl.makro[1:r]}
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
prognose.svensson = function(x, h = 1, tau1 = 0.0609, tau2 = 0.5, VAR = FALSE, makro = FALSE, makro.var = cbind(CU,FFR,INFL)){
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
### --- RMSE --- ######################
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
rmse.1.s = function(x, tau1 = 0.0609, tau2 = 0.5, VAR = FALSE, makro = FALSE, start.training = 100){
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
rmse.6.s = function(x, tau1 = 0.0609, tau2 = 0.5, VAR = FALSE, makro = FALSE, start.training = 100){ 
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
rmse.12.s = function(x, tau1 = 0.0609, tau2 = 0.5, VAR = FALSE, makro = FALSE, start.training = 100){ 
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
### --- Tests --- #############################################################
# f
sink("RMSE-MSE-f.txt")
# VAR = FALSE, makro = FALSE
rmse.1.f(yields, VAR = FALSE, makro = FALSE) #mse 0.05836523, rmse 0.241589
rmse.6.f(yields, VAR = FALSE, makro = FALSE) #mse 0.5719575, rmse 0.7562787
rmse.12.f(yields, VAR = FALSE, makro = FALSE) #mse 1.411788, rmse 1.188187
# VAR = TRUE, makro = FALSE
rmse.1.f(yields, VAR = TRUE, makro = FALSE) #mse 0.06155781, rmse 0.2481085
rmse.6.f(yields, VAR = TRUE, makro = FALSE) #mse 0.6212923, rmse 0.788221
rmse.12.f(yields, VAR = TRUE, makro = FALSE) #mse 1.49806, rmse 1.223953
# VAR = FALSE, makro = TRUE
rmse.1.f(yields, VAR = FALSE, makro = TRUE) #mse 0.05668167, rmse 0.2380791
rmse.6.f(yields, VAR = FALSE, makro = TRUE) #mse 0.494982, rmse 0.7035496
rmse.12.f(yields, VAR = FALSE, makro = TRUE) #mse 1.098745, rmse 1.04821
# VAR = TRUE, makro = TRUE
rmse.1.f(yields, VAR = TRUE, makro = TRUE) #mse 0.08644342, rmse 0.2940126
rmse.6.f(yields, VAR = TRUE, makro = TRUE) #mse 1.040285, rmse 1.019944
rmse.12.f(yields, VAR = TRUE, makro = TRUE) #mse 2.362424, rmse 1.537018
sink()
# dl
sink("RMSE-MSE-dl.txt")
# VAR = FALSE, makro = FALSE
rmse.1.dl(yields, VAR = FALSE, makro = FALSE) #mse 0.06433368, rmse 0.2536409
rmse.6.dl(yields, VAR = FALSE, makro = FALSE) #mse 1.641818, rmse 1.281335
rmse.12.dl(yields, VAR = FALSE, makro = FALSE) #mse 5.359948, rmse 2.315156
# VAR = TRUE, makro = FALSE
rmse.1.dl(yields, VAR = TRUE, makro = FALSE) #mse 0.06425708, rmse 0.2534898
rmse.6.dl(yields, VAR = TRUE, makro = FALSE) #mse 2.910559, rmse 1.706036
rmse.12.dl(yields, VAR = TRUE, makro = FALSE) #mse 8.810468, rmse 2.968243
# VAR = FALSE, makro = TRUE
rmse.1.dl(yields, VAR = FALSE, makro = TRUE) #mse 0.1684597, rmse 0.4104385
rmse.6.dl(yields, VAR = FALSE, makro = TRUE) #mse 9.175558, rmse 3.029118
rmse.12.dl(yields, VAR = FALSE, makro = TRUE) #mse 19.24672, rmse 4.387109
# VAR = TRUE, makro = TRUE
rmse.1.dl(yields, VAR = TRUE, makro = TRUE) #mse 0.9318019, rmse 0.9652989
rmse.6.dl(yields, VAR = TRUE, makro = TRUE) #mse 35.89161, rmse 5.990961
rmse.12.dl(yields, VAR = TRUE, makro = TRUE) #mse 71.40152, rmse 8.449942
sink()
# dl
sink("RMSE-MSE-s.txt")
# VAR = FALSE, makro = FALSE
rmse.1.s(yields, VAR = FALSE, makro = FALSE) #mse 0.06472953, rmse 0.25442
rmse.6.s(yields, VAR = FALSE, makro = FALSE) #mse 1.151443, rmse 1.073053
rmse.12.s(yields, VAR = FALSE, makro = FALSE) #mse 3.789512, rmse 1.946667
# VAR = TRUE, makro = FALSE
rmse.1.s(yields, VAR = TRUE, makro = FALSE) #mse 0.1560899, rmse 0.3950822
rmse.6.s(yields, VAR = TRUE, makro = FALSE) #mse 3.322557, rmse 1.822788
rmse.12.s(yields, VAR = TRUE, makro = FALSE) #mse 8.140711, rmse 2.853193
# VAR = FALSE, makro = TRUE
rmse.1.s(yields, VAR = FALSE, makro = TRUE) #mse 0.1196051, rmse 0.3458396
rmse.6.s(yields, VAR = FALSE, makro = TRUE) #mse 6.607175, rmse 2.570443
rmse.12.s(yields, VAR = FALSE, makro = TRUE) #mse 15.93473, rmse 3.991833
# VAR = TRUE, makro = TRUE
rmse.1.s(yields, VAR = TRUE, makro = TRUE) #mse 1.469275, rmse 1.212136
rmse.6.s(yields, VAR = TRUE, makro = TRUE) #mse 34.01371, rmse 5.832127
rmse.12.s(yields, VAR = TRUE, makro = TRUE) #mse 133.9556, rmse 11.57392
sink()


### Warum die DieboldLi und Svensson Prognosen für Makrovariablen so schlecht sind (für h > 1):
dieboldli.6 = prognose.dieboldli(yields, makro = TRUE, h = 6)
dieboldli.6$theta.dl
dieboldli.6$prog.faktor

svensson.6 = prognose.svensson(yields, makro = TRUE, h = 6, VAR = TRUE)
svensson.6$theta.s
svensson.6$prog.faktor

faktoranalyse.6 = prognose.faktoranalyse(yields, makro = TRUE, h = 6)
faktoranalyse.6$theta.f
faktoranalyse.6$prog.faktor

plot(CU)
plot(FFR)
plot(INFL)







plot(yields[56,], type = "l")
plot(yields[60,], lwd = 1, ylim = c(7.6,7.95), pch = 4)
dl = diebold.li(yields)$zinsen.dl[60,]
p.dl = dl[c(1,4,7,10,13,16,19,22,28,34,46,58,70,82,94,106,118)]
lines(p.dl, lwd = 2, lty = 2)
s = svensson(yields)$zinsen.s[60,]
p.s = s[c(1,4,7,10,13,16,19,22,28,34,46,58,70,82,94,106,118)]
lines(p.s, lwd = 2, lty = 3)
f = faktoranalyse(yields)$zinsen.f[60,]
lines(f)
plot(f, type = "l")
lines(yields[60,], pch = 4)
