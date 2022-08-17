## Cindy Data: processing ##
library(R2OpenBUGS)
library(fUnitRoots)
library(coda)
library(tseries)
# packet data stores the demographic information
pdata.all <- read.table('packet_data.txt', as.is = T)
pdata.all[pdata.all=="."] <- NA
pdata.id <- unique(pdata.all[,1])
# daily data
ddata.all <- read.table('daily_data.txt',as.is = T)
ddata.all[ddata.all=="."] <- NA
ddata.id <- unique(ddata.all[,1])
index <- which(pdata.id %in% ddata.id)
id <- pdata.id[index] # the ids in both packet data and daily data

age <- pdata.all[pdata.id %in% id, 4]
gender <- pdata.all[pdata.id %in% id,2]
race <- pdata.all[pdata.id %in% id,5]
edu <- pdata.all[pdata.id %in% id,6]

ddata <- ddata.all[ddata.id %in% id,]
#race.long <- rep(race,each=56)
#edu.long <- rep(edu,each=56)
#ddata <- data.frame(ddata,race.long,edu.long)
#colnames(ddata) <- c('id','day','NA1','PA1','Str1','NA3','PA3','Str3','NA5','PA5','Str5','race','edu')
#write.table(ddata,'ddata.txt',row.names=F,col.names=F)

## NA, PA, and stress at the first burst ##
NA1 <- as.vector(ddata[,3])
PA1 <- as.vector(ddata[,4])
str1 <- as.vector(ddata[,5])
NA1t <- PA1t <- str1t <- NULL

# 56 is the number of days measured for each individual; long format to wide format;
for (i in 1:length(id)) {
  
  NA1t <- rbind(NA1t, t(NA1[(56*i-55):(56*i)]))
  PA1t <- rbind(PA1t, t(PA1[(56*i-55):(56*i)]))
  str1t <- rbind(str1t, t(str1[(56*i-55):(56*i)]))
  
}

# choose only cases with complete data
dat.use <- cbind(age, gender, race, edu, NA1t, PA1t, str1t)
dat.use <- dat.use[complete.cases(dat.use),]
age <- as.numeric(dat.use[,1])
gender <- as.numeric(dat.use[,2])
race <- as.numeric(dat.use[,3])
edu <- as.numeric(dat.use[,4])
y <- apply(dat.use[,5:60],2,as.numeric)
y <- as.matrix(as.data.frame(y))
m <- apply(dat.use[,61:116],2,as.numeric)
m <- as.matrix(as.data.frame(m))
x <- apply(dat.use[,117:172],2,as.numeric)
x <- as.matrix(as.data.frame(x))

#KPSS test of trend stationary for original data;
p.x=p.m=p.y=NULL
for (i in 1:nrow(x)) {
  p.x = c(p.x, suppressWarnings(kpss.test(x[i,],null="Trend"))$p.value)
  p.m = c(p.m, suppressWarnings(kpss.test(m[i,],null="Trend"))$p.value)
  p.y = c(p.y, suppressWarnings(kpss.test(y[i,],null="Trend"))$p.value)
}

# detrend x, m, and y
t <- 1:56
res.x <- res.m <- res.y <- NULL
for (j in 1:nrow(x)) {
  res.x <- rbind(res.x, residuals(lm(x[j,]~t)))
  res.m <- rbind(res.m, residuals(lm(m[j,]~t)))
  res.y <- rbind(res.x, residuals(lm(y[j,]~t)))
}

#KPSS test of level stationary for detrended data;
p.rx=p.rm=p.ry=NULL
for (i in 1:nrow(x)) {
  p.rx = c(p.rx, suppressWarnings(kpss.test(res.x[i,]))$p.value)
  p.rm = c(p.rm, suppressWarnings(kpss.test(res.m[i,]))$p.value)
  p.ry = c(p.ry, suppressWarnings(kpss.test(res.y[i,]))$p.value)
}

x.nonstat <- which(p.rx<.05)
m.nonstat <- which(p.rm<.05)
y.nonstat <- which(p.ry<.05)
ind.nonstat <- unique(c(x.nonstat,m.nonstat,y.nonstat))

# Use individuals that show stationarity with detrended data;
x <- res.x[-ind.nonstat,]
m <- res.m[-ind.nonstat,]
y <- res.y[-ind.nonstat,]

n <- nrow(x)
t <- 56
Z <- gender

## Bayesian estimation 
R <- diag(6)
B <- 10000 # burn-in period
iter <- 20000 # number of iterations
dataun <- list("n","t", "x", "m", "y","R","Z")

# generate initial values
initun <-function() {
inits <- list(taux=.1,taum=.1,tauy=.1)
}

# model specification

model <-function() {

# model's likelihood
for (i in 1:n) {
	for (j  in 2:t) {
		x[i, j] ~ dnorm (mux[i,j] , taux)
		mux[i,j] <- b[i,1]*x[i,j-1]
		m[i, j] ~ dnorm(mum[i,j], taum)
		mum[i,j] <- b[i,2]*m[i,j-1]+b[i,4]*x[i, j-1]
		y[i, j] ~ dnorm(muy[i,j], tauy)
		muy[i,j] <- b[i,3]*y[i,j-1]+b[i,5]*m[i,j-1]+b[i,6]*x[i, j-1]
	}
}

for (i  in 1:n)
{
b[i,1:6] ~ dmnorm(mur[i,1:6], taur[1:6, 1:6])
mur[i,1] <- dx
mur[i,2] <- dm
mur[i,3] <- dy
mur[i,4] <- a1+a2*Z[i]
mur[i,5] <- b1+b2*Z[i]
mur[i,6] <- C1+C2*Z[i]
}

# priors
dx ~ dunif(-1,1)
dm~ dunif(-1,1)
dy ~ dunif(-1,1)
a1 ~ dunif(-3,3)
a2 ~ dunif(-3,3)
b1 ~ dunif(-3,3)
b2 ~ dunif(-3,3)
C1 ~ dunif(-3,3)
C2 ~ dunif(-3,3)
taux ~ dgamma(.01, .01)
taum ~ dgamma(.01, .01)
tauy ~ dgamma(.01, .01)
taur[1:6, 1:6] ~ dwish(R[1:6,1:6], 6)

# parameters of interest
sigma2r[1:6,1:6] <- inverse(taur[1:6, 1:6])
par[1] <- 1/taux
par[2] <- 1/taum
par[3] <- 1/tauy
par[4] <- dx # dx
par[5] <- dm # dm
par[6] <- dy # dy
par[7] <- a1 # a1
par[8] <- a2 # a2
par[9] <- b1 # b1
par[10] <- b2 # b2
par[11] <- C1 # c1
par[12] <- C2 # c2
par[13] <- sigma2r[1,1] # sigma(dxi)
par[14] <- sigma2r[2,2] # sigma(dmi)
par[15] <- sigma2r[3,3] # sigma^2(dyi)
par[16] <- sigma2r[4,4] # sigma^2(a)
par[17] <- sigma2r[5,5] # sigma^2(b)
par[18] <- sigma2r[6,6] # sigma^2(c)
par[19] <- sigma2r[1,2] # cov(dxi, dmi)
par[20] <- sigma2r[1,3] # cov(dxi, dyi)
par[21] <- sigma2r[1,4] # cov(dxi, ai)
par[22] <- sigma2r[1,5] # cov(dxi, bi)
par[23] <- sigma2r[1,6] # cov(dxi, ci)
par[24] <- sigma2r[2,3] # cov(dmi, dyi)
par[25] <- sigma2r[2,4] # cov(dmi, ai)
par[26] <- sigma2r[2,5] # cov(dmi, bi)
par[27] <- sigma2r[2,6] # cov(dmi, ci)
par[28] <- sigma2r[3,4] # cov(dyi, ai)
par[29] <- sigma2r[3,5] # cov(dyi, bi)
par[30] <- sigma2r[3,6] # cov(dyi, ci)
par[31] <- sigma2r[4,5] # cov(ai, bi)
par[32] <- sigma2r[4,6] # cov(ai, ci)
par[33] <- sigma2r[5,6] # cov(bi, ci)
par[34] <- (a1+a2)*(b1+b2)+par[31] # E(ai*bi);
par[35] <- pow((b1+b2),2)*par[16] +pow((a1+a2),2)*par[17]+par[16]*par[17]+2*(a1+a2)*(b1+b2)*par[31]+pow(par[31],2) # var(ai*bi);
}
para <- c("par")

# run the analysis
output <- bugs(dataun, initun, para, model,
                   n.chains = 1, n.burnin = B, DIC=T, n.iter = iter)
results <- output$summary
mcmc <- output$sims.matrix
write.table(results,'est.txt',row.names=F, col.names=F)
write.table(mcmc,'mcmc.txt',row.names=F, col.names=F)
