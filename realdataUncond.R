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
dat.use <- cbind(gender, race, edu, NA1t, PA1t, str1t)
dat.use <- dat.use[complete.cases(dat.use),]
gender <- as.numeric(dat.use[,1])
race <- as.numeric(dat.use[,2])
edu <- as.numeric(dat.use[,3])
y <- apply(dat.use[,4:59],2,as.numeric)
y <- as.matrix(as.data.frame(y))
m <- apply(dat.use[,60:115],2,as.numeric)
m <- as.matrix(as.data.frame(m))
x <- apply(dat.use[,116:171],2,as.numeric)
x <- as.matrix(as.data.frame(x))

# person-mean centering;
onescol <- rep(1,ncol(x))
xctr <- as.matrix(apply(x,1,mean))
mctr <- as.matrix(apply(m,1,mean))
yctr <- as.matrix(apply(y,1,mean))
x <- x-xctr%*%onescol
m <- m-mctr%*%onescol
y <- y-yctr%*%onescol

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

## Bayesian estimation 
R <- diag(6)
B <- 10000 # burn-in period
iter <- 20000 # number of iterations
dataun <- list("n","t", "x", "m", "y","R")

# generate initial values
initun <-function() {
inits <- list(dx = .1, dm = .1, dy = .1, a=.1, b=.1, c=.1, taux=.1,taum=.1,tauy=.1)
}

# model specification

model <-function() {

# model's likelihood
for (i in 1:n) {
	for (j  in 2:t) {
		x[i, j] ~ dnorm (mux[i,j] , taux)
		mux[i,j] <- b[i,1]*x[i,j-1]
		m[i, j] ~ dnorm(mum[i,j], taum)
		mum[i,j] <- b[i,2]*m[i,j-1]+b[i, 4]*x[i, j-1]
		y[i, j] ~ dnorm(muy[i,j], tauy)
		muy[i,j] <- b[i,3]*y[i,j-1]+b[i,5]*m[i,j-1]+b[i,6]*x[i, j-1]
	}
}

for (i  in 1:n)
{
b[i,1:6] ~ dmnorm(mur[1:6], taur[1:6, 1:6])
}

# priors
mur[1] ~ dunif(-1,1)
mur[2] ~ dunif(-1,1)
mur[3] ~ dunif(-1,1)
mur[4] ~ dunif(-30,30)
mur[5] ~ dunif(-30,30)
mur[6] ~ dunif(-30,30)
taux ~ dgamma(.01, .01)
taum ~ dgamma(.01, .01)
tauy ~ dgamma(.01, .01)
taur[1:6, 1:6] ~ dwish(R[1:6,1:6], 6)

# parameters of interest
sigma2r[1:6,1:6] <- inverse(taur[1:6, 1:6])
par[1] <- 1/taux
par[2] <- 1/taum
par[3] <- 1/tauy
par[4] <- mur[1] # dx
par[5] <- mur[2] # dm
par[6] <- mur[3] # dy
par[7] <- mur[4] # a
par[8] <- mur[5] # b
par[9] <- mur[6] # c
par[10] <- sigma2r[1,1] # sigma(dxi)
par[11] <- sigma2r[2,2] # sigma(dmi)
par[12] <- sigma2r[3,3] # sigma^2(dyi)
par[13] <- sigma2r[4,4] # sigma^2(a)
par[14] <- sigma2r[5,5] # sigma^2(b)
par[15] <- sigma2r[6,6] # sigma^2(c)
par[16] <- sigma2r[1,2] # cov(dxi, dmi)
par[17] <- sigma2r[1,3] # cov(dxi, dyi)
par[18] <- sigma2r[1,4] # cov(dxi, ai)
par[19] <- sigma2r[1,5] # cov(dxi, bi)
par[20] <- sigma2r[1,6] # cov(dxi, ci)
par[21] <- sigma2r[2,3] # cov(dmi, dyi)
par[22] <- sigma2r[2,4] # cov(dmi, ai)
par[23] <- sigma2r[2,5] # cov(dmi, bi)
par[24] <- sigma2r[2,6] # cov(dmi, ci)
par[25] <- sigma2r[3,4] # cov(dyi, ai)
par[26] <- sigma2r[3,5] # cov(dyi, bi)
par[27] <- sigma2r[3,6] # cov(dyi, ci)
par[28] <- sigma2r[4,5] # cov(ai, bi)
par[29] <- sigma2r[4,6] # cov(ai, ci)
par[30] <- sigma2r[5,6] # cov(bi, ci)
par[31] <- par[7]*par[8]+par[28] # E(ai*bi);
for (i  in 1:n)
{
  med[i] <- b[i,4]*b[i,5]
}

}

para <- c("par","med")

# run the analysis
output <- bugs(dataun, initun, para, model,
                   n.chains = 1, n.burnin = B, DIC=T, n.iter = iter)
results <- output$summary
mcmc <- as.mcmc(output$sims.matrix)
hpd <- as.matrix(HPDinterval(mcmc[,1:125],prob=.95))
ci95 <- matrix(0,125,2)
# 95% percentile interval;
for(j in 1:125) {
ci95[j,] <- quantile(mcmc[,j],c(.025,.975))
}
NumNonsigMed <- sum((hpd[32:125,1]<0)&(hpd[32:125,2]>0)) # number of non-significant mediation effect among all the individuals;
TableParaIndex <- c(7,8,9,28,31,13,14,15)
Z <- geweke.diag(mcmc[,TableParaIndex])$z
write.table(Z,'Z.txt',row.names=F, col.names=F)
write.table(results[TableParaIndex,],'estUNMAMMpm.txt',row.names=F, col.names=F)
write.table(mcmc[TableParaIndex,],'mcmcUNMAMMpm.txt',row.names=F, col.names=F)
write.table(hpd[TableParaIndex,],'hpdUNMAMMpm.txt',row.names=F, col.names=F)
write.table(ci95[TableParaIndex,],'ci95UNMAMMpm.txt',row.names=F,col.names=F)