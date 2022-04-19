#Fitting (beta-)binomial AR(1) models to a given time series

#Time series about price stability in EA17 countries, i.e.,
n <- 17

#Transition probabilities of binomial AR(1) model:
tpbin <- function(k,l,p,rho){
	beta <- p*(1-rho)
	alpha <- beta+rho
	
    tp <- 0
    for(j in c(max(0,k+l-n):min(k,l))){
		tp <- tp + dbinom(j,l,alpha)*dbinom(k-j,n-l,beta)
    }
tp
}


#Transition probabilities of beta-binomial AR(1) model:
tpbbin <- function(k,l,p,rho,phi){
	beta <- p*(1-rho)
	alpha <- beta+rho
	phif <- (1-phi)/phi
	
    tp <- 0
    for(j in c(max(0,k+l-n):min(k,l))){
		tp <- tp + choose(l,j) * beta(j+phif*alpha, l-j+phif*(1-alpha)) / beta(phif*alpha, phif*(1-alpha))   *   choose(n-l,k-j) * beta(k-j+phif*beta, n-l-k+j+phif*(1-beta)) / beta(phif*beta, phif*(1-beta))
    }
tp
}


#Log-likelihood of binomial AR(1) model:
llbar1 <- function(par,data){
#par is vector (p,rho)
T <- length(data)
value <- -log(dbinom(data[1], n, par[1])) #full likelihood, otherwise use 0 here

for(t in c(2:T)) {
	value <- value-log(tpbin(data[t], data[t-1], par[1], par[2]))
}
value
}


#Log-likelihood of beta-binomial AR(1) model:
llbbar1 <- function(par,data){
#par is vector (p,rho,phi)
T <- length(data)

#For the first observation, we require the marginal distribution, which we have to compute numerically:

#Transition matrix of BBAR1 model:
tpm <- array(0,c(n+1,n+1))
for(k in c(0:n)){
	for(l in c(0:n)){
		tpm[k+1,l+1] <- tpbbin(k,l, par[1], par[2], par[3])
	}
}

#Vector of marginal probabilities has to be determined numerically from Markov property:
pmarg <- eigen(tpm, symmetric = FALSE)$vectors[,1]
pmarg <- Re(pmarg/sum(pmarg))

#Full likelihood:
value <- -log(pmarg[(data[1]+1)]) #Makes computations quite slow.
for (t in c(2:T)) {
	value <- value-log(tpbbin(data[t], data[t-1], par[1], par[2], par[3]))
}
value
}






#Data from Jan. 2000 to Aug. 2012:
data0 <- scan(file.choose()) #PriceStability.txt
Tlen0 <- length(data0) #152

#Data from Jan. 2000 to Dec. 2006:
Tlen <- 84
data <- data0[1:Tlen]

maxval <- max(data)
maxval #11

plot(data, type="b", pch=19, cex=0.5, xlab = "t", ylab = expression("Stability counts  x"[t]), cex.axis=0.85)

plot(data0, type="b", pch=19, cex=0.5, xlab = "t", ylab = expression("Stability counts  x"[t]), cex.axis=0.85, col=gray(0.5), xaxp=c(0,144,12), yaxp=c(0,15,3))
points(data, type="b", col=gray(0), pch=19, cex=0.5)
abline(v=Tlen+0.5, lty=2)


acf(data)

pacf(data)
#AR(1)-like dependence structure, with SACF(1) given by:
rho1 <- acf(data, plot=FALSE)[[1]][2]
rho1 #0.658172


#Analysis of marginal distribution
hist(data-0.5)

#Observations' mean
barX <- mean(data)
barX
#4.27381

#"Success" probability
pmm <- barX/n
pmm
#0.2514006

#Absolute frequencies
absfreq <- tabulate(data+1) #+1 to include 0
plot(0:maxval, absfreq/Tlen, type="h", xlab = "k", ylab = expression(paste("estimated P(X"[t],"=k)")), lwd=4, xlim=c(0,n), ylim=c(0,0.25))
points((0:n)+.25, dbinom(0:n, n, pmm), type="h", lwd=4, col=gray(0.5))


#Observations' variance
sX <- var(data)
sX
#4.924125

(Tlen-1)/Tlen*sX
#4.865505


#Testing for equidispersion:
IDbin <- n*(Tlen-1)/Tlen*sX/barX/(n-barX)
IDbin #1.520769

meanID <- 1-1/Tlen*(1-1/n)*(1+rho1)/(1-rho1)
meanID #0.9456482

sdID <- sqrt(2/Tlen*(1-1/n)*(1+rho1^2)/(1-rho1^2))
sdID #0.2380369

#Upper critical value:
meanID+qnorm(0.95)*sdID
#1.337184

#Upper-sided P-value:
1-pnorm(IDbin, mean=meanID, sd=sdID)
#0.007843893

#Zero frequency:
# p0 <- length(data[data==0])/Tlen
p0 <- absfreq[1]/Tlen
p0 #0.02380952

(1-pmm)^n
#0.007281846





#Binomial AR1
estml <- suppressWarnings(optim(c(0.5,0.65), llbar1, method="L-BFGS-B", lower=c(0.0001,0.0001), upper=c(0.9999,0.9999), control=list(ndeps=c(1e-4,1e-4)), data=foo[,1], hessian=TRUE))

pestml <- estml$par[[1]]
rhoestml <- estml$par[[2]]
ofiest <- estml$hessian #inverse covariance
neglmax <- estml$value
estcov <- solve(ofiest)

# #Estimates:
# c(pestml,rhoestml)
# #0.2555437 0.5784506
# 
# #Thinning parameters:
# betaestml <- pestml*(1-rhoestml)
# betaestml #0.1077243
# alphaestml <- betaestml+rhoestml
# alphaestml #0.6861749
# 
# (1-pestml)^n #0.006626226
# 
# #Estimated standard errors:
# c(sqrt(diag(estcov)))
# #0.02207082 0.05799249
# 
# #AIC and BIC:
# AIC <- 2*neglmax+2*2
# BIC <- 2*neglmax+log(Tlen)*2
# c(neglmax, AIC, BIC)
# #161.6514 327.3028 332.1645
# 
# 
# #Pearson residuals
# cMeanB <- function(l,p,rho){
# 	beta <- p*(1-rho)
# 	rho*l+n*beta
# }
# 
# cVarB <- function(l,p,rho){
# 	beta <- p*(1-rho)
# 	rho*(1-rho)*(1-2*p)*l+n*beta*(1-beta)
# }
# 
# 
# res <- (data[2:Tlen]-cMeanB(data[1:(Tlen-1)], pestml,rhoestml))/sqrt(cVarB(data[1:(Tlen-1)], pestml,rhoestml))
# acf(res)
# plot(res)
# mean(res)
# #-0.03048619
# var(res)
# #1.274702
# min(res)
# #-2.852998
max(res)
#3.717468



#PIT:

#Matrix of all CDFs:
# allcdfs <- array(0, c(maxval+2,maxval+1)) #additional zero row corresponding to F(-1)
# 
# for(l in c(0:maxval)){
# 	cpmf <- rep(0, (maxval+1))
# 	for(k in c(0:maxval)){
# 		cpmf[k+1] <- tpbin(k,l,pestml,rhoestml)
# 	}
# 	allcdfs[(2:(maxval+2)),l+1] <- cumsum(cpmf)
# }
# 
# nobins <- 5#10
# PIT <- array(0, c(2,nobins+1))
# 
# for(j in c(1:nobins)){
# 	u <- j/nobins
# 	pitval <- 0
# 	
# 	for(t in c(2:Tlen)){
# 		if(allcdfs[(data[t]+1), (data[t-1]+1)]<u){
# 			if(allcdfs[(data[t]+2), (data[t-1]+1)]<u){
# 				pitval <- pitval+1
# 			}else{
# 				pitval <- pitval+ (u-allcdfs[(data[t]+1), (data[t-1]+1)])/(allcdfs[(data[t]+2), (data[t-1]+1)]-allcdfs[(data[t]+1), (data[t-1]+1)])
# 			}
# 		}
# 	}
# 	PIT[1,j+1] <- pitval/(Tlen-1)
# 	PIT[2,j+1] <- PIT[1,j+1]-PIT[1,j]
# }
# PIT
#      # [,1]      [,2]       [,3]       [,4]       [,5]       [,6]      [,7]
# # [1,]    0 0.1276227 0.21309378 0.30103886 0.38374755 0.48191942 0.5918308
# # [2,]    0 0.1276227 0.08547107 0.08794508 0.08270869 0.09817187 0.1099114
#           # [,8]       [,9]      [,10]     [,11]
# # [1,] 0.7128437 0.79867680 0.88949965 1.0000000
# # [2,] 0.1210129 0.08583314 0.09082285 0.1105003
# 
# 
# PIT.freq= as.vector(rep(((1:nobins)-0.5)/nobins, PIT[2,2:(nobins+1)]*1000))
# PIT.hist <- hist(PIT.freq, plot=FALSE, breaks=nobins)
# PIT.hist$counts <- PIT.hist$counts/1000
# plot(PIT.hist, freq=TRUE, main="PIT histogram", ylab="PIT histogram", xlab="u", col="gray")
# 
# 
# 
# 
# #Forecasting based on the fitted model:
# #Last observation:
# data[Tlen] #7

#The h-step-ahead forecasting distribution is obtained from tppois(k,l,lambda,alpha) by inserting mu*(1-alpha^h) and alpha^h.

hmax <- 1
maxval <- 1
forecasts <- array(0, c(hmax+1,1+1))
cforecasts <- array(0, c(hmax+1,maxval+1))
for(h in c(1:hmax)){
	for(k in c(0:maxval)){
		forecasts[h,k+1] <- tpbin(k, foo[280, 1],pestml,rhoestml^h)
	}
	cforecasts[h,] <- cumsum(forecasts[h,])
	}
#Stationary marginal distribution:
forecasts[hmax+1,] <- dbinom(c(0:maxval),n,pestml)
cforecasts[hmax+1,] <- cumsum(forecasts[hmax+1,])


cforecasts[1,]
 # [1] 9.589683e-05 1.679414e-03 1.314191e-02 6.083263e-02 1.865306e-01
 # [6] 4.048929e-01 6.579032e-01 8.532021e-01 9.538778e-01 9.894025e-01
# [11] 9.982147e-01 9.997799e-01

cforecasts[hmax+1,]
 # [1] 0.006626226 0.045293293 0.151476844 0.333720865 0.552672086 0.748082074
 # [7] 0.882235816 0.954599967 0.985649819 0.996308060 0.999234918 0.999874259
 
plot(0:maxval,forecasts[1,], type="h", xlab = "k", ylab = expression(paste("P(X "[T+h],"=k | X"[T],"=7)")), lwd=4, ylim=c(0,0.35))
points((0:maxval)+.25,forecasts[hmax+1,], type="h", lwd=4, col=gray(0.5))







#Beta-binomial AR1
phimm <- 1/(1+(n-IDbin)/(IDbin-1)*(1-2*pmm*(1-pmm)*(1-rho1))/(1+rho1))
phimm #0.06017108

estml <- suppressWarnings(optim(c(pmm,rho1,phimm), llbbar1, method = "L-BFGS-B", lower=c(0.0001,0.0001,0.0001), upper=c(0.9999,0.9999,0.9999), control=list(ndeps=c(1e-4,1e-4,1e-4)), data=data, hessian=TRUE))

pestml <- estml$par[[1]]
rhoestml <- estml$par[[2]]
phiestml <- estml$par[[3]]
ofiest <- estml$hessian #inverse covariance
neglmax <- estml$value
estcov <- solve(ofiest)

#Estimates:
c(pestml,rhoestml,phiestml)
#0.26014861 0.62237572 0.03669981

#Estimated standard errors:
c(sqrt(diag(estcov)))
#0.02766035 0.06352620 0.02793572

#AIC and BIC:
AIC <- 2*neglmax+2*3
BIC <- 2*neglmax+log(Tlen)*3
c(neglmax, AIC, BIC)
#160.3385 326.6770 333.9695



#Properties of this model:
#Observations' mean
n*pestml #4.422526

#Observations' index of Dispersion:
1+ (n-1)*(1-2*pestml*(1-pestml)*(1-rhoestml)) / ( (1/phiestml-1)*(1+rhoestml) + 1-2*pestml*(1-pestml)*(1-rhoestml) ) #1.314791

#Thinning parameters:
betaestml <- pestml*(1-rhoestml)
betaestml #0.09823843
alphaestml <- betaestml+rhoestml
alphaestml #0.7206141

#Transition matrix of fitted model:
tpmest <- array(0,c(n+1,n+1))
for(k in c(0:n)){
	for(l in c(0:n)){
		tpmest[k+1,l+1] <- tpbbin(k,l, pestml,rhoestml,phiestml)
	}
}

#Vector of marginal probabilities has to be determined numerically from Markov property:
pmargest <- eigen(tpmest, symmetric = FALSE)$vectors[,1]
pmargest <- Re(pmargest/sum(pmargest))
 # [1] 1.261966e-02 5.326287e-02 1.153869e-01 1.699597e-01 1.901986e-01
 # [6] 1.712620e-01 1.282557e-01 8.146600e-02 4.437236e-02 2.081810e-02
# [11] 8.404629e-03 2.900655e-03 8.447237e-04 2.031052e-04 3.892053e-05
# [16] 5.599557e-06 5.399563e-07 2.627496e-08
freq <- rep(0, n+1)
freq[1:(maxval+1)] <- absfreq/Tlen
plot(0:n, freq, type="h", xlab = "k", ylab = expression(paste("estimated P(X"[t],"=k)")), lwd=4, ylim=c(0,0.3))
points((0:n)+.25, pmargest, type="h", lwd=4, col=gray(0.5))

#Observations' zero probability:
pmargest[1] #0.01261966


#Pearson residuals
cMeanBB <- function(l,p,rho){
	beta <- p*(1-rho)
	rho*l+n*beta
}

cVarBB <- function(l,p,rho,phi){
	beta <- p*(1-rho)
	alpha <- beta+rho
	phi*(alpha*(1-alpha)+beta*(1-beta))*l^2 + (rho*(1-rho)*(1-2*p)*(1-phi) - 2*n*beta*(1-beta)*phi)*l + n*beta*(1-beta)*(1+phi*(n-1))
}

res <- (data[2:Tlen]-cMeanBB(data[1:(Tlen-1)], pestml,rhoestml))/sqrt(cVarBB(data[1:(Tlen-1)], pestml,rhoestml,phiestml))
acf(res)
plot(res)
mean(res)
#-0.03693761
var(res)
#1.039959
min(res)
#-2.603386
max(res)
#3.367156



#PIT:

#Matrix of all CDFs:
allcdfs <- array(0, c(maxval+2,maxval+1)) #additional zero row corresponding to F(-1)

for(l in c(0:maxval)){
	cpmf <- rep(0, (maxval+1))
	for(k in c(0:maxval)){
		cpmf[k+1] <- tpbbin(k,l,pestml,rhoestml,phiestml)
	}
	allcdfs[(2:(maxval+2)),l+1] <- cumsum(cpmf)
}

nobins <- 5#10
PIT <- array(0, c(2,nobins+1))

for(j in c(1:nobins)){
	u <- j/nobins
	pitval <- 0
	
	for(t in c(2:Tlen)){
		if(allcdfs[(data[t]+1), (data[t-1]+1)]<u){
			if(allcdfs[(data[t]+2), (data[t-1]+1)]<u){
				pitval <- pitval+1
			}else{
				pitval <- pitval+ (u-allcdfs[(data[t]+1), (data[t-1]+1)])/(allcdfs[(data[t]+2), (data[t-1]+1)]-allcdfs[(data[t]+1), (data[t-1]+1)])
			}
		}
	}
	PIT[1,j+1] <- pitval/(Tlen-1)
	PIT[2,j+1] <- PIT[1,j+1]-PIT[1,j]
}
PIT
     # [,1]      [,2]       [,3]       [,4]       [,5]      [,6]      [,7]
# [1,]    0 0.1184296 0.18696447 0.27340969 0.36575405 0.4703274 0.5908643
# [2,]    0 0.1184296 0.06853489 0.08644522 0.09234436 0.1045733 0.1205369
          # [,8]       [,9]     [,10]      [,11]
# [1,] 0.7223583 0.81004776 0.9165624 1.00000000
# [2,] 0.1314940 0.08768949 0.1065146 0.08343762


PIT.freq= as.vector(rep(((1:nobins)-0.5)/nobins, PIT[2,2:(nobins+1)]*1000))
PIT.hist <- hist(PIT.freq, plot=FALSE, breaks=nobins)
PIT.hist$counts <- PIT.hist$counts/1000
plot(PIT.hist, freq=TRUE, main="PIT histogram", ylab="PIT histogram", xlab="u", col="gray")




#Forecasting based on the fitted model:
#Last observation:
data[Tlen] #7

#The h-step-ahead forecasting distribution is obtained by an MC approximation.

#Transition matrix of BBAR1 model:
tpmestml <- array(0,c(n+1,n+1))
for(k in c(0:n)){
	for(l in c(0:n)){
		tpmestml[k+1,l+1] <- tpbbin(k,l, pestml,rhoestml,phiestml)
	}
}

#Perron-Frobenius theory:
eigen(tpmestml)$values
# 1.000000000 0.622375718 0.397991434 ...
#rhoestml^c(0:n) #1.0000000000 0.6223757180 0.3873515344 ...
hmax <- 10
forecasts <- array(0, c(hmax+1,n+1))
cforecasts <- array(0, c(hmax+1,n+1))

tpmat <- diag(1, n+1)
for(h in c(1:hmax)){
	tpmat <- tpmat %*% tpmestml
	forecasts[h,] <- tpmat[,data[Tlen]+1]
	cforecasts[h,] <- cumsum(forecasts[h,])
	}
#Stationary marginal distribution:
forecasts[hmax+1,] <- pmargest
cforecasts[hmax+1,] <- cumsum(forecasts[hmax+1,])

cforecasts[1,]
# [1] 0.000251857 0.002951644 0.017141887 0.065010751 0.178722389 0.374534020
 # [7] 0.616118625 0.819550416 0.930237931 0.977116670 0.993576512 0.998465872
# [13] 0.999693936 0.999950694 0.999993956 0.999999497 0.999999979 1.000000000
cforecasts[2,]
 # [1] 0.00188692 0.01528405 0.06175966 0.16556736 0.33139798 0.53087721
 # [7] 0.71688118 0.85409670 0.93579508 0.97587416 0.99229054 0.99792718
# [13] 0.99954033 0.99991873 0.99998920 0.99999904 0.99999996 1.00000000
cforecasts[3,]
 # [1] 0.004542386 0.030256574 0.102427130 0.235462051 0.415526859 0.605292187
 # [7] 0.766602152 0.879759136 0.946249703 0.979263520 0.993155510 0.998093993
# [13] 0.999562004 0.999919833 0.999988990 0.999998987 0.999999953 1.000000000
cforecasts[4,]
 # [1] 0.007070266 0.042372774 0.130970685 0.279222110 0.464325283 0.647267885
 # [7] 0.795464785 0.895993600 0.953827841 0.982219723 0.994118970 0.998354199
# [13] 0.999619201 0.999929739 0.999990267 0.999999097 0.999999958 1.000000000
cforecasts[5,]
 # [1] 0.008987401 0.050821141 0.149579258 0.306396826 0.493768088 0.672370819
 # [7] 0.812897545 0.906060114 0.958716830 0.984225000 0.994811759 0.998553800
# [13] 0.999666238 0.999938508 0.999991487 0.999999210 0.999999963 1.000000000
cforecasts[6,]
 # [1] 0.01029958 0.05636679 0.16141290 0.32330274 0.51185175 0.68771800
 # [7] 0.82358165 0.91228103 0.96177742 0.98550112 0.99526107 0.99868595
# [13] 0.99969806 0.99994457 0.99999235 0.99999929 0.99999997 1.00000000
cforecasts[hmax+1,]
 # [1] 0.01261966 0.06588253 0.18126938 0.35122908 0.54142770 0.71268966
 # [7] 0.84094533 0.92241134 0.96678370 0.98760180 0.99600643 0.99890708
# [13] 0.99975181 0.99995491 0.99999383 0.99999943 0.99999997 1.00000000

plot(0:n,forecasts[1,], type="h", xlab = "k", ylab = expression(paste("P(X"[T+h],"=k | X"[T],"=7)")), lwd=4, ylim=c(0,0.35))
points((0:n)+.25,forecasts[hmax+1,], type="h", lwd=4, col=gray(0.5))

plot(0:n,forecasts[1,], type="h", xlab = "k", ylab = expression(paste("P(X"[T+h],"=k | X"[T],"=7)")), lwd=4, ylim=c(0,0.35))
points((0:n)+.25,forecasts[2,], type="h", lwd=4, col=gray(0.3))
points((0:n)+.5,forecasts[hmax+1,], type="h", lwd=4, col=gray(0.6))


#1-step-ahead quantiles using conditional cdfs:
ctpmestml <- tpmestml
for(l in c(0:n)){
	ctpmestml[,l+1] <- cumsum(tpmestml[,l+1])
}

#5%
q5 <- rep(0, n+1)
for(l in c(0:n)){
	qval <- 0
	while(ctpmestml[qval+1,l+1]<0.05){
		qval <- qval+1
		}
	q5[l+1] <- qval
}
q5
#[1] 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8

#50%
q50 <- rep(0, n+1)
for(l in c(0:n)){
	qval <- 0
	while(ctpmestml[qval+1,l+1]<0.5){
		qval <- qval+1
		}
	q50[l+1] <- qval
}
q50
#[1]  1  2  3  3  4  5  5  6  7  7  8  9  9 10 11 11 12 12

#95%
q95 <- rep(0, n+1)
for(l in c(0:n)){
	qval <- 0
	while(ctpmestml[qval+1,l+1]<0.95){
		qval <- qval+1
		}
	q95[l+1] <- qval
}
q95
#[1]  5  5  6  6  7  8  8  9  9 10 11 11 12 13 14 14 15 16

#Retrospective application:
plot(1:Tlen,data, type="b", pch=19, cex=0.5, xlab = "t", ylab = expression("Stability counts  x"[t]), cex.axis=0.85, axes = FALSE)
axis(side=1, at=seq(1,Tlen, by=12))
axis(side=2, at=c(seq(0,n, by=5),n))
points(2:Tlen, q50[data0[1:(Tlen-1)]+1], type="l", col=gray(0.5), lty=2)
points(2:Tlen, q5[data0[1:(Tlen-1)]+1], type="l", col=gray(0.5))
points(2:Tlen, q95[data0[1:(Tlen-1)]+1], type="l", col=gray(0.5))

#Forecasting:
plot(Tlen:Tlen0,data0[Tlen:Tlen0], type="b", pch=19, cex=0.5, xlab = "t", ylab = expression("Stability counts  x"[t]), cex.axis=0.85, xaxp=c(84,144,5), yaxp=c(0,15,3))
abline(v=Tlen+0.5, lty=2)
points((Tlen+1):Tlen0, q50[data0[Tlen:(Tlen0-1)]+1], type="l", col=gray(0.5), lty=2)
points((Tlen+1):Tlen0, q5[data0[Tlen:(Tlen0-1)]+1], type="l", col=gray(0.5))
points((Tlen+1):Tlen0, q95[data0[Tlen:(Tlen0-1)]+1], type="l", col=gray(0.5))



#Conditional probability for staying in zero:
tpmestml[1,1] #0.2563017

#Probability for 6 zeros in a run:
tpmestml[1,1]^5 #0.001106007
pmargest[1]*tpmestml[1,1]^5 #1.395742e-05


