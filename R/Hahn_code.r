require(glmnet)
require(ROCR)
require(survival)
require(data.table)


# ***********************************************
# simulation of genotype and methylation datasets
# ***********************************************

# https://en.wikipedia.org/wiki/Balding%E2%80%93Nichols_model
drawFreq <- function(Fst=0.001,maf_min=0,maf_max=0.1) {
	maf_value <- runif(1,min=maf_min,max=maf_max)
	freq <- rbeta(1, maf_value*(1-Fst)/Fst, (1-maf_value)*(1-Fst)/Fst)
	return(freq)
}

# simulate genotype and methylation datasets
simdata <- function(n=100,nGeno=100,nMeth=100,overlap=50,h=0.3,R2=c(0,0.1),ageRange=c(50,100)) {
	vareps <- 1-h**2
	g <- NULL
	m <- NULL
	counter <- 0
	# draw geno and y vector within given R2 range
	while(counter<overlap) {
		while(TRUE) {
			maf <- drawFreq()
			genos <- rbinom(n,2,p=maf)
			if(any(genos!=0)) break
		}
		betas <- sqrt(h**2/(2*maf*(1-maf)))
		eps <- rnorm(n,0,sqrt(vareps))
		y <- betas*genos + eps
		cor2 <- cor(genos,y)**2
		if(is.na(cor2)) browser()
		if(R2[1]<=cor2 & cor2<R2[2]) {
			g <- cbind(g,genos)
			m <- cbind(m,y)
			counter <- counter+1
		}
	}
	# add genos
	for(i in 1:(nGeno-overlap)) {
		maf <- drawFreq()
		genos <- rbinom(n,2,p=maf)
		g <- cbind(g,genos)
	}
	g <- unname(g)
	# add methylation
	for(i in 1:(nMeth-overlap)) {
		eps <- rnorm(n,0,sqrt(vareps))
		m <- cbind(m,eps)
	}
	m <- unname(m)
	# generate response, add both g+m, compute row sums
	temp <- as.numeric(cbind(g,m) %*% runif(nGeno+nMeth))
	# generate age from "temp" by scaling to ageRange
	onset <- (temp-min(temp))/(max(temp)-min(temp)) * (ageRange[2]-ageRange[1]) + ageRange[1]
	T <- runif(n,min=ageRange[1],max=ageRange[2])
	y <- as.numeric(onset<=T)
	return(list(T=T,y=y,X1=g,X2=m,X=rbind(g,m)))
}



# *******************************************
# likelihood and survival functions, coxlasso
# *******************************************

# vector b (of length ncol(Z)), vector of times T, vector of censoring indicators delta, matrix Z with one row per subject
l1 <- function(b,T,delta,Z) {
	# indicator T matrix, where each row i contains the indicators T>=T[i]
	indicator <- matrix(T,byrow=TRUE,nrow=length(T),ncol=length(T))>=T
	# multiply each row with exp(b*Z[j,]) for j=1,2,...
	temp <- t( t(indicator) * exp(as.numeric(Z %*% b)) )
	res <- delta * (as.numeric(Z %*% b) - log(rowSums(temp)))
	return(sum(res))
}

# regularized lasso as in the paper of Zhang and Lu (2007)
f1 <- function(b,T,delta,Z,lambda,b_hat) {
	n <- length(T)
	-1/n * l1(b,T,delta,Z) + lambda * sum(abs(b)/abs(b_hat))
}

# numerical integration
integral <- function(x,y) {
	res <- numeric(length(x))
	for(i in 2:length(x)) {
		res[i] <- (x[i]-x[i-1])*min(y[i],y[i-1]) + 0.5*(x[i]-x[i-1])*abs(y[i]-y[i-1])
	}
	return(res)
}

# normalize vector
vectornormalize <- function(x,option=1) {
	if(option==1) return( x/sqrt(sum(x**2)) )
	if(option==2) return( (x-min(x))/(max(x)-min(x)) )
}

# numerical differentiation to compute basehaz from integrated basehaz
numdiff <- function(m) {
	y <- m[,1]
	x <- m[,2]
	y <- diff(y)#/diff(x)
	x <- x[-length(x)]
	cbind(y,x)
}

# Cox-lasso of Hahn (2023)
coxlasso <- function(Z,T,delta,validation,lambda,maxit=1000,normalize=TRUE) {
	b_hat <- optim(par=runif(ncol(Z)), fn=function(b) -l1(b,T,delta,Z), method="SANN", gr=NULL, control=list(maxit=maxit))$par
	if(normalize) b_hat <- vectornormalize(b_hat)
	b <- optim(par=runif(ncol(Z)), fn=function(b) f1(b,T,delta,Z,lambda,b_hat), method="SANN", gr=NULL, control=list(maxit=maxit))$par
	if(normalize) b <- vectornormalize(b)
	# compute baseline hazard and multiply
	coxfit <- coxph(Surv(T,delta==1)~0)
	bhest <- numdiff(basehaz(coxfit))
	xv <- bhest[,2]
	hv <- sapply(1:nrow(validation), function(i) bhest[,1]*exp(sum(b*validation[i,])) )
	# calculate survival
	sv <- sapply(1:ncol(hv), function(i) exp(-cumsum(integral(xv,hv[,i]))) )
	return(list(xv=xv,hv=hv,sv=sv))
}


# ***************************************************
# survfit and Desikan
# ***************************************************

# converts model$strata to matrix of values
strataToMatrix <- function(x) {
	res <- NULL
	for(i in 1:length(x)) {
		s <- names(x[i])
		s <- paste0(s,",")
		# find equality sign and comma
		indices1 <- which(strsplit(s,"")[[1]]=="=")
		indices2 <- which(strsplit(s,"")[[1]]==",")
		# extract value combinations and store in matrix "res"
		temp <- NULL
		for(j in 1:length(indices1)) {
			temp <- c(temp,as.numeric(substr(s,indices1[j]+1,indices2[j]-1)))
		}
		res <- rbind(res,temp)
	}
	return(unname(res))
}

# expand a vector pair of t (times) and y (values) to some reference time tref
expandpair <- function(t,y,tref) {
	res <- NULL
	for(i in 1:length(t)) {
		w <- which(t[i]==tref)
		res <- c(res,rep(y[i],w))
		tref <- tref[-(1:w)]
	}
	res <- c(res,rep(0,length(tref)))
	return(res)
}

# find closest row in matrix to vector
closestRow <- function(X,v) {
	temp <- sapply(1:nrow(X), function(i) sum((X[i,]-v)**2) )
	which.min(temp)[1]
}

# survfit
survfitpackage <- function(Z,T,delta,validation) {
	inputdata <- data.frame(Z1=Z[,1],Z2=Z[,2],Z3=Z[,3],Z4=Z[,4],Z5=Z[,5],Z6=Z[,6],Z7=Z[,7],Z8=Z[,8],Z9=Z[,9],Z10=Z[,10],Z11=Z[,11],Z12=Z[,12])
	n <- 2
	model <- survfit(Surv(T,delta==1)~Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11+Z12,data=inputdata,type="kaplan-meier")
	m <- strataToMatrix(model$strata)
	xv <- sort(unique(model$time))
	sv <- matrix(0,nrow=nrow(validation),ncol=length(xv))
	for(i in 1:nrow(validation)) {
		# find model
		index <- closestRow(m,validation[i,1:n])
		sv[i,] <- expandpair(model[index]$time,model[index]$surv,xv)
	}
	return(list(xv=xv,sv=t(sv)))
}

# method of Desikan (2017)
desikan <- function(Z,T,delta,validation,maxit=1000,normalize=TRUE) {
	# define likelihood
	logL <- function(b) {
		numerator <- as.numeric(Z %*% b)
		# indicator T matrix, where each row i contains the indicators T<=T[i]
		indicator <- matrix(T,byrow=TRUE,nrow=length(T),ncol=length(T))<=T
		# multiply each row with exp(b*Z[j,]) for j=1,2,...
		temp <- t( t(indicator) * exp(as.numeric(Z %*% b)) )
		res <- numerator - log(rowSums(temp))
		return(sum(res))
	}
	# fit coefficients b
	fn <- function(z) logL(z)
	b <- optim(par=runif(ncol(Z)), fn=fn, method="SANN", gr=NULL, control=list(maxit=maxit))$par
	if(normalize) b <- vectornormalize(b)
	# compute baseline hazard and multiply
	coxfit <- coxph(Surv(T,delta==1)~0)
	bhest <- numdiff(basehaz(coxfit))
	xv <- bhest[,2]
	hv <- sapply(1:nrow(validation), function(i) bhest[,1]*exp(sum(b*validation[i,])) )
	# calculate survival
	sv <- sapply(1:ncol(hv), function(i) exp(-cumsum(integral(xv,hv[,i]))) )
	return(list(xv=xv,hv=hv,sv=sv))
}

# find closest index to given value in vector
closestIndex <- function(x,value) {
	which.min(abs(x-value))[1]
}



# ***************************************************
# cooperative learning
# ***************************************************

# likelihood of the cooperative cox-lasso
f2 <- function(b,T,delta,Z1,Z2,Z3,lambda,b_hat,alpha) {
	n <- length(T)
	indices1 <- 1:ncol(Z1)
	b1 <- b[indices1]
	indices2 <- (ncol(Z1)+1):(ncol(Z1)+ncol(Z2))
	b2 <- b[indices2]
	if(is.null(Z3)) {
		return( -1/n*l1(b1,T,delta,Z1) -1/n*l1(b2,T,delta,Z2) + lambda*sum(abs(b)/abs(b_hat)) + alpha*sum((Z1 %*% b1 - Z2 %*% b2)**2) )
	}
	else {
		indices3 <- (ncol(Z1)+ncol(Z2)+1):(ncol(Z1)+ncol(Z2)+ncol(Z3))
		b3 <- b[indices3]
		return( -1/n*l1(b1,T,delta,Z1) -1/n*l1(b2,T,delta,Z2) -1/n*l1(b3,T,delta,Z3) + lambda*sum(abs(b)/abs(b_hat)) + alpha*sum((Z1 %*% b1 - Z2 %*% b2)**2) + alpha*sum((Z2 %*% b2 - Z3 %*% b3)**2) )
	}
}

# full cooperative cox-lasso
cooplasso <- function(Z1,Z2,Z3=NULL,T,validation,delta,lambda,alpha,maxit=1000,normalize=TRUE) {
	Z <- cbind(Z1,Z2,Z3)
	b_hat <- optim(par=runif(ncol(Z)), fn=function(b) -l1(b,T,delta,Z), method="SANN", gr=NULL, control=list(maxit=maxit/10))$par
	if(normalize) b_hat <- vectornormalize(b_hat)
	b <- optim(par=runif(ncol(Z)), fn=function(b) f2(b,T,delta,Z1,Z2,Z3,lambda,b_hat,alpha), method="SANN", gr=NULL, control=list(maxit=maxit))$par
	if(normalize) b <- vectornormalize(b)
	# compute baseline hazard and multiply
	coxfit <- coxph(Surv(T,delta==1)~0)
	bhest <- numdiff(basehaz(coxfit))
	xv <- bhest[,2]
	hv <- sapply(1:nrow(validation), function(i) bhest[,1]*exp(sum(b*validation[i,])) )
	# calculate survival
	sv <- sapply(1:ncol(hv), function(i) exp(-cumsum(integral(xv,hv[,i]))) )
	return(list(xv=xv,hv=hv,sv=sv, b=b))
}

# neural net via Python3 script
neuralnet <- function(X,T,validation,delta) {
	write.table(X,"temp1.csv", row.names=FALSE, col.names=FALSE)
	write.table(T,"temp2.csv", row.names=FALSE, col.names=FALSE)
	write.table(validation,"temp3.csv", row.names=FALSE, col.names=FALSE)
	system("python3 neuralnet1.py")
	res <- fread("temp4.csv")$V1
	res <- vectornormalize(res)
	# compute baseline hazard and multiply
	coxfit <- coxph(Surv(T,delta==1)~0)
	bhest <- numdiff(basehaz(coxfit))
	xv <- bhest[,2]
	hv <- sapply(1:nrow(validation), function(i) bhest[,1]*exp(sum(res[i])) )
	# calculate survival
	sv <- sapply(1:ncol(hv), function(i) exp(-cumsum(integral(xv,hv[,i]))) )
	return(list(xv=xv,hv=hv,sv=sv))
}

# sksurv via Python3 script
sksurv <- function(X,T,validation,delta) {
	write.table(X,"temp1.csv", row.names=FALSE, col.names=FALSE)
	write.table(cbind(delta,T),"temp2.csv", row.names=FALSE, col.names=FALSE)
	write.table(validation,"temp3.csv", row.names=FALSE, col.names=FALSE)
	system("python3 neuralnet3.py")
	res <- read.table("temp4.csv",header=FALSE)
	res <- as.matrix(unname(res))
	# check
	n <- nrow(res)/3
	xv <- res[1,]
	for(i in 1:n) {
		if(!all(res[i,]==xv)) stop()
	}
	sv <- res[(n+1):(2*n),]
	hv <- res[(2*n+1):(3*n),]
	return(list(xv=xv,hv=t(hv),sv=t(sv)))
}



# *******************************************************
# cross validation to determine suitable lambda for Lasso
# *******************************************************

crossval <- function(T,y,X1,X2,X3=NULL,lambdaVector=seq(0.1,1,0.1),folds=10,maxit=100,alpha=0.5,makeplot=FALSE) {
	resL1 <- matrix(0,nrow=folds,ncol=length(lambdaVector))
	resAUC <- matrix(0,nrow=folds,ncol=length(lambdaVector))
	resC <- matrix(0,nrow=folds,ncol=length(lambdaVector))
	resIBS <- matrix(0,nrow=folds,ncol=length(lambdaVector))
	pool <- sample(1:nrow(X1))
	size <- nrow(X1) %/% folds
	for(f in 1:folds) {
		if(makeplot) print(c(f,folds))
		# sex, apoe, pcs + geno
		m <- min(size,length(pool))
		indices <- pool[1:m]
		pool <- pool[-(1:m)]

		Z1 <- X1[indices,]
		Z2 <- X2[indices,]
		Z3 <- X3[indices,]
		Ztrain <- cbind(Z1,Z2,Z3)
		Ttrain <- T[indices]
		Zvalid <- cbind(X1[-indices,],X2[-indices,],X3[-indices,])
		Tvalid <- T[-indices]
		# delta is I(T<=C) i.e. if failure time is before censoring time, meaning failure (Alzheimer) has occurred
		delta <- y[indices]
		status <- y[-indices]

		# loop
		for(i in 1:length(lambdaVector)) {
			lambda <- lambdaVector[i]
			if(makeplot) print(lambda)
			res <- cooplasso(Z1,Z2,Z3,Ttrain,Zvalid,delta,lambda=lambda,alpha=alpha,maxit=maxit)

			# get survival prediction
			temp <- numeric(nrow(Zvalid))
			for(j in 1:nrow(Zvalid)) {
				index <- closestIndex(res$xv,Tvalid[j])
				temp[j] <- res$sv[index,j]
			}

			# metrics
			indices <- is.na(temp)
			temp[indices] <- 0
			resL1[f,i] <- sum(abs(status-temp))/length(status)
			p <- prediction(temp,status)
			resAUC[f,i] <- performance(p,"auc")@y.values[[1]]
			resC[f,i] <- concordance(y~x,data=data.frame(y=status,x=temp))$concordance
			resIBS[f,i] <- sum((status-temp)**2)/length(status)
		}
	}
	metric <- colMeans(resIBS)
	if(makeplot) {
		plot(lambdaVector,metric,type="p",pch="+",xlab="lambda",ylab="Integrated Brier Score",xaxt="n",xlim=c(0,1))
		axis(side=1,at=(0:10)/10,labels=(0:10)/10)
	}
	else {
		return(lambdaVector[which.min(metric)])
	}
}



# *********************************************************************
# testing-validation split of input data and application of all methods
# *********************************************************************

maincomparison <- function(T,y,X1,X2,X3=NULL,lambda=NA,alpha=0.5,subsample=50,maxit=100,nMethods=5) {
	if(is.null(subsample)) indices <- 1:nrow(X)
	else indices <- sample(1:length(y), size=subsample, replace=FALSE)
	# prepare datasets
	Z1 <- X1[indices,]
	Z2 <- X2[indices,]
	Z3 <- X3[indices,]
	Ztrain <- cbind(Z1,Z2,Z3)
	Ttrain <- T[indices]
	Zvalid <- cbind(X1[-indices,],X2[-indices,],X3[-indices,])
	Tvalid <- T[-indices]
	# delta is I(T<=C) i.e. if failure time is before censoring time, meaning failure (Alzheimer) has occurred
	delta <- y[indices]
	status <- y[-indices]
	# Lasso, Desikan, Coop, Neuralnet, sksurv
	resTime <- numeric(nMethods)
	if(is.na(lambda)) lambda <- crossval(Ttrain,delta,Z1,Z2,Z3)
	resTime[1] <- system.time( resL <- coxlasso(Ztrain,Ttrain,delta,Zvalid,lambda=lambda,maxit=maxit) )[[3]]
	resTime[2] <- system.time( resD <- desikan(Ztrain,Ttrain,delta,Zvalid,maxit=maxit) )[[3]]
	resTime[3] <- system.time( resC <- cooplasso(Z1,Z2,Z3,Ttrain,Zvalid,delta,lambda=lambda,alpha=alpha,maxit=maxit) )[[3]]
	resTime[4] <- system.time( resN <- neuralnet(Ztrain,Ttrain,Zvalid,delta) )[[3]]
	resTime[5] <- system.time( resS <- sksurv(Ztrain,Ttrain,Zvalid,delta) )[[3]]
	# plot hazard or survival
	if(FALSE) {
		matplot(resL$xv,resL$hv,type="l",pch="+",xlab="age",ylab="hazard",col="black")
# 		matplot(resL$xv,resL$sv,type="l",pch="+",xlab="age",ylab="survival",col="black")
# 		matplot(resC$xv,resC$hv,type="l",pch="+",xlab="age",ylab="hazard",col="black")
# 		matplot(resC$xv,resC$sv,type="l",pch="+",xlab="age",ylab="survival",col="black")
# 		matplot(resD$xv,resD$hv,type="l",pch="+",xlab="age",ylab="hazard",col="black")
# 		matplot(resD$xv,resD$sv,type="l",pch="+",xlab="age",ylab="survival",col="black")
# 		matplot(resN$xv,resN$hv,type="l",pch="+",xlab="age",ylab="hazard",col="black")
# 		matplot(resN$xv,resN$sv,type="l",pch="+",xlab="age",ylab="survival",col="black")
# 		matplot(resS$xv,resS$hv,type="l",pch="+",xlab="age",ylab="hazard",col="black")
# 		matplot(resS$xv,resS$sv,type="l",pch="+",xlab="age",ylab="survival",col="black")
	 stop("End of script")
	}
	# evaluation metric
	temp <- matrix(0,nrow=nrow(Zvalid),ncol=nMethods)
	for(i in 1:nrow(Zvalid)) {
		# get closest indices for all three approaches
		indexL <- closestIndex(resL$xv,Tvalid[i])
		indexD <- closestIndex(resD$xv,Tvalid[i])
		indexC <- closestIndex(resC$xv,Tvalid[i])
		indexN <- closestIndex(resN$xv,Tvalid[i])
		indexS <- closestIndex(resS$xv,Tvalid[i])
		# get survival prediction
		temp[i,] <- c(resL$sv[indexL,i], resD$sv[indexD,i], resC$sv[indexC,i], resN$sv[indexN,i], resS$sv[indexS,i])
	}
	# L1, AUC, C, IBS
	resL1 <- numeric(ncol(temp))
	resAUC <- numeric(ncol(temp))
	resC <- numeric(ncol(temp))
	resIBS <- numeric(ncol(temp))
	for(i in 1:ncol(temp)) {
		indices <- is.na(temp[,i])
		temp[indices,i] <- 0
		resL1[i] <- sum(abs(status-temp[,i]))/length(status)
		p <- prediction(temp[,i],status)
		resAUC[i] <- performance(p,"auc")@y.values[[1]]
		resC[i] <- concordance(y~x,data=data.frame(y=status,x=temp[,i]))$concordance
		resIBS[i] <- sum((status-temp[,i])**2)/length(status)
	}
	return(list(resL1=resL1,resAUC=resAUC,resC=resC,resIBS=resIBS,resTime=resTime))
}
