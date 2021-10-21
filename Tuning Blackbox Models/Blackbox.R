#install.packages("Rlab")
library(Rlab)
# these functions will allow the user to enter necessary values 
# make sure to run these one at a time otherwise R will put NA 
n <- as.numeric(readline("Enter an integer for number of Y obs  "))
p <- as.numeric(readline("Enter an integer for dimension of estimator  "))
K <- as.numeric(readline("Enter an integer for number of CV folds  "))
crit <- as.character(readline("Enter (in lowercase) mse or mae  "))
mmode <- as.numeric(readline("Enter 1 for Dropout, 2 for NoiseAdd, 3 for Robust  "))
if (mmode==1 || mmode==2){
  M <- as.numeric(readline("Enter an integer for number of simulations  "))
}
# This is the main function that will return a model based on the input
# set c to 0 if this isn't a robust model
# if robust set it to a px1 matrix (same as beta) or a constant
# set learning algorithm as function per the problem
claim <- function(X,Y,c,learna){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  c <- as.matrix(c)
  # calculating criteria 
  calcms <- function(X,Y,pfun){
    sse <- 0
    for (i in 1:n){
      pred <- pfun(X)
      sse <- sse+ (Y[i] - pred)^2
    }
    return(sse/n)
  }
  calcma <- function(X,Y,pfun){
    sad <- 0
    for (i in 1:n){
      pred <- pfun(X)
      sad <- sad + abs(Y[i]-pred)
    }
    return(sad/n)
  }
  #different modification/regularzation functions 
  modif <- function(X,mas){
    if (mmode==1){
      for (j in 1:M){
        for (i in 1:n){
          Z = rbern(1,1-mas)
          X[i,] <- (X[i,]*Z)/(1-mas)
        }
      }
      return(X)
    }
    if (mmode==2){
      for (j in 1:M){
        X <- X + rnorm(n,0,mas)
      }
      return(X)
    }
    if (mmode==3){
      X + mas
      return(X)
    }
  }
  # tuning function for noise and dropout 
  # params be a possible range of parameters (be a fx1 matrix for an integer f)
  tun <- function(learna,X,Y,params){
    # keeping track of the best performance as the loops iterate
    bestp <- 0 
    # setting high initial average so the first iteration will set the..
    # ..first benchmark for average performance
    besta <- 10e99
    for (j in 1:length(params)){
      para <- params[j]
      # storing variables
      currp <- 0
      #shuffle the data
      X <-X[sample(nrow(X)),]
      Y <- Y[sample(nrow(Y)),]
      
      # K folds
      folds <- cut(seq(1,nrow(X)),breaks=K,labels=FALSE)
      
      #K Cross-validation 
      for(i in 1:K){
        #Segmenting and modifying data 
        testi <- which(folds==i,arr.ind=TRUE)
        testx <- X[testi, ]
        trainx <- X[-testi, ]
        testy <- Y[testi, ]
        trainy <- Y[-testi, ]
        # imposing the modifying method
        testx <- modif(testx,para)
        trainx <- modif(trainx,para)
        #using alg to make prediction 
        predf <- learna(trainx,trainy)
        if (crit=="mse"){
          currp <- currp + calcms(testx,testy,predf)
        }else{
          currp <- currp + calcma(testx,testy,predf)
        }
      }
      # average performance 
      avp <- currp/K
      if (avp < besta){
        besta <- avp 
        bestp <- para
      }
    }
    return(bestp)
  }
  # robust tuning 
  robtun <- function(learna,X,Y){
    # storing variable 
    worstp <- 0
    deltmax <- matrix()
    # iterate this a high number of times since the process is random 
    for (iter in 1:200){
      deltamat <- matrix(,nrow = n,ncol = 1)
      for (d in 1:length(c)){
        colb <- c[d]
        # creating random distributional variables
        randv <- as.matrix(runif(n))
        rands <- sum(randv)
        scaledv <- randv/rands
        scaledv <- scaledv + (colb^2)
        for (i in 1:n){
          sqroot <- sqrt(randv[i])
          # randomly change sign 
          if (runif(1)<.5){
            sqroot <- sqroot*(-1)
          }
          deltamat[i] <- sqroot 
        }
      }
      # now basically do the reverse of the regular tuning function
      X <- modif(X)
      currp <- 0
      X <-X[sample(nrow(X)),]
      Y <- Y[sample(nrow(Y)),]
      
      # K folds
      folds <- cut(seq(1,nrow(X)),breaks=K,labels=FALSE)
      
      #K Cross-validation 
      for(i in 1:K){
        #Segmenting and modifying data 
        testi <- which(folds==i,arr.ind=TRUE)
        testx <- X[testi, ]
        trainx <- X[-testi, ]
        testy <- Y[testi, ]
        trainy <- Y[-testi, ]
        # imposing the modifying method
        testx <- modif(testx,deltamat)
        trainx <- modif(trainx,deltamat)
        #using alg to make prediction 
        predf <- learna(trainx,trainy)
        if (crit=="mse"){
          currp <- currp + calcms(testx,testy,predf)
        }else{
          currp <- currp + calcma(testx,testy,predf)
        }
      }
      # average performance 
      avp <- currp/K
      if (avp > worstp){
        worstp <- avp 
        deltmax <- deltamat
      }
    }
    return(deltmax)
  }
  tunep <- 0 
  # these if conditions will take a range of parameters and use the tuning..
  #..functions to collect the best one 
  if (mmode==1){
    paramr <- matrix(0,nrow = 20,ncol=1)
    for (i in 1:20){
      paramr[i] <- .5*i
    }
    tunep <- tun(learna,X,Y,paramr)
  }
  if (mmode==2){
    paramr <- matrix(0,nrow = 50,ncol = 1)
      for (i in 1:50){
        paramr[i] <- .1*i
      }
    tunep <- tune(learna,X,Y,paramr)
  }
  if (mmode==3){
    tunep <- robtun(learna,X,Y)
  }
  # return your given algo with the tuned parameter 
  return(learna(modif(X,tunep),Y))
}
