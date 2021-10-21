# gradient function, takes beta_hat, and sigma matrix as input
compute_min <- function (beta_hat,sig, maxit=500, tol=1e-8){
  fn = function (b){
    return (as.numeric(
      t(b-beta_hat) %*% sig %*% (b-beta_hat)))
  }
  gr = function (b){
    return(as.numeric(2*sig%*%(b-beta_hat)))
  }
  beta_cur <- beta_hat
  direction_cur <- gr (beta_cur)
  
  # Line search objective
  ls_obj = function (alpha){
    return (fn(beta_cur - alpha*direction_cur))
  }
  
  alpha_opt <- optimize (ls_obj, interval=c(0,1))$minimum
  beta_next <- beta_cur - alpha_opt*gr(beta_cur)
  for (k in 1:maxit){
    direction_next <- gr (beta_next)
    dx = beta_next - beta_cur
    if (max(abs(dx)) < tol){
      break
    }
    dg <- direction_next - direction_cur
    alpha_opt <- sum(dg*dx)/sum(dg^2)
    beta_cur <- beta_next
    direction_cur <- direction_next
    beta_next <- beta_next - alpha_opt*direction_next
  }
  return (beta_next)
  
  
  
}

###################
# this is the general case 
# this assumes x is defined correctly 
# we assume a specific structure of omega
# however each omega will be different 
p <- as.numeric(readline("enter value for p  "))
y <- as.numeric(readline("enter scalar for y  "))
B <- as.numeric(readline("enter large integer for B  "))
genav <- function(x){
  bomg <- 0 
  for (i in 1:B){
    # we draw K in a uniform distribution centered at log(N)
    k <- runif(1,1,2*ceiling(log(p)))
    omg <- matrix(0,nrow = k,ncol = p)
    for (c in 1:p){
      for (d in 1:k){
        omg[d,c] <- rnorm(1,0,1)
      }
    }
    q <- omg%*%x
    bt <- t(solve(t(q)%*%q)%*%t(q)*y)
    sig <- q%*%t(q)
    bomg1 <- t(omg)%*%compute_min(bt,sig)
    bomg <- (bomg*(i-1))/i + bomg1/i
    rm(omg)
  }
  return(bomg)
}
far <- matrix(0,nrow=2,ncol = 2)
dim(far)
dim

###########

#example with a specific x, omega, y, B, p 
# B from the problem 
B <- 1000
# bomg is running average, see code from #2 
bomg <- 0 
# p is fixed far greater than what k could be 
p <- 100
# y is a scalar 
y <- 10
for (i in 1:B){
  # we draw K in a uniform distribution centered at log(N), per Dr. Laber 
  k <- runif(1,1,10)
  x <- as.matrix(runif(p,5,15))
  omg <- matrix(0,nrow = k,ncol = p)
  for (c in 1:p){
    for (d in 1:k){
      omg[d,c] <- rnorm(1,0,1)
    }
  }
  q <- omg%*%x
  bt <- t(solve(t(q)%*%q)%*%t(q)*y)
  sig <- q%*%t(q)
  bomg1 <- t(omg)%*%compute_min(bt,sig)
  bomg <- (bomg*(i-1))/i + bomg1/i
  rm(omg)
}
# the next line returns the target 
bomg