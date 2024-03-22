## Function to generate recruits from spawners
## spawn: number of spawners
## sigma: standard deviation of random error
## alpha: productivity parameter
## beta: density dependence parameter
## rho: temporal autocorrelation
## last.eps: recruitment residual year y-1
ricker<-function(spawn,sigma,alpha,beta,rho,last.eps) {
## normal random error with bias correction for the mean
delta<-rnorm(1,mean=-(sigma^2)/2,sd=sigma) 
## residual (year y) in log-space with temporal autocorrelation
eps<-rho*last.eps+sqrt(1-rho^2)*delta 
## ln(recruits/spawner) with autocorrelated residuals
lnrs<-log(alpha)-beta*spawn+eps 
## recruits
rec<-exp(lnrs)*spawn 
return(list=c(rec,eps))
}
