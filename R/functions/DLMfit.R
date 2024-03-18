##========================================================================##
##                                                                        ##
## Function for fitting DLMs with time-varying intercept/slope parameters ##
##                                                                        ##
##========================================================================##
## data: contains columns 'Esc' (escapement), 'Rec' (recruitment)
## var_alpha: if TRUE alpha estimated as time-varying parameter
## var_beta: if TRUE beta estimated as time-varying parameter
##========================================================================##
## applied to a linearized Ricker with alpha/beta parameters
DLMfit<-function(data,var_alpha,var_beta) {
lnRS<-log(data$Rec/data$Esc) ## ln(recruits/spawner)
alpha_y<-beta_y<-NULL
mod<-dlmModReg(data$Esc) ## specify a linear model
npara<-3+sum(c(var_alpha,var_beta)) ## number of parameters
build_mod<-function(parm) { ## build model with variance structure
mod$V<-exp(parm[1]) ## variance observation error 
mod$W[1,1]<-mod$W[2,2]<-0 ## variance process error
if(var_alpha){ mod$W[1,1]=exp(parm[2]) }
if(var_beta){ mod$W[2,2]=exp(parm[2]) }
if(var_alpha & var_beta){ mod$W[1,1]=exp(parm[2]);mod$W[2,2]=exp(parm[3])}
return(mod)
}
dlm_out<-suppressWarnings(dlmMLE(y=lnRS,build=build_mod,parm=c(-.1,-.1,-.1),method="L-BFGS-B"))
## maximum likelihood optimization ('L-BFGS-B' or 'Nelder-Mead')
lls<-dlm_out$value ## log-likelihood
dlmMod<-build_mod(dlm_out$par) ## model with variance structure
sigma<-sqrt(dlmMod$V)
outsFilter<-suppressWarnings(dlmFilter(y=lnRS,mod=dlmMod)) ## added Kalman filter
outsSmooth<-dlmSmooth(outsFilter) ## added backward recursive smoothing
alpha_y<-signif(outsSmooth$s[-1,1],4) ## alpha
beta_y<-signif(outsSmooth$s[-1,2],5) ## beta
AICc<-2*lls+2*npara+(2*npara*(npara+1)/(length(data$Rec)-npara-1)) ## AICc
output<-list(results=cbind(data,alpha_y,beta_y),AICc=AICc,sigma=sigma) 
return(output) 
} 
##========================================================================##