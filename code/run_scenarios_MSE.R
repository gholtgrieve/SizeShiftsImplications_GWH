##========================================================================##
##                                                                        ##
##                           Run scenarios for MSE                        ##
##                                                                        ##
##========================================================================##
# rm(list=ls()) 
##===============================================================## packages
# pkgs<-c("here","readxl","xlsx","dplyr","gtools","faraway","RColorBrewer","gsl","nlme","R2admb","matrixcalc","data.table","R2jags","rstan","sampling","dlm")
pkgs<-c("here","readxl","openxlsx","dplyr","gtools","faraway","RColorBrewer","gsl","nlme","matrixcalc","data.table","dlm")
if(length(setdiff(pkgs,rownames(installed.packages())))>0) { install.packages(setdiff(pkgs,rownames(installed.packages())), dependencies=TRUE) }
invisible(lapply(pkgs,library,character.only=T))
# setup_admb() 
##============================================================## directories
homedir<-here();codedir<-here("code");plotdir<-here("plots");outdir<-here("out");scendir<-here("scenarios")
##==============================================================## functions
setwd(codedir)
source("run_model_MSE.R") ## code to run model for MSE (use large ny)
source("selectivity.R") ## age-selectivity
source("ricker.R") ## generates recruits when beta is constant
source("agecomp.R") ## age composition of recruits
source("reprod_output.R") ## fecundity/egg mass based on size
source("DLMfit.R") ## when using DLMs with time-varying parameters
##==============================================================## scenarios
setwd(homedir)
scen<-data.frame(read_excel("scenarios.xlsx")) ## list of all scenarios
nscen<-dim(scen)[1] ## number of scenarios
##===============================================================## settings
niter<-3 ## iterations per scenario
##=====================================================## estimated run time
## smaller beta means larger  abundances and more time (esp. below 5e-5)
## also depends on ny, escgoalrev, harvmgmt, and some other parameters
time_est<-niter*nscen/20
# if(exists("est_method")) if(est_method=="JAGS") time_est<-time_est*35
if(time_est<60) print(paste0("est.time ",round(time_est,2)," min"))
if(time_est>=60) print(paste0("est.time ",round(time_est/60,2)," hours"))
##===========================================================## output lists
output.list<-replicate(nscen,list(replicate(niter,list())))
para.list<-sr_sim.list<-fec.list<-egg.list<-S_msy.list<-data.list<-obs.list<-output.list
##==========================================## loop scenarios and iterations
# seed.list<-sample(seq(1e6),niter*nscen,replace=F) ## random seed
seed.list<-seq(niter*nscen) ## reproducible seed for comparing test runs
start.time<-Sys.time()
for(j in 1:nscen) { 
for(k in 1:niter) {
seednum<-seed.list[k] ## same iteration seeds for each scenario
print(paste0("scenario=",j," iteration=",k," seed=",seednum))
setwd(codedir);parms.list<-source("run_parameters_MSE.R")$value
mod.out<-try(run_model(parms.list=parms.list))
if(class(mod.out)!="try-error") {
para.list[[j]][[k]]<-mod.out$para ## true alpha and beta parameters
sr_sim.list[[j]][[k]]<-mod.out$sr_sim ## SR estimates on simulated data
fec.list[[j]][[k]]<-mod.out$fec	## change in fecundity over time
egg.list[[j]][[k]]<-mod.out$egg	## change in eggmass over time
S_msy.list[[j]][[k]]<-mod.out$S_msy ## S_msy estimates
data.list[[j]][[k]]<-mod.out$data ## simulated data
obs.list[[j]][[k]]<-mod.out$obs	## observations
} ## end if statement
} ## end k loop
} ## end j loop
##==========================================================## true run time
end.time<-Sys.time()
elapsed<-end.time-start.time
print(round(elapsed,2)) 
time.save<-gsub(":","-",end.time)
##============================================================## save output
setwd(outdir)
save.image(paste0("run_",time.save,".Rdata"),compress=T)
saveRDS(parms.list,paste0("run_",time.save,"_parms.Rdata"))
write.xlsx(scen,file=paste0("run_",time.save,"_scen.xlsx"))
##========================================================================##
##========================================================================##
##========================================================================##