##========================================================================##
##                                                                        ##
##                           Run scenarios for MSE                        ##
##                                                                        ##
##========================================================================##
##===============================================================## packages
pkgs<-c("here","readxl","openxlsx","dplyr","gtools","faraway","RColorBrewer","gsl","nlme","matrixcalc","data.table","dlm","progress")
if(length(setdiff(pkgs,rownames(installed.packages())))>0) { install.packages(setdiff(pkgs,rownames(installed.packages())), dependencies=TRUE) }
invisible(lapply(pkgs,library,character.only=T))
##==============================================================## functions
fxn<-list.files(paste0(here::here(),"/R/functions/"))[fxn!="parameters.R"]
invisible(sapply(FUN=source,paste0(here::here(),"/R/functions/",fxn)))
##==============================================================## scenarios
scen<-data.frame(read_excel("R/scenarios.xlsx")) ## list of all scenarios
nscen<-dim(scen)[1] ## number of scenarios
niter<-100 ## iterations per scenario
print(paste0("time estimate: ",round((niter*nscen/20)/60,2)," hours"))
##===========================================================## output lists
output.list<-replicate(nscen,list(replicate(niter,list())))
para.list<-sr_sim.list<-fec.list<-egg.list<-S_msy.list<-data.list<-obs.list<-output.list
##==========================================## loop scenarios and iterations
pb<-progress_bar$new(total=niter*nscen)
# seed.list<-sample(seq(1e6),niter*nscen,replace=F) ## random seed
seed.list<-seq(niter*nscen) ## reproducible seed for comparing test runs
start.time<-Sys.time()
for(j in 1:nscen) { 
  for(k in 1:niter) {
    seednum<-seed.list[k] ## same iteration seeds for each scenario
    # print(paste0("scenario=",j," iteration=",k," seed=",seednum))
    parms.list<-source("R/functions/parameters.R")$value
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
    pb$tick()
  } ## end k loop
} ## end j loop
##==========================================================## true run time
end.time<-Sys.time()
print(round(end.time-start.time,2)) 
time.save<-gsub(":","-",end.time)
##============================================================## save output
save.image(paste0("R/out/run_",time.save,".Rdata"),compress=T)
saveRDS(parms.list,paste0("R/out/run_",time.save,"_parms.Rdata"))
write.xlsx(scen,file=paste0("R/out/run_",time.save,"_scen.xlsx"))
##========================================================================##
##========================================================================##
##========================================================================##