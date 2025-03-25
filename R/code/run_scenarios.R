##==================================================================##
##                                                                  ##
##                      Run scenarios for MSE                       ##
#                                                                   ##
##==================================================================##
##=========================================================## packages
pkgs<-c("here","readxl","openxlsx","dplyr","gtools","faraway","gsl", "nlme","matrixcalc","data.table","dlm","progress")
if(length(setdiff(pkgs,rownames(installed.packages())))>0){install.packages(setdiff(pkgs,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkgs,library,character.only=T))
home<-here::here()
##========================================================## functions
fxn<-list.files(paste0(home,"/R/functions/"))
fxn<-fxn[fxn!="parameters.R"]
invisible(sapply(FUN=source,paste0(home,"/R/functions/",fxn)))
##========================================================## scenarios
scen<-data.frame(read_excel(paste0(home,"/R/scenarios.xlsx"))) 
# scen<-scen[c(4,8,12),]
# scen<-scen[c(3,7,11),]
# scen<-scen %>% filter(selectivity=="unselective") %>% filter(factorMSY==1) 
nscen<-dim(scen)[1] ## number of scenarios
niter<-1000 ## iterations per scenario
print(paste0("time estimate: ",round((niter*nscen/20)/60,2)," hours"))
##=====================================================## output lists
output.list<-replicate(nscen,list(replicate(niter,list())))
para.list<-sr_sim.list<-fec.list<-egg.list<-S_msy.list<-data.list<-obs.list<-ret_age.list<-meanSaA.list<-propfemale.list<-select_age.list<-MSY_Goals.list<-impl_errors.list<-output.list
##====================================## loop scenarios and iterations
pb<-progress_bar$new(total=niter*nscen)
# seed.list<-sample(seq(1e6),niter*nscen,replace=F) ## random seed
seed.list<-seq(niter*nscen) ## reproducible seed
start.time<-Sys.time()
for(j in 1:nscen) { 
  for(k in 1:niter) {
    seednum<-seed.list[k] ## same iteration seeds for each scenario
    parms.list<-source(paste0(home,"/R/functions/parameters.R"))$value
    mod.out<-try(run_model(parms.list=parms.list))
    if(class(mod.out)!="try-error") {
      para.list[[j]][[k]]<-mod.out$para ## true alpha & beta parameters
      sr_sim.list[[j]][[k]]<-mod.out$sr_sim ## SR est on simulated data
      fec.list[[j]][[k]]<-mod.out$fec	## change in fecundity over time
      egg.list[[j]][[k]]<-mod.out$egg	## change in eggmass over time
      S_msy.list[[j]][[k]]<-mod.out$S_msy ## S_msy estimates
      data.list[[j]][[k]]<-mod.out$data ## simulated data
      obs.list[[j]][[k]]<-mod.out$obs	## observations
      ret_age.list[[j]][[k]]<-mod.out$ret_by_age	## return age comp
      meanSaA.list[[j]][[k]]<-mod.out$meanSaA	## mean size-at-age
      propfemale.list[[j]][[k]]<-mod.out$propfemale	## prop female
      select_age.list[[j]][[k]]<-mod.out$selectivities_by_age
      MSY_Goals.list[[j]][[k]]<-mod.out$MSY_Goals
      impl_errors.list[[j]][[k]]<-mod.out$impl_errors
    } ## end if statement
    pb$tick()
  } ## end k loop
} ## end j loop
##====================================================## true run time
end.time<-Sys.time()
print(round(end.time-start.time,2)) 
time.save<-gsub(":","-",end.time)
##======================================================## save output
path<-paste0(home,"/R/out/")
dir.create(file.path(path),showWarnings=F)
setwd(file.path(paste0(path)))
save.image(paste0(home,"/R/out/run_",time.save,".Rdata"),compress=T)
saveRDS(parms.list,paste0(home,"/R/out/run_",time.save,"_parms.Rdata"))
write.xlsx(scen,file=paste0(home,"/R/out/run_",time.save,"_scen.xlsx"))
##==================================================================##
##==================================================================##
##==================================================================##