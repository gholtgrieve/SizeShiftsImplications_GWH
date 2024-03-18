##========================================================================##
##                                                                        ##
##                     Plot simulation model results                      ##
##                                                                        ##
##========================================================================##
# rm(list=ls()) 
##===============================================================## packages
pkgs<-c("here","scales","RColorBrewer","readxl","fishualize","gsl","vioplot", "caroline")
if(length(setdiff(pkgs,rownames(installed.packages())))>0) {install.packages(setdiff(pkgs,rownames(installed.packages())),dependencies=TRUE)}
invisible(lapply(pkgs,library,character.only=T))
'%!in%'<-function(x,y)!('%in%'(x,y))

##============================================================## directories
homedir<-here()
plotdir<-here("plots")
outdir<-here("out")
scendir<-here("scenarios")
##========================================================================##
##===========================================================## load results
##========================================================================##
setwd(outdir)
my.names<-dir();myfile<-my.names[length(my.names)] ## last file saved
load(myfile)
my.names<-dir();myfile<-my.names[length(my.names)] ## last file saved
timestamp<-substr(myfile,5,nchar(myfile)-6)
##-------------------------------------------------------------## parameters
parameters<-readRDS(paste0("run_",timestamp,"_parms.Rdata"))
nyi<-parameters$nyi
ny<-parameters$ny
year_index<-11:(nyi+ny-10) 
nyrec<-length(year_index) ## number of reconstructed years
##----------------------------------------------------## scenarios from file
setwd(scendir)
scenarios<-data.frame(read_excel("scenarios.xlsx"))
nscen<-dim(scenarios)[1]		
##---------------------------------------------------------## scenario names
if(length(unique(scenarios$mgmt))==1){
scenario_names<-paste0(scenarios$trends,"\n",scenarios$selectivity)
}
if(length(unique(scenarios$selectivity))==1){
scenario_names<-paste0(scenarios$trends,"\n",scenarios$mgmt)
}
trend_names<-unique(scenarios$trends)
selectivity_names<-unique(scenarios$selectivity)
mgmt_names<-unique(scenarios$mgmt)
##-----------------------------------------------------## make new directory
setwd(plotdir)
dir.create(file.path(timestamp),showWarnings=F)
setwd(file.path(timestamp))

##=========================================================## results arrays
dnames<-list(paste0("scen=",seq(nscen)),paste0("iter=",seq(niter)))
res_array<-array(dim=c(nscen,niter),dimnames=dnames)
probs<-c(0.05,0.25,0.5,0.75,0.95)
##=================================================================## colors
colors<-c("goldenrod","chocolate2","forestgreen","darkorchid3","mediumblue")
colors<-fish(5,begin=0.1,end=0.9)

##========================================================================##
##==========================## change in overall mean size in simulated data
##========================================================================##
nyref<-10  ## length of early/late periods 
##====================================================## mean size in return
mean_size_trends_mm<-array(NA,dim=c(nscen,niter))
colnames(mean_size_trends_mm)<-paste("iter",seq(niter))
rownames(mean_size_trends_mm)<-paste("scen",seq(nscen))
mean_size_trends_pct<-array(NA,dim=c(nscen,niter))
colnames(mean_size_trends_pct)<-paste("iter",seq(niter))
rownames(mean_size_trends_pct)<-paste("scen",seq(nscen))
##------------------------------------------## loop scenarios and iterations
probs<-c(0.05,0.5,0.95)
for(j in 1:nscen) { 
	for(i in 1:niter) { 
	if(length(unlist(data.list[[j]][i]))>0) {
	sizes_ij<-data.frame(data.list[[j]][i])$meanSret
	earlymean<-mean(sizes_ij[1:nyref],na.rm=T)
	latemean<-mean(sizes_ij[(nyrec-nyref+1):nyrec],na.rm=T)
	mean_size_trends_mm[j,i]<-round((latemean-earlymean),2)
	mean_size_trends_pct[j,i]<-round((latemean-earlymean)*100/earlymean,2)
	}
	}
} 
##--------------------------------------------------## medians and quartiles
mean_size_trends_mm_qs_ret<-apply(mean_size_trends_mm,1,function(x) quantile(x,probs=probs,na.rm=T))
mean_size_trends_pct_qs_ret<-apply(mean_size_trends_pct,1,function(x) quantile(x,probs=probs,na.rm=T))
## currently, scenarios 1, 4, and 7 do not have any demographic trends 
## scenarios 2/3, 5/6, and 8/9 differ in sex trends so are the same here
## for the return fishery selectivity does not matter (but see below)

##================================================## mean size in escapement
mean_size_trends_mm<-array(NA,dim=c(nscen,niter))
colnames(mean_size_trends_mm)<-paste("iter",seq(niter))
rownames(mean_size_trends_mm)<-paste("scen",seq(nscen))
mean_size_trends_pct<-array(NA,dim=c(nscen,niter))
colnames(mean_size_trends_pct)<-paste("iter",seq(niter))
rownames(mean_size_trends_pct)<-paste("scen",seq(nscen))
##------------------------------------------## loop scenarios and iterations
probs<-c(0.05,0.5,0.95)
for(j in 1:nscen) { 
	for(i in 1:niter) { 
	if(length(unlist(data.list[[j]][i]))>0) {
	sizes_ij<-data.frame(data.list[[j]][i])$meanSesc
	earlymean<-mean(sizes_ij[1:nyref],na.rm=T)
	latemean<-mean(sizes_ij[(nyrec-nyref+1):nyrec],na.rm=T)
	mean_size_trends_mm[j,i]<-round((latemean-earlymean),2)
	mean_size_trends_pct[j,i]<-round((latemean-earlymean)*100/earlymean,2)
	}
	}
} 
##--------------------------------------------------## medians and quartiles
mean_size_trends_mm_qs_esc<-apply(mean_size_trends_mm,1,function(x) quantile(x,probs=probs,na.rm=T))
mean_size_trends_pct_qs_esc<-apply(mean_size_trends_pct,1,function(x) quantile(x,probs=probs,na.rm=T))
## now fishery selectivity matters for the change in mean size

##========================================================================##
##========================## resulting trends in mean fecundity and egg mass
##========================================================================##
nyrs<-10 ## number of years for calculating diffference (3,5,10)
if(nyrs==3) usecol<-1;if(nyrs==5) usecol<-2;if(nyrs==10) usecol<-3
fec_trend_mean<-egg_trend_mean<-array(NA,dim=c(nscen,3))
fec_trends<-egg_trends<-array(NA,dim=c(nscen,niter))
##------------------------------------------------------------## array names
colnames(egg_trend_mean)<-colnames(fec_trend_mean)<-c("3yrs","5yrs","10yrs")
colnames(egg_trends)<-colnames(fec_trends)<-paste("iter",seq(niter))
rownames(egg_trend_mean)<-rownames(egg_trend_mean)<-rownames(egg_trends)<-rownames(egg_trends)<-paste("scen",seq(nscen))
##---------------------------------------------------------## loop scenarios
for(j in 1:nscen) { 
fec_trend_scen<-data.frame(t(array(unlist(fec.list[[j]]),dim=c(3,niter))))
egg_trend_scen<-data.frame(t(array(unlist(egg.list[[j]]),dim=c(3,niter))))
fec_trend_mean[j,]<-apply(fec_trend_scen,2,mean)
egg_trend_mean[j,]<-apply(egg_trend_scen,2,mean)
fec_trends[j,]<-fec_trend_scen[,usecol]
egg_trends[j,]<-egg_trend_scen[,usecol]
}

##========================================================================##
##=======================## plot changes in fecundity and egg mass over time
##========================================================================##
ylim<-c(-60,20)
##====================================================## change in fecundity
myplot<-fec_trends 
select_scen<-c(1:nscen) ## select scenarios
n_scen_new<-length(select_scen)
myplot<-myplot[select_scen,]
cols<-colorRampPalette(brewer.pal(9,"Set1"))(n_scen_new)
ylim<-c(min(myplot),max(myplot))
##---------------------------------------------------## boxplots
pdf("percent_change_fecundity.pdf",height=5,width=3+0.25*nscen)
par(mar=c(8.5,3.5,1,1),mgp=c(2,.5,0),tck=-0.02,pch=16,cex.lab=1,lheight=.9)
plot(NA,NA,xaxt="n",xlab="",ylab="Change in mean fecundity per spawner (%)", xlim=c(0.5,.5+n_scen_new),ylim=ylim)
abline(h=0,lty=3)
for(j in 1:n_scen_new) { 
#boxplot(as.numeric(myplot[j,]),at=j,add=T,boxwex=0.9,cex=0.25,lty=1,lwd=0.5,yaxt="n",col=cols[j],outline=F) 
vioplot(as.numeric(myplot[j,]),at=j,add=T,col=cols[j],pchMed=16,colMed=1,wex=0.8,border=NA)
}
axis(side=1,at=seq(n_scen_new),labels=scenario_names,las=2,cex.axis=0.9)
dev.off()
##=====================================================## change in egg mass
myplot<-egg_trends 
select_scen<-c(1:nscen) ## select scenarios
n_scen_new<-length(select_scen)
myplot<-myplot[select_scen,]
cols<-colorRampPalette(brewer.pal(9,"Set1"))(n_scen_new)
ylim<-c(min(myplot),max(myplot))
##---------------------------------------------------## boxplots
pdf("percent_change_eggmass.pdf",height=5,width=3+0.25*nscen)
par(mar=c(8.5,3.5,1,1),mgp=c(2,.5,0),tck=-0.02,pch=16,cex.lab=1,lheight=.9)
plot(NA,NA,xaxt="n",xlab="",ylab="Change in mean egg mass per spawner (%)", xlim=c(0.5,.5+n_scen_new),ylim=ylim)
abline(h=0,lty=3)
for(j in 1:n_scen_new) {
# boxplot(as.numeric(myplot[j,]),at=j,add=T, boxwex=0.9,cex=0.25,lty=1,lwd=0.5,yaxt="n",col=cols[j],outline=F) 
vioplot(as.numeric(myplot[j,]),at=j,add=T,col=cols[j],pchMed=16,colMed=1,wex=0.8,border=NA)
}
axis(side=1,at=seq(n_scen_new),labels=scenario_names,las=2,cex.axis=0.9)
dev.off()

##========================================================================##
##=================================================## equilibrium escapement
##========================================================================##
rep_units<-c("spawners","fecundity","eggmass")
nreps<-length(rep_units)
period_names<-c("first 10 years","all years","last 10 years")
nprds<-length(period_names)
scen_names<-paste("scen",seq(nscen))
iter_names<-paste("iter",seq(niter))

##************************************************************************##
##*******************************## difference between early and late period
##************************************************************************##
myarray<-array(NA,dim=c(nscen,niter))
dimnames(myarray)<-list(scen_names,iter_names)
S_eq_fec_diff<-S_eq_egg_diff<-myarray ## % change late vs early 
##---------------------------------------------------------## loop scenarios
for(j in 1:nscen) { 
S_eq_scen<-array(unlist(S_eq.list[[j]]),dim=c(3,3,niter))
dimnames(S_eq_scen)<-list(rep_units,period_names,iter_names)
S_eq_fec_diff[j,]<-(S_eq_scen[2,3,]-S_eq_scen[2,1,])*100/S_eq_scen[2,1,]
S_eq_egg_diff[j,]<-(S_eq_scen[3,3,]-S_eq_scen[3,1,])*100/S_eq_scen[3,1,]
}

##========================================================================##
##==========================## plot difference in S_msy early vs late period
##========================================================================##
select_scen<-c(1:nscen) ## select scenarios
n_scen_new<-length(select_scen)
cols<-rep(rev(colors[1:4]),n_scen_new/4)
probs<-c(0.05,0.25,0.5,0.75,0.95)
pchs<-rep(21,4)
if(length(unique(scenarios$mgmt))==1) plot_names<-selectivity_names
if(length(unique(scenarios$selectivity))==1) plot_names<-mgmt_names
##---------------------------------------------------------------## fecudity
if(sim_recruits!="eggmass") {
myplot<-S_eq_fec_diff 
myplot<-myplot[select_scen,]
quants<-apply(myplot,1,function(x) quantile(x,probs=probs,na.rm=T))
ylim<-c(min(quants),1.1*max(quants)) #;ylim<-c(-5,85)
pdf("S_msy__early_vs_late_period_fecundity.pdf",height=4,width=2+0.25*nscen)
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2,0.5,0),tck=-0.025,cex.lab=1.1)
plot(NA,NA,xaxt="n",xlab="",ylab="% S_msy change early to late period",xlim=c(0.5,0.5+n_scen_new),ylim=ylim)
abline(h=0,lty=3)
axis(side=1,at=0.5+c(2,6,10),labels=plot_names,cex.axis=0.9)
xx<-seq(n_scen_new);adj<-rep(c(0.6,0.2,-0.2,-0.6),3)
segments(xx+adj,quants[1,],xx+adj,quants[5,],lwd=0.5)
segments(xx+adj,quants[2,],xx+adj,quants[4,],lwd=1.5)
points(xx+adj,quants[3,],pch=pchs,bg=cols,cex=1.2)	
legend("topright",trend_names,pch=pchs,pt.bg=cols,cex=0.8,bty="n")
dev.off()
}
##----------------------------------------------------------------## eggmass
if(sim_recruits!="fecundity") {
myplot<-S_eq_egg_diff 
myplot<-myplot[select_scen,]
quants<-apply(myplot,1,function(x) quantile(x,probs=probs,na.rm=T))
ylim<-c(min(quants),1.1*max(quants)) #;ylim<-c(-10,90)
pdf("S_msy__early_vs_late_period_eggmass.pdf",height=4,width=2+0.25*nscen)
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2,0.5,0),tck=-0.025,cex.lab=1.1)
plot(NA,NA,xaxt="n",xlab="",ylab="% S_msy change early to late period", xlim=c(0.5,0.5+n_scen_new),ylim=ylim)
abline(h=0,lty=3)
axis(side=1,at=0.5+c(2,6,10),labels=plot_names,cex.axis=0.9)
xx<-seq(n_scen_new);adj<-rep(c(0.6,0.2,-0.2,-0.6),3)
segments(xx+adj,quants[1,],xx+adj,quants[5,],lwd=0.5)
segments(xx+adj,quants[2,],xx+adj,quants[4,],lwd=1.5)
points(xx+adj,quants[3,],pch=pchs,bg=cols,cex=1.2)	
legend("topright",trend_names,pch=pchs,pt.bg=cols,cex=0.8,bty="n")
dev.off()
}

##************************************************************************##
## difference spawner model Smsy based on Ricker vs equilibirum calculations
##************************************************************************##
## sr_rep_out.list only for spawner abundance model to allow this comparison
myarray<-array(NA,dim=c(nscen,niter))
dimnames(myarray)<-list(scen_names,iter_names)
S_eq_msy_diff<-myarray
for(j in 1:nscen) {
S_eq_scen<-array(unlist(S_eq.list[[j]]),dim=c(nreps,nprds,niter))
dimnames(S_eq_scen)<-list(rep_units,period_names,iter_names)
S_eq_spawnrs<-as.numeric(S_eq_scen[1,2,]) ## spawner model all years
S_msy_spawnrs<-round(array(unlist(sr_rep_out.list[[j]]),dim=c(6,niter))[5,])
S_eq_msy_diff[j,]<-(S_eq_spawnrs-S_msy_spawnrs)*100/S_msy_spawnrs
} ## end j loop
S_eq_msy_diff_by_scen<-apply(S_eq_msy_diff,1,mean) ## mean diff by scenario
if(sum(S_eq_msy_diff)==0) { print("Same estimates for S_msy and S_eq") } else { print("Different estimates for S_msy and S_eq") }
## confirmation that the equilibrium calculations yield the same result for  spawner abundance model when using the estimated Ricker parameters directly
## difference when using DLM as estimation method

##************************************************************************##
##***************## compare spawner model Smsy to Smsy from time-varying DLM
##************************************************************************##
if(explore_DLM){
##------------------------------------------------------## compute difference
nydlm<-10 ## number of years averaged (e.g. first/last 10)
S_msy_diff_early<-S_msy_diff_all<-S_msy_diff_late<-array(NA,dim=c(nscen,niter),dimnames=list(scen_names,iter_names))
# S_msys_dlm<-array(NA,dim=c(nscen,niter,ny-nyi),dimnames=list(scen_names,iter_names,seq(ny-nyi)))
for(j in 1:nscen) { 
  # print(paste("j=",j))
  for(k in 1:niter) { 
    # print(k)
    S_msy_dlm<-NA
    S_msy_dlm<-data.frame(DLM.list[[j]][[k]])$S_msy_dlm
    if(is.null(S_msy_dlm)) { next } else {
    S_msy<-round(data.frame(sr_rep_out.list[[j]][[k]])$S_msy)
    S_msy_diff_all[j,k]<-100*(mean(S_msy_dlm)-S_msy)/S_msy
    S_msy_diff_early[j,k]<-100*(mean(mean(S_msy_dlm[1:nydlm]))-S_msy)/S_msy
    S_msy_diff_late[j,k]<-100*(mean(mean(S_msy_dlm[(ny-nyi-nydlm+1):(ny-nyi)]))-S_msy)/S_msy
    }## end if statement
  } ## end k-loop
} ## end j-loop
all_list<-list(S_msy_diff_early,S_msy_diff_all,S_msy_diff_late)

##========================## plot alpha trajectories by scenario/selectivity
pdf("DLM_alpha_trends.pdf",height=6,width=6)
layout(matrix(1:nscen,nrow=3,ncol=3,byrow=T)) 
par(mar=c(1,1,1,1),oma=c(3,3,0,0),mgp=c(2,0.25,0),tck=-0.025,cex.axis=0.8)
for(j in 1:nscen) { 
  plot(NA,NA,xlab="",ylab="",xlim=c(0,ny-nyi),ylim=c(0,20))
  text<-paste0(scenarios$selectivity[j],"\n",scenarios$trend[j])
  mtext(text,side=3,cex=0.6,line=-2,font=2)
  for(k in 1:niter) { 
    S_msy_dlm<-NA
    alpha_dlm_jk<-data.frame(DLM.list[[j]][[k]])$alpha_dlm
    lines(alpha_dlm_jk,lwd=0.1)
  } ## end k-loop
} ## end j-loop
mtext("Year",side=1,cex=1,line=1,font=1,outer=T)
mtext("alpha",side=2,cex=1,line=1,font=1,outer=T)
dev.off()

##================================## plot difference by scenario/selectivity
pdf("DLM_Smsy_difference.pdf",height=6,width=3)
layout(matrix(1:np,ncol=1,byrow=T)) 
par(mar=c(0.5,5,0,0),oma=c(6,0,1,1),mgp=c(2,0.5,0),tck=-0.025)
##---------------------------------------------## loop periods
for(p in 1:np) {
if(p==1) select_scen<-c(1:3) 
if(p==2) select_scen<-c(4:6) 
if(p==3) select_scen<-c(7:9) 
myplot1<-all_list[[1]][select_scen,]
myplot2<-all_list[[2]][select_scen,]
myplot3<-all_list[[3]][select_scen,]
n_scen_new<-length(select_scen)
ps<-c(0.05,0.25,0.5,0.75,0.95) ## median and 50% and 90% quantiles
myquants1<-apply(myplot1,1,function(x) quantile(x,probs=ps,na.rm=T))
myquants2<-apply(myplot2,1,function(x) quantile(x,probs=ps,na.rm=T))
myquants3<-apply(myplot3,1,function(x) quantile(x,probs=ps,na.rm=T))
quants<-c(myquants1,myquants2,myquants3)
##---------------------------------------------## color/names/axes
my_names<-scenarios$trend[select_scen]
my_names<-gsub(" ","\n",my_names)
panel_name<-unique(scenarios$selectivity[select_scen])
cols<-c("goldenrod1","chocolate1","firebrick1")
##---------------------------------------------## plot
plot(NA,NA,xaxt="n",xlab="",ylab="% difference in S_msy\ncompared to spawner model",xlim=c(0.5,0.5+n_scen_new),ylim=c(min(quants),max(quants)))
abline(h=0,lty=3)
if(p==np) axis(side=1,at=seq(n_scen_new),labels=my_names,cex.axis=1,las=2)
legend("topleft",rev(period_names),pch=c(25,21,24),cex=0.8,bty="n")
mtext(paste0(panel_name," "),side=3,cex=0.6,line=-1.2,font=2)
xx<-seq(n_scen_new)
##---------------------------------------------------## first 10 yrs
adj<--0.25
segments(xx+adj,myquants1[1,],xx+adj,myquants1[5,],col="darkgray")
segments(xx+adj,myquants1[2,],xx+adj,myquants1[4,],lwd=2)
points(xx+adj,myquants1[3,],pch=24,bg=cols,cex=1.5)
##-------------------------------------------------------## all yrs
adj<-0
segments(xx+adj,myquants2[1,],xx+adj,myquants2[5,],col="darkgray")
segments(xx+adj,myquants2[2,],xx+adj,myquants2[4,],lwd=2)
points(xx+adj,myquants2[3,],pch=21,bg=cols,cex=1.5)
##---------------------------------------------------## last 10 yrs
adj<-0.25
segments(xx+adj,myquants3[1,],xx+adj,myquants3[5,],col="darkgray")
segments(xx+adj,myquants3[2,],xx+adj,myquants3[4,],lwd=2)
points(xx+adj,myquants3[3,],pch=25,bg=cols,cex=1.5)
##-----------------------------------------------## end loop over p
}
dev.off()
} ## end if(explore_DLM) statement 

##************************************************************************##
## plot difference fecundity or egg mass models compared to spawner model ##
##************************************************************************##
## difference in S_msy compared to spawner model given trends for all years
myarray<-array(NA,dim=c(nscen,niter))
dimnames(myarray)<-list(scen_names,iter_names)
##==================================================## equilibrium escapement
S_eq_fec_diff_first<-S_eq_fec_diff_all<-S_eq_fec_diff_last<-myarray
S_eq_egg_diff_first<-S_eq_egg_diff_all<-S_eq_egg_diff_last<-myarray
##---------------------------------------## loop scenarios
for(j in 1:nscen) { 
S_eq_scen<-array(unlist(S_eq.list[[j]]),dim=c(nreps,nprds,niter))
dimnames(S_eq_scen)<-list(rep_units,period_names,iter_names)
##---------------------------------------## difference fecundity vs spawners
S_eq_fec_diff_first[j,]<-(S_eq_scen[2,1,]-S_eq_scen[1,2,])*100/S_eq_scen[1,2,]
S_eq_fec_diff_all[j,]<-(S_eq_scen[2,2,]-S_eq_scen[1,2,])*100/S_eq_scen[1,2,]
S_eq_fec_diff_last[j,]<-(S_eq_scen[2,3,]-S_eq_scen[1,2,])*100/S_eq_scen[1,2,]
##----------------------------------------## difference egg mass vs spawners
S_eq_egg_diff_first[j,]<-(S_eq_scen[3,1,]-S_eq_scen[1,2,])*100/S_eq_scen[1,2,]
S_eq_egg_diff_all[j,]<-(S_eq_scen[3,2,]-S_eq_scen[1,2,])*100/S_eq_scen[1,2,]
S_eq_egg_diff_last[j,]<-(S_eq_scen[3,3,]-S_eq_scen[1,2,])*100/S_eq_scen[1,2,]
} ## end j loop
##=====================================================## equilibrium harvest
H_eq_fec_diff_first<-H_eq_fec_diff_all<-H_eq_fec_diff_last<-myarray
H_eq_egg_diff_first<-H_eq_egg_diff_all<-H_eq_egg_diff_last<-myarray
##---------------------------------------## loop scenarios
for(j in 1:nscen) { 
H_eq_scen<-array(unlist(H_eq.list[[j]]),dim=c(nreps,nprds,niter))
dimnames(H_eq_scen)<-list(rep_units,period_names,iter_names)
##---------------------------------------## difference fecundity vs spawners
H_eq_fec_diff_first[j,]<-(H_eq_scen[2,1,]-H_eq_scen[1,2,])*100/H_eq_scen[1,2,]
H_eq_fec_diff_all[j,]<-(H_eq_scen[2,2,]-H_eq_scen[1,2,])*100/H_eq_scen[1,2,]
H_eq_fec_diff_last[j,]<-(H_eq_scen[2,3,]-H_eq_scen[1,2,])*100/H_eq_scen[1,2,]
##----------------------------------------## difference egg mass vs spawners
H_eq_egg_diff_first[j,]<-(H_eq_scen[3,1,]-H_eq_scen[1,2,])*100/H_eq_scen[1,2,]
H_eq_egg_diff_all[j,]<-(H_eq_scen[3,2,]-H_eq_scen[1,2,])*100/H_eq_scen[1,2,]
H_eq_egg_diff_last[j,]<-(H_eq_scen[3,3,]-H_eq_scen[1,2,])*100/H_eq_scen[1,2,]
} ## end j loop

##========================================================================##
##==========## plot difference in S_msy by period for each reproductive unit
##========================================================================##
if(sim_recruits=="fecundity") rep_unit_names<-"fecundity"
if(sim_recruits=="eggmass") rep_unit_names<-"eggmass"
metrics<-c("S_msy","H_msy")
##-----------------------------------------------## loop metrics (Smsy/Hmsy)
for(m in 1:length(metrics)) {
  metric<-metrics[m]
  ##----------------------------------------------## loop reproductive units
  for(r in 1:length(rep_unit_names)) {
    use<-rep_unit_names[r]
    np<-3 ## number of plots
    if(metric=="S_msy") {
      if(use=="fecundity"){ 
        myplot_1<-S_eq_fec_diff_first
        myplot_2<-S_eq_fec_diff_all
        myplot_3<-S_eq_fec_diff_last 
      }
      if(use=="eggmass"){ 
        myplot_1<-S_eq_egg_diff_first
        myplot_2<-S_eq_egg_diff_all
        myplot_3<-S_eq_egg_diff_last 
      }
    }
    if(metric=="H_msy") {
      if(use=="fecundity"){ 
        myplot_1<-H_eq_fec_diff_first
        myplot_2<-H_eq_fec_diff_all
        myplot_3<-H_eq_fec_diff_last 
      }
      if(use=="eggmass"){ 
        myplot_1<-H_eq_egg_diff_first
        myplot_2<-H_eq_egg_diff_all
        myplot_3<-H_eq_egg_diff_last 
      }	
    }
    ##---------------------------------------------------## start plot
    for(d in 1:2) { ## orientation (vertical|horizontal)
      if(np==1) { 
        pdf_name<-paste0(metric,"_diff_",use,"_vs_spawner_model_",d,".pdf")
        pdf(pdf_name,height=4.5,width=2+0.5*nscen) 
      }
      if(np!=1) { 
        if(d==1) {
          pdf_name<-paste0(metric,"_diff_",use,"_vs_spawner_model_",d,".pdf")
          pdf(pdf_name,height=6,width=3)
          layout(matrix(1:np,ncol=1,byrow=T)) 
          par(mar=c(0.5,5,0,0),oma=c(7,0,1,1),mgp=c(2,0.5,0),tck=-0.025)
        }
        if(d==2) {
          pdf_name<-paste0(metric,"_diff_",use,"_vs_spawner_model_",d,".pdf")
          pdf(pdf_name,height=3,width=7)
          layout(matrix(1:np,ncol=np,byrow=T)) 
          par(mar=c(0.5,0.5,0,0),oma=c(7,5,1,1),mgp=c(2,0.5,0),tck=-0.025)
        }
      }
      for(p in 1:np) {
        if(p==1) select_scen<-c(1:4) ## select scenarios
        if(p==2) select_scen<-c(5:8) ## select scenarios
        if(p==3) select_scen<-c(9:12) ## select scenarios
        if(p==1 & np==1) select_scen<-c(1:nscen) ## select scenarios
        myplot1<-myplot_1[select_scen,]
        myplot2<-myplot_2[select_scen,]
        myplot3<-myplot_3[select_scen,]
        n_scen_new<-length(select_scen)
        ps<-c(0.05,0.25,0.5,0.75,0.95) ## median and 50% and 90% quantiles
        myquants1<-apply(myplot1,1,function(x) quantile(x,probs=ps,na.rm=T))
        myquants2<-apply(myplot2,1,function(x) quantile(x,probs=ps,na.rm=T))
        myquants3<-apply(myplot3,1,function(x) quantile(x,probs=ps,na.rm=T))
        quants<-c(myquants1, myquants2, myquants3)
        ##---------------------------------------------## color/names/axes
        my_names<-scenarios$trend[select_scen]
        my_names<-gsub(" ","\n",my_names)
        if(length(unique(scenarios$mgmt))==1) {
          panel_name<-panel_name<-unique(scenarios$selectivity[select_scen])
        }
        if(length(unique(scenarios$selectivity))==1) {
          panel_name<-panel_name<-unique(scenarios$mgmt[select_scen]) 
        }
        cols<-colors[c(4,3,2,1)]
        #cols<-c("goldenrod1","chocolate1","firebrick1")
        if(m==1) ylim<-c(-40,100); if(m==2) ylim<-c(-50,50)
        ##===========================================================## plot
        if(m==1) y_lab<-"% difference S_msy from\nspawner model all years"
        if(m==2) y_lab<-"% difference H_eq from\nspawner model all years"
        if(d==1){ 
          xlim<-c(0.5,0.5+n_scen_new)
          plot(NA,NA,xaxt="n",xlim=xlim,ylim=ylim,xlab="",ylab=ylab)
          if(p==np) axis(side=1,at=seq(n_scen_new),labels=my_names,las=2)
        }
        if(d==2) {
          xlim<-c(0.5,0.5+n_scen_new)
          plot(NA,NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=xlim,ylim=ylim) 
          if(p==1) {
            axis(side=2,at=seq(-100,100,20),labels=T,cex.axis=1,las=2)
            mtext(y_lab,side=2,cex=0.9,line=2.2)
          }
          axis(side=1,at=seq(n_scen_new),labels=my_names,cex.axis=1,las=2)
        }
        legend("topleft",rev(period_names),pch=c(25,21,24),cex=0.8,bty="n")
        mtext(paste0(panel_name," "),side=3,cex=0.6,line=-1.2,font=2)
        abline(h=0,lty=3)
        probs<-c(0.05,0.25,0.5,0.75,0.95) ## median and 50% and 95% quantiles
        xx<-seq(n_scen_new)
        ##---------------------------------------------------## first 10 yrs
        adj<--0.25
        segments(xx+adj,myquants1[1,],xx+adj,myquants1[5,],col="darkgray")
        segments(xx+adj,myquants1[2,],xx+adj,myquants1[4,],lwd=2)
        points(xx+adj,myquants1[3,],pch=24,bg=cols,cex=1.5)
        ##-------------------------------------------------------## all yrs
        adj<-0
        segments(xx+adj,myquants2[1,],xx+adj,myquants2[5,],col="darkgray")
        segments(xx+adj,myquants2[2,],xx+adj,myquants2[4,],lwd=2)
        points(xx+adj,myquants2[3,],pch=21,bg=cols,cex=1.5)
        ##---------------------------------------------------## last 10 yrs
        adj<-0.25
        segments(xx+adj,myquants3[1,],xx+adj,myquants3[5,],col="darkgray")
        segments(xx+adj,myquants3[2,],xx+adj,myquants3[4,],lwd=2)
        points(xx+adj,myquants3[3,],pch=25,bg=cols,cex=1.5)
        ##-----------------------------------------------## end loop over p
      }
      ##---------------------------------------------------## save plot
      dev.off()
    } ## end loop over direction/orientation
  } ## end loop over reproductive units
} ## end loop over metrics

##************************************************************************##
##******************## S_msy of spawner model compared to alternative models
##************************************************************************##
myarray<-array(NA,dim=c(nscen,niter))
dimnames(myarray)<-list(scen_names,iter_names)
S_eq_diff_spawner_vs_fecundity<-S_eq_diff_spawner_vs_eggmass<-myarray
S_eq_bias_spawner_vs_fecundity<-S_eq_bias_spawner_vs_eggmass<-myarray
##---------------------------------------------------------## loop scenarios
for(j in 1:nscen) { 
S_eq_scen<-array(unlist(S_eq.list[[j]]),dim=c(3,3,niter))
dimnames(S_eq_scen)<-list(rep_units,period_names,iter_names)
##---------------------------------------## difference fecundity vs spawners
S_eq_diff_spawner_vs_fecundity[j,]<-(S_eq_scen[1,2,]-S_eq_scen[2,2,])/S_eq_scen[2,2,]
S_eq_bias_spawner_vs_fecundity[j,]<-(S_eq_scen[1,2,]-S_eq_scen[2,2,])*100/S_eq_scen[2,2,] ## same but as percent bias
##----------------------------------------## difference egg mass vs spawners
S_eq_diff_spawner_vs_eggmass[j,]<-(S_eq_scen[1,2,]-S_eq_scen[3,2,])/S_eq_scen[3,2,]
S_eq_bias_spawner_vs_eggmass[j,]<-(S_eq_scen[1,2,]-S_eq_scen[3,2,])*100/S_eq_scen[3,2,] ## same but as percent bias
}

##========================================================================##
##====## probability of underestimating S_msy by x% when using spawner model
##========================================================================##
thresholds<-c(10,25) ## probability across iterations for each scenario
nts<-length(thresholds)
select_scen<-c(1:nscen) ## select scenarios
n_scen_new<-length(select_scen)
#cols<-rep(c("goldenrod1","chocolate1","firebrick1"),3)
cols<-rep(rev(colors[1:4]),n_scen_new/4)
if(length(unique(scenarios$mgmt))==1) plot_names<-selectivity_names
if(length(unique(scenarios$selectivity))==1) plot_names<-mgmt_names

##========================## prob to under-estimate S_msy using spawner model
for(t in 1:nts){
threshold<-thresholds[t]/100
##-----------------------------------------------## compared to eggmass model
if(sim_recruits!="fecundity") {
pdf(paste0("S_msy_spawner_model_prob_underest_by_",thresholds[t],"+pct.pdf"), height=4,width=2+0.25*nscen)
par(mar=c(3.5,4.5,1,1),mgp=c(2,0.5,0),tck=-0.025,cex.lab=1,cex.axis=0.9)
prob_underest_smsy<-apply(S_eq_diff_spawner_vs_eggmass,1, function(x) sum(x<=-threshold,na.rm=T)/length(x[!is.na(x)]))
prob_plot<-prob_underest_smsy[select_scen]
x<-seq(n_scen_new);xlim<-c(0.5,.5+n_scen_new)
ylim=c(0,1.1*max(prob_plot))
ylab<-paste0("Probability to underestimate S_msy by >",thresholds[t],"%", "\nusing spawner compared to eggmass model")
plot(NA,NA,xaxt="n",xlab="",ylab=ylab,xlim=xlim,ylim=ylim)
axis(side=1,at=0.5+c(2,6,10),labels=plot_names)
xx<-seq(n_scen_new);adj<-rep(c(0.6,0.2,-0.2,-0.6),3)
points(xx+adj,prob_plot,pch=21,bg=cols,cex=1.5)	
legend("bottomleft",trend_names,pch=21,pt.bg=cols,cex=0.8,bty="n")
dev.off()
}
##--------------------------------------------------## fecundity
if(sim_recruits!="eggmass") {
pdf(paste0("S_msy_spawner_model_prob_underest_by_",thresholds[t],"+pct.pdf"),height=4,width=2+0.25*nscen)
par(mar=c(3.5,4.5,1,1),mgp=c(2,0.5,0),tck=-0.025,cex.lab=1,cex.axis=0.9)
prob_underest_smsy<-apply(S_eq_diff_spawner_vs_fecundity,1, function(x) sum(x<=-threshold,na.rm=T)/length(x[!is.na(x)]))
prob_plot<-prob_underest_smsy[select_scen]
x<-seq(n_scen_new);xlim<-c(0.5,.5+n_scen_new)
ylim=c(0,1.1*max(prob_plot))
ylab<-paste0("Probability to underestimate S_msy by >",thresholds[t],"%", "\nusing spawner compared to fecundity model")
plot(NA,NA,xaxt="n",xlab="",ylab=ylab,xlim=xlim,ylim=ylim)
axis(side=1,at=0.5+c(2,6,10),labels=plot_names)
xx<-seq(n_scen_new);adj<-rep(c(0.6,0.2,-0.2,-0.6),3)
points(xx+adj,prob_plot,pch=21,bg=cols,cex=1.5)	
legend("bottomleft",trend_names,pch=21,pt.bg=cols,cex=0.8,bty="n")
dev.off()
}
##--------------------------------------------## end loop over thresholds
} 

##========================================================================##
##=## S_msy poterior uncertainty in spawner model (Bayesian estimation only)
##========================================================================##
if(est_method %in% c("JAGS","STAN")) {
myarray<-array(NA,dim=c(nscen,5))
quantiles<-c("5%","25%","50%","75%","95%")
dimnames(myarray)<-list(scen_names,quantiles)
S_msy_post_quant<-myarray ## S_msy posterior quantiles of spawner model
##---------------------------------------------------------## loop scenarios
for(j in 1:nscen) { 
S_msy_scen<-array(unlist(S_msy_posterior.list[[j]]),dim=c(5,niter))
dimnames(S_msy_scen)<-list(quantiles,iter_names)
S_msy_post_quant[j,]<-rowMeans(S_msy_scen)
}
## average % uncertainty in S_msy across iterations for the spawner model
##=================================================## plot for each scenario
trend_names<-c("no trend","age-size trend","age-size-sex trend")
pdf("S_msy_relative_uncertainty_in_spawner_model.pdf",height=5,width=4)
par(mar=c(9,5,1,1),mgp=c(2,0.5,0),tck=-0.025,cex.lab=1.2)
select_scen<-c(1:3) ## select scenarios (same for all selectivities)
n_scen_new<-length(select_scen)
xx<-seq(n_scen_new)
plot_quants<-S_msy_post_quant[select_scen,]
plot(NA,NA,xlab="",ylab="Uncertainty in S_msy\nin spawner model (%)", xaxt="n",xlim=c(0.5,3.5),ylim=c(min(plot_quants),max(plot_quants)))
axis(side=1,at=seq(n_scen_new),labels=trend_names,cex.axis=1.2,las=2)
abline(h=0,lty=3)
xx<-seq(n_scen_new)
segments(xx,plot_quants[,1],xx,plot_quants[,5],lwd=1,col="darkgray")
segments(xx,plot_quants[,2],xx,plot_quants[,4],lwd=2)
points(xx,plot_quants[,3],pch=16,bg=cols,cex=1.2)
dev.off()
}

##========================================================================##
##=======## Fecundtiy / egg mass S_msy as percentile of spawner model (JAGS)
##========================================================================##
if(est_method %in% c("JAGS","STAN")){
S_msy_percentiles_fec<-S_msy_percentiles_eggm<-array(NA,dim=c(nscen,niter))
##=========================================================## loop scenarios
for(j in 1:nscen) { 
S_msy_percentile_scen<-array(unlist(percentile.list[[j]]),dim=c(2,niter))
dimnames(S_msy_percentile_scen)<-list(rep_units[2:3],iter_names)
S_msy_percentiles_fec[j,]<-S_msy_percentile_scen[1,]
S_msy_percentiles_eggm[j,]<-S_msy_percentile_scen[2,]
}
##===============================================## quantiles of percentiles
probs<-c(0.05,0.25,0.5,0.75,0.95)
S_msy_percentiles_fec_quants<-apply(S_msy_percentiles_fec,1,function(x) quantile(x,prob=probs))
S_msy_percentiles_eggm_quants<-apply(S_msy_percentiles_eggm,1,function(x) quantile(x,prob=probs))
##=========================================================## fecundity plot
pdf("S_msy_fecundity_model_as_percentile_of_S_msy_in_spawner_model.pdf", height=5,width=4)
par(mar=c(8,5,2,1),mgp=c(2,0.5,0),tck=-0.025,cex.lab=1.2)
select_scen<-c(1:nscen) ## select scenarios (same for all selectivities)
n_scen_new<-length(select_scen)
myplot<-S_msy_percentiles_fec_quants[,select_scen]
##-------------------------------------------------------## main
plot(NA,NA,xlab="",ylab="Median S_msy in fecundity model as\npercentile of spawner model S_msy",xaxt="n",xlim=c(0.5,n_scen_new+0.5), ylim=c(min(c(50,myplot)),max(myplot)))
axis(side=1,at=seq(n_scen_new),labels=rep(trend_names,3),cex.axis=1,las=2)
abline(h=50,lty=3)
abline(v=c(3.5,6.5),lty=1,lwd=0.5)
lnames<-c("8.5 inch gillnet    unselective    6 inch gillnet ")
mtext(lnames,side=3,cex=0.8,line=0.5,font=2)
##---------------------------------------------------------## segemnts
xx<-seq(n_scen_new)
segments(xx,myplot[1,],xx,myplot[5,],lwd=1,col="darkgray")
segments(xx,myplot[2,],xx,myplot[4,],lwd=2)
points(xx,myplot[3,],pch=16,bg=cols,cex=1.5)
##---------------------------------------------------------## or boxplot
# for(j in 1:n_scen_new) { boxplot(as.numeric(myplot[,j]),at=j,add=T, boxwex=0.9,cex=0.25,lty=1,lwd=0.5,yaxt="n",col=cols[j],outline=F) }
##---------------------------------------------------------## save plot
dev.off()

##=========================================================## egg mass plot
pdf("S_msy_eggmass_model_as_percentile_of_S_msy_in_spawner_model.pdf", height=5,width=4)
par(mar=c(8,5,2,1),mgp=c(2,0.5,0),tck=-0.025,cex.lab=1.2)
select_scen<-c(1:nscen) ## select scenarios (same for all selectivities)
n_scen_new<-length(select_scen)
myplot<-S_msy_percentiles_eggm_quants[,select_scen]
##-------------------------------------------------------## main
plot(NA,NA,xlab="",ylab="Median S_msy in egg mass model as\npercentile of spawner model S_msy",xaxt="n",xlim=c(0.5,n_scen_new+0.5), ylim=c(min(c(50,myplot)),max(myplot)))
axis(side=1,at=seq(n_scen_new),labels=rep(trend_names,3),cex.axis=1,las=2)
abline(h=50,lty=3)
abline(v=c(3.5,6.5),lty=1,lwd=0.5)
lnames<-c("8.5 inch gillnet    unselective    6 inch gillnet ")
mtext(lnames,side=3,cex=0.8,line=0.5,font=2)
##---------------------------------------------------------## segemnts
xx<-seq(n_scen_new)
segments(xx,myplot[1,],xx,myplot[5,],lwd=1,col="darkgray")
segments(xx,myplot[2,],xx,myplot[4,],lwd=2)
points(xx,myplot[3,],pch=16,bg=cols,cex=1.5)
##---------------------------------------------------------## or boxplot
# for(j in 1:n_scen_new) { boxplot(as.numeric(myplot[,j]),at=j,add=T, boxwex=0.9,cex=0.25,lty=1,lwd=0.5,yaxt="n",col=cols[j],outline=F) }
##---------------------------------------------------------## save plot
dev.off()
##================================================## end method if statement
} 

##========================================================================##
##============================================## simulated vs observed S_msy
##========================================================================##
## S_msy calculated using Scheuerell 2016 formula
smsys_obs<-smsys_sim<-res_array
alphas_obs<-alphas_sim<-res_array
betas_obs<-betas_sim<-res_array
##--------------------------## loop cases for estimated using simulated data
for(i in 1:nscen) { 
	for(j in 1:niter) { 
		if(length(sr_sim.list[[i]][[j]])==0) j=j+1 else {
		alphas_sim[i,j]<-alpha<-sr_sim.list[[i]][[j]]$alpha
		betas_sim[i,j]<-beta<-sr_sim.list[[i]][[j]]$beta
		smsys_sim[i,j]<-(1-lambert_W0(exp(1-log(alpha))))/beta
		}
	} 
} 
##---------------------------## loop cases for estimated using observed data
for(i in 1:nscen) { 
	for(j in 1:niter) {
		if(length(sr_obs.list[[i]][[j]])==0) j=j+1 else {
		alphas_obs[i,j]<-alpha<-sr_obs.list[[i]][[j]]$alpha
		betas_obs[i,j]<-beta<-sr_obs.list[[i]][[j]]$beta
		smsys_obs[i,j]<-(1-lambert_W0(exp(1-log(alpha))))/beta
		}
	}
} 
##------------------------------------## S_msy ratio of observed / simulated
smsys_ratio<-smsys_obs/smsys_sim
quants_smsys_ratio<-apply(smsys_ratio,1, function(x) quantile(x,prob=probs,na.rm=T))
##------------------------------------------------## relative error in S_msy
RE_smsy<-abs((smsys_obs-smsys_sim)/smsys_sim)
RE_smsys_quants<-apply(RE_smsy,1, function(x) quantile(x,prob=probs,na.rm=T))
##----------------------------------------------------------## bias in S_msy
bias_smsy<-(smsys_obs-smsys_sim)/smsys_sim
bias_smsys_quants<-apply(bias_smsy,1, function(x) quantile(x,prob=probs,na.rm=T))

##=======================================================================##
##=============================## plot ratio of observed vs simulated S_smy
##=======================================================================##
myplot<-smsys_ratio 
select_scen<-c(1:nscen) ## select scenarios
n_scen_new<-length(select_scen)
myplot<-myplot[select_scen,]
quants<-apply(myplot,1,function(x) quantile(x,probs=c(0.01,0.99),na.rm=T))
cols<-colorRampPalette(brewer.pal(9,"Set1"))(n_scen_new)
##---------------------------------------------------## boxplots
pdf("observed_vs_simulated_smsy_ratio.pdf",width=5.5,height=2.5+0.25*nscen)
par(mar=c(3.5,8.5,1,1),mgp=c(2,.5,0),tck=-0.02,pch=16,cex.lab=1.1,lheight=.9)
xlim<-c(-1.1*max(log(abs(quants))),1.1*max(log(abs(quants))))
plot(NA,NA,yaxt="n",ylab="",xlab="log(Smsy ratio)",ylim=c(0.5,.5+n_scen_new),xlim=xlim)
abline(v=0,lty=3)
for(j in 1:n_scen_new) { 
#boxplot(log(as.numeric(myplot[j,])),at=j,add=T,lty=1,boxwex=0.9,cex=0.25,lwd=0.5,yaxt="n",col=cols[j],outline=F,horizontal=T) 
vioplot(log(as.numeric(myplot[j,])),at=j,add=T,col=cols[j],pchMed=16,colMed=1,wex=0.8,horizontal=T,border=NA)
}
axis(side=2,at=seq(n_scen_new),labels=scenario_names,las=2,cex.axis=0.9)
dev.off()

##========================================================================##
##=====## autocorrelation and residuals of estimates - simulated vs observed
##========================================================================##
sigmas_obs<-sigmas_sim<-res_array
rhos_obs<-rhos_sim<-res_array
ntot<-niter*nscen
##==========================## loop cases for estimated using simulated data
for(i in 1:nscen) { 
	for(j in 1:niter) {
		if(length(sr_sim.list[[i]][[j]])==0) j=j+1 else {
		rhos_sim[i,j]<-data.frame(sr_sim.list[[i]][[j]])$rho
		sigmas_sim[i,j]<-data.frame(sr_sim.list[[i]][[j]])$sigma
		}
	} 
} 
##===========================## loop cases for estimated using observed data
for(i in 1:nscen) { 
	for(j in 1:niter) {
		if(length(sr_obs.list[[i]][[j]])==0) j=j+1 else {
		rhos_obs[i,j]<-data.frame(sr_obs.list[[i]][[j]])$rho
		sigmas_obs[i,j]<-data.frame(sr_obs.list[[i]][[j]])$sigma
		}
	}
} 
##--------------------------------------------------------## residual error
## ratio of model residual errors based on simulated versus observed
sigma_sim_quants<-apply(sigmas_sim,1, function(x) quantile(x,prob=probs,na.rm=T))
sigma_obs_quants<-apply(sigmas_obs,1, function(x) quantile(x,prob=probs,na.rm=T))
sigma_ratio<-sigmas_obs/sigmas_sim 
quants_sigma_ratio<-apply(sigma_ratio,1, function(x) quantile(x,prob=probs,na.rm=T))
##--------------------------------------------------------## autocorrelation
## ratio of estimated autocorrelation based on simulated versus observed
rho_sims_quants<-apply(rhos_sim,1, function(x) quantile(x,prob=probs,na.rm=T))
rho_obs_quants<-apply(rhos_obs,1, function(x) quantile(x,prob=probs,na.rm=T))
rho_ratio<-rhos_obs/rhos_sim 
quants_rho_ratio<-apply(rho_ratio,1, function(x) quantile(x,prob=probs,na.rm=T))

##========================================================================##
##==============================================## plot residual error ratio
##========================================================================##
myplot<-sigma_ratio 
drop<-length(myplot[myplot<=0])
myplot[myplot<=0]<-0.1
pdrop<-drop/ntot ## proportion negative not included in log-ratio
select_scen<-c(1:nscen) ## select scenarios
n_scen_new<-length(select_scen)
myplot<-myplot[select_scen,]
quants<-apply(myplot,1,function(x) quantile(x,probs=c(0.01,0.99),na.rm=T))
cols<-colorRampPalette(brewer.pal(9,"Set1"))(n_scen_new)
##---------------------------------------------------## boxplots
pdf("observed_vs_simulated_residual_error_ratio.pdf",width=5.5,height=2.5+0.25*nscen)
par(mar=c(3.5,8.5,1,1),mgp=c(2,.5,0),tck=-0.02,pch=16,cex.lab=1.1,lheight=.9)
xlim<-c(-1.1*max(log(abs(quants))),1.1*max(log(abs(quants))))
plot(NA,NA,yaxt="n",ylab="",xlab="log(residual error ratio)",ylim=c(0.5,.5+n_scen_new),xlim=xlim)
abline(v=0,lty=3)
for(j in 1:n_scen_new) { 
#boxplot(log(as.numeric(myplot[j,])),at=j,add=T,lty=1,boxwex=0.9,cex=0.25,lwd=0.5,yaxt="n",col=cols[j],outline=F,horizontal=T) 
vioplot(log(as.numeric(myplot[j,])),at=j,add=T,col=cols[j],pchMed=16,colMed=1,wex=0.8,horizontal=T,border=NA)
}
axis(side=2,at=seq(n_scen_new),labels=scenario_names,las=2,cex.axis=0.9)
dev.off()

##=======================================================================##
##=========================================## plot ratio of autocorrelation 
##=======================================================================##
myplot<-rho_ratio
drop<-length(myplot[myplot<=0])
myplot[myplot<=0]<-0.1
pdrop<-drop/ntot ## proportion negative not included in log-ratio
select_scen<-c(1:nscen) ## select scenarios
n_scen_new<-length(select_scen)
myplot<-myplot[select_scen,]
quants<-apply(myplot,1,function(x) quantile(x,probs=c(0.01,0.99),na.rm=T))
cols<-colorRampPalette(brewer.pal(9,"Set1"))(n_scen_new)
##---------------------------------------------------## boxplots
pdf("observed_vs_simulated_autocorrelation_ratio.pdf",width=5.5,height=2.5+0.25*nscen)
par(mar=c(3.5,8.5,1,1),mgp=c(2,.5,0),tck=-0.02,pch=16,cex.lab=1.1,lheight=.9)
xlim<-c(-1.1*max(log(abs(quants))),1.1*max(log(abs(quants))))
plot(NA,NA,yaxt="n",ylab="",xlab="log(autocorrelation ratio)",ylim=c(0.5,.5+n_scen_new),xlim=xlim)
abline(v=0,lty=3)
for(j in 1:n_scen_new) { 
#boxplot(log(as.numeric(myplot[j,])),at=j,add=T,lty=1,boxwex=0.9,cex=0.25,lwd=0.5,yaxt="n",col=cols[j],outline=F,horizontal=T) 
vioplot(log(as.numeric(myplot[j,])),at=j,add=T,col=cols[j],pchMed=16,colMed=1,wex=0.8,horizontal=T,border=NA)
}
axis(side=2,at=seq(n_scen_new),labels=scenario_names,las=2,cex.axis=0.9)
dev.off()

##========================================================================##
##========================================================================##
##========================================================================##