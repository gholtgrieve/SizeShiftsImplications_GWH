##========================================================================##
##                                                                        ##
##                 Plot simulation model results for MSE                  ##
##                                                                        ##
##========================================================================##
##===============================================================## packages
pkgs<-c("here","tidyverse","dplyr","scales","RColorBrewer","readxl","gsl","vioplot", "caroline","ggplot2","viridis","forcats","gridExtra","lemon","fmsb","MetBrewer","ggradar","showtext","DataCombine","ggridges","ggExtra","remotes")
if(length(setdiff(pkgs,rownames(installed.packages())))>0) {install.packages(setdiff(pkgs,rownames(installed.packages())),dependencies=TRUE)}
invisible(lapply(pkgs,library,character.only=T))
## install.packages("ggpubr");library(ggpubr)
## remotes::install_github("ricardo-bion/ggradar");library(ggradar)

##==============================================================## functions
'%!in%'<-function(x,y)!('%in%'(x,y))

##============================================================## directories
homedir<-here();plotdir<-here("plots");outdir<-here("out");scendir<-here("scenarios")

##========================================================================##
##========================================================================##
##===========================================================## load results
##========================================================================##
##========================================================================##
setwd(outdir)
my.names<-dir();myfile<-my.names[length(my.names)] ## last file saved
load(myfile)
my.names<-dir();myfile<-my.names[length(my.names)] ## last file saved
timestamp<-substr(myfile,5,nchar(myfile)-6)
##============================================================## directories
## need to reset after switching to new computer as it was saved in output
homedir<-here();plotdir<-here("plots");outdir<-here("out");scendir<-here("scenarios")
##=============================================================## parameters
parameters<-readRDS(paste0("run_",timestamp,"_parms.Rdata"))
npara<-length(parameters)
for(i in 1:npara) { assign(names(parameters)[i],unlist(parameters[i])) }
##==============================================================## scenarios
setwd(scendir);scenarios<-read_excel("scenarios.xlsx")
# scenarios<-read_excel(paste0("run_",timestamp,"_scen.xlsx"))[,-1] ## saved
scenarios_all<-data.frame(scenarios)
##==============================================================## directory
setwd(plotdir)
dir.create(file.path(timestamp),showWarnings=F)
setwd(file.path(timestamp))

##========================================================================##
##========================================================================##
##========================================================================##
##=============## 1st set of scenarios on estimation methods and selectivity
##========================================================================##
##========================================================================##
##========================================================================##
scenarios<-scenarios_all[scenarios_all$factorMSY==1,]
nscen<-dim(scenarios)[1]		
##---------------------------------------------------------## scenario names
scenario_names<-paste0(scenarios$trends," ",scenarios$mgmt," ", scenarios$selectivity)
##--------------------------------------------------------## names of trends
trend_names<-unique(scenarios$trends)
n_trends<-length(trend_names)
##-------------------------------## names of estimation methods used in mgmt
mgmt_names<-unique(scenarios$mgmt)
n_mgmt<-length(mgmt_names)
##-------------------------------------------------## names on selectivities
selectivity_names<-unique(scenarios$selectivity)
n_select<-length(selectivity_names)

##========================================================================##
##========================================================================##
##=============## simulated mean egg mass in each stochastic run by scenario
##========================================================================##
##========================================================================##
nyobs<-dim(data.frame(obs.list[[1]][[1]]))[1] ## number of reconstructed yrs
y_obs_all<-seq(nyobs)
##--------------------------------------------------## array to save results
scen_names<-paste0("scen.",seq(nscen))
iter_names<-paste0("iter.",seq(niter))
year_names<-paste0("year.",seq(nyobs))
dnames<-list(scen_names,iter_names,year_names)
myarray<-array(dim=c(nscen,niter,nyobs),dimnames=dnames)
av_eggmass_list<-myarray ## scenario/iteration/year
##------------------------------------------## loop scenarios and iterations
for(j in 1:nscen) { 
  for(k in 1:niter) { 
    dat_jk<-data.frame(obs.list[[j]][[k]])$AvPerSpawnerEggmass[y_obs_all] 
    if(is.null(dat_jk)) { next } else { av_eggmass_list[j,k,]<-dat_jk }
  } 
} 
##--------------------------------------------------------## select scenario
ind<-as.numeric(scenarios$No[scenarios$selectivity=="unselective" & scenarios$mgmt=="smsy_goal" & scenarios$trends=="age-sex-length trends"])
av_eggmass_plot<-data.frame(av_eggmass_list[ind,,1:nyh])
##-------------------------------------------------------------------## plot
pdf("change_in_mean_egg_mass_per_spawner_iterations.pdf",width=4,height=4)
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2,0.5,0),cex.lab=1.2,cex.axis=1) 
plot(NA,NA,xlim=c(0,nyh),ylim=c(0,1.25e3),xlab="Year",ylab="Mean egg mass per spawner (g)")
cols<-c("goldenrod1","chocolate2","firebrick2","darkorchid3","midnightblue", "dodgerblue2","forestgreen")
matlines(t(av_eggmass_plot),lwd=0.05,lty=1,col=alpha(cols,1))
#lines(seq(nyobs),av_eggmass_plot[1,],lwd=2,lty=1,col=1)
dev.off()

pdf("change_in_mean_egg_mass_per_spawner_single_draw.pdf",width=4,height=4)
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2,0.5,0),cex.lab=1.2,cex.axis=1) 
plot(NA,NA,xlim=c(0,nyh),ylim=c(0,1.25e3),xlab="Year",ylab="Mean egg mass per spawner (g)")
cols<-c("goldenrod1","chocolate2","firebrick2","darkorchid3","midnightblue", "dodgerblue2","forestgreen")
lines(seq(nyh),av_eggmass_plot[1,],lwd=2,lty=1,col=1)
dev.off()

##========================================================================##
##========================================================================##
##==================================## trends in mean fecundity and egg mass
##========================================================================##
##========================================================================##
nyrs<-10 ## number of years for calculating difference (3,5,10)
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

##==================================## function to plot median and quantiles
## using 90% CIs (5th-95th quantile)
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) }
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) }

##=================================================================## ggplot
df_scen<-dplyr::select(scenarios,trends,selectivity,mgmt)
df_egg_mass<-data.frame(cbind(df_scen,egg_trends))
df_egg_mass<-dplyr::filter(df_egg_mass,mgmt=="smsy_goal")
df_egg_mass<-dplyr::filter(df_egg_mass,trends!="continuing trends")
df_egg_mass$trends<-as.factor(as.character(df_egg_mass$trends))
plot_egg_mass<-df_egg_mass %>% pivot_longer(!c(trends,selectivity,mgmt), names_to="iteration", values_to="value") %>% data.frame()
##-----------------------------------------------------------## selectivity
s_levels<-c("6 inch gillnet","unselective","8.5 inch gillnet")
#plot_egg_mass$sel_ordered<-factor(plot_egg_mass$selectivity,levels=s_levels)
s_labels<-c("small-mesh","unselective","large-mesh")
plot_egg_mass$sel_ordered<-factor(plot_egg_mass$selectivity,levels=s_levels,labels=s_labels)
plot_egg_mass<-dplyr::select(plot_egg_mass,-iteration,-mgmt,-selectivity)
##----------------------------------## function to plot median and quantiles
## using 90% CIs (5th-95th quantile)
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) }
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) }
##-----------------------------------------------------------------## ggplot
jig<-position_dodge(width=0.6)
p<-plot_egg_mass %>% ggplot(aes(x=fct_inorder(trends),y=value,fill=sel_ordered))+
  geom_violin(position=jig,trim=F,lwd=0.1,scale="width",width=0.6)+ 
  #geom_hline(yintercept=-49,linewidth=0.1)+ ## Staton et al. 2021
  #geom_hline(yintercept=-28,linewidth=0.1)+ ## Ohlberger et al. 2020
  geom_hline(yintercept=0,linetype="dashed",linewidth=0.2)+ 
  geom_vline(xintercept=1.5,linetype="solid",linewidth=0.1)+
  geom_vline(xintercept=2.5,linetype="solid",linewidth=0.1)+
  stat_summary(fun.data=summary_CI90,position=jig,linewidth=0.2) +
  stat_summary(fun.data=summary_CI50,position=jig,linewidth=0.6) +
  scale_fill_manual(values=c("firebrick","chocolate3","goldenrod1"))+
  scale_y_continuous(breaks=seq(-100,100,20))+
  coord_cartesian(ylim=c(-70,40))+
  theme_classic() +
  labs(x="",y="Change in mean egg mass per spawner (%)")+ 
  labs(fill="Selectivity")+
  theme(strip.background=element_blank(),
        axis.line=element_line(linewidth=0.1),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.4,'cm'),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        legend.position=c(0.82,0.86)
  )
p<-p+guides(fill=guide_legend(override.aes=list(size=0.3)))
ggsave(paste0("changes_in_eggmass_",nyrs,"yrs.pdf"),p,width=4.2,height=5.5,units="in")

##=====================================================## ggplot all methods
df_scen<-dplyr::select(scenarios,trends,selectivity,mgmt)
df_egg_mass<-data.frame(cbind(df_scen,egg_trends))
df_egg_mass<-dplyr::filter(df_egg_mass,trends!="continuing trends")
df_egg_mass$trends<-as.factor(as.character(df_egg_mass$trends))
plot_egg_mass<-df_egg_mass %>% pivot_longer(!c(trends,selectivity,mgmt), names_to="iteration", values_to="value") %>% data.frame()
# s_levels<-c("6 inch gillnet", "unselective", "8.5 inch gillnet")
#plot_egg_mass$sel_ordered<-factor(plot_egg_mass$selectivity,levels=s_levels)
s_labels<-c("small-mesh","unselective","large-mesh")
plot_egg_mass$sel_ordered<-factor(plot_egg_mass$selectivity,levels=s_levels,labels=s_labels)
plot_egg_mass<-dplyr::select(plot_egg_mass,-iteration,-selectivity)
##----------------------------------## function to plot median and quantiles
## using 90% CIs (5th-95th quantile)
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) }
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) }
##-----------------------------------------------------------------## ggplot
jig<-position_dodge(width=0.6)
p<-plot_egg_mass %>% ggplot(aes(x=fct_inorder(trends),y=value,fill=sel_ordered))+
  geom_violin(position=jig,trim=F,lwd=0.1,scale="width",width=0.6)+ 
  geom_hline(yintercept=0,linetype="dashed",size=0.2)+ 
  geom_vline(xintercept=1.5,linetype="solid",size=0.1)+
  geom_vline(xintercept=2.5,linetype="solid",size=0.1)+
  stat_summary(fun.data=summary_CI90,position=jig,size=0.2) +
  stat_summary(fun.data=summary_CI50,position=jig,size=0.6) +
  scale_fill_manual(values=c("firebrick","chocolate3","goldenrod1"))+
  scale_y_continuous(breaks=seq(-100,100,20))+
  coord_cartesian(ylim=c(-70,40))+
  theme_classic() +
  labs(x="",y="Change in mean egg mass per spawner (%)")+ 
  labs(fill="Selectivity")+
  theme(strip.background=element_blank(),
        axis.line=element_line(linewidth=0.1),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.4,'cm'),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        #legend.position=c(0.82,0.86)
  )
col.labs<-c("YPR","DLM","TRM")
names(col.labs)<-c("s_eq_goal","smsy_dlm_goal","smsy_goal")
g<-p+facet_grid(cols=vars(mgmt),scale="free",labeller=labeller(mgmt=col.labs))+theme(strip.text.x=element_text(size=12,face="bold"),legend.position=c(0.942,0.86))+guides(fill=guide_legend(override.aes=list(size=0.3)))
ggsave(paste0("changes_in_eggmass_",nyrs,"yrs_all_methods.pdf"),g,width=11, height=5.5,units="in")


##========================================================================##
##========================================================================##
##=========================================## differences in S_msy estimates
##========================================================================##
##========================================================================##
review_years<-seq((nyi+20),ny,goalfreq)
nrev<-length(review_years)
myarray<-array(NA,dim=c(nscen,niter,nrev))
S_msy_scen<-myarray
for(j in 1:nscen) { 
  for(k in 1:niter) { 
    s_msy_est<-as.vector(unlist(S_msy.list[[j]][[k]]))
    if(length(s_msy_est)!=0) S_msy_scen[j,k,]<-s_msy_est
  } ## end j-loop
} ## end i-loop
##------------------------------------------------------------## time period
## evaluated after 50 years, i.e. after the historical period
S_msy_hist<-S_msy_scen[,,which(review_years==nyh)]
##------------------------------------------------------## S_msy differences
## DLM and YPR analysis compared to traditional model
ref_msy_ind<-which(grepl("smsy_goal",scenario_names,fixed=T)) 
ypr_msy_ind<-which(grepl("s_eq_goal",scenario_names,fixed=T)) 
dlm_msy_ind<-which(grepl("smsy_dlm_goal",scenario_names,fixed=T))  
n_diffs<-length(ref_msy_ind) ## number of scenarios to compare
S_msy_diff_YPR<-S_msy_diff_DLM<-array(NA,dim=c(n_diffs,niter))
for(i in 1:n_diffs){
ref_msys<-S_msy_hist[ref_msy_ind[i],]
ypr_msys<-S_msy_hist[ypr_msy_ind[i],]
dlm_msys<-S_msy_hist[dlm_msy_ind[i],]
S_msy_diff_YPR[i,]<-100*(ypr_msys-ref_msys)/ref_msys
S_msy_diff_DLM[i,]<-100*(dlm_msys-ref_msys)/ref_msys
}
##-----------------------------------------------## add scenario information
use_all<-dplyr::select(scenarios,trends,selectivity,mgmt)
S_msy_df<-data.frame(cbind(use_all,S_msy_hist))
S_msy_df<-dplyr::filter(S_msy_df,trends!="ASL trends continued")
S_msy_df<-dplyr::filter(S_msy_df,trends!="continuing trends")
S_msy_df$trends[S_msy_df$trends=="ASL trends stabilized"]<-"age-sex-length trends"
use_ypr<-use_all[ypr_msy_ind,]
S_msy_df_YPR<-data.frame(cbind(use_ypr,S_msy_diff_YPR))
S_msy_df_YPR<-dplyr::filter(S_msy_df_YPR,trends!="ASL trends continued")
S_msy_df_YPR<-dplyr::filter(S_msy_df_YPR,trends!="continuing trends")
S_msy_df_YPR$trends[S_msy_df_YPR$trends=="ASL trends stabilized"]<-"age-sex-length trends"
use_dlm<-use_all[dlm_msy_ind,]
S_msy_df_DLM<-data.frame(cbind(use_dlm,S_msy_diff_DLM))
S_msy_df_DLM<-dplyr::filter(S_msy_df_DLM,trends!="ASL trends continued")
S_msy_df_DLM<-dplyr::filter(S_msy_df_DLM,trends!="continuing trends")
S_msy_df_DLM$trends[S_msy_df_DLM$trends=="ASL trends stabilized"]<-"age-sex-length trends"
n_diff<-dim(S_msy_df_DLM)[1]
##--------------------------------------------------------## ordered factors
s_levels<-c("6 inch gillnet", "unselective", "8.5 inch gillnet")
S_msy_df$sel_ordered<-factor(S_msy_df$selectivity,levels=s_levels)
S_msy_df_YPR$sel_ordered<-factor(S_msy_df_YPR$selectivity,levels=s_levels)
S_msy_df_DLM$sel_ordered<-factor(S_msy_df_DLM$selectivity,levels=s_levels)
S_msy_df<-dplyr::select(S_msy_df,-selectivity)
S_msy_df_YPR<-dplyr::select(S_msy_df_YPR,-selectivity)
S_msy_df_DLM<-dplyr::select(S_msy_df_DLM,-selectivity)
##-------------------------------------------## add which alternative method
S_msy_df_YPR$method<-"YPR"
S_msy_df_DLM$method<-"DLM"
S_msy_diffs<-data.frame(rbind(S_msy_df_YPR,S_msy_df_DLM))
S_msy_diffs<-dplyr::select(S_msy_diffs,-mgmt)

##============================================================## select data
plot_smsy<-S_msy_diffs %>% pivot_longer(!c(trends,sel_ordered,method), names_to="iteration", values_to="value") %>% data.frame()
plot_smsy<-plot_smsy[complete.cases(plot_smsy),] ## <1%
plot_smsy<-dplyr::select(plot_smsy,-iteration)

##================================================================## medians
S_msy_scen<-dplyr::select(S_msy_diffs,sel_ordered,method,trends)
S_msy_iter<-dplyr::select(S_msy_diffs,-sel_ordered,-method,-trends)
S_msy_med<-apply(S_msy_iter,1,function(x) median(x,na.rm=T))
S_msy_diff_med<-data.frame(cbind(S_msy_scen,median=S_msy_med))

##===========================================================## for plotting
method_colors<-c("dodgerblue2", "goldenrod1")
jig<-position_dodge(width=0.75)
##----------------------------------## function to plot median and quantiles
## using 90% CIs (5th-95th quantile)
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) }
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) }

##=========================================================## plot quantiles
method_colors<-c("dodgerblue2", "goldenrod1")
jig<-position_dodge(width=0.5)
colors<-rep(method_colors,n_diff)
ymin<-quantile(plot_smsy$value,probs=0.01,na.rm=T)
ymax<-quantile(plot_smsy$value,probs=0.99,na.rm=T)
p<-plot_smsy %>% ggplot(aes(x=fct_inorder(trends),y=value,fill=fct_inorder(method)))+
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  stat_summary(fun.data=summary_CI90,position=jig,size=0.2,color=colors)+
  stat_summary(fun.data=summary_CI50,position=jig,size=0.6,color=colors)+
  # scale_y_continuous(limits=c(ymin,ymax),breaks=seq(-50,150,25))+
  # ylim(-20,30)+
  ## Setting limits affects medians/CIs as it drops data from stat_summary
  ## Do not use 'ylim' or scale with 'limits' > only use coord_cartesian
  scale_y_continuous(breaks=seq(-100,150,25))+
  coord_cartesian(ylim=c(-55,80))+
  theme_classic()+
  labs(fill="Method")+
  labs(x="",y="Difference in S_MSY (%)\nfrom time-invariant model")+ 
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=10),
        panel.border=element_rect(fill=NA,size=1),
        legend.position="none",
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8)
  )
col.labs<-c("large-mesh","unselective","small-mesh")
names(col.labs)<-c("8.5 inch gillnet","unselective","6 inch gillnet")
g<-p+facet_grid(cols=vars(sel_ordered),scale="free",labeller=labeller(sel_ordered=col.labs))+theme(strip.text.x=element_text(size=10))
ggsave("estimated_S_msy_by_method_quantiles.pdf",g,width=5,height=4,units="in")  

##===========================================================## plot medians
cols<-c("dodgerblue2","goldenrod1") ## when using fct_inorder(method)
ymin<-min(S_msy_diff_med$median,na.rm=T)*1.1
ymax<-max(S_msy_diff_med$median,na.rm=T)*1.1
p<-S_msy_diff_med %>% 
  filter(method=="YPR") %>%
  ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=fct_inorder(method)))+
 geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=0.5) +
  #geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point(size=2.5) + 
  scale_y_continuous(breaks=seq(-100,150,10))+
  coord_cartesian(ylim=c(-20,30))+
  scale_colour_manual(values=cols)+
  theme_classic()+
  labs(x="",y="Difference in S_MSY (%)\nfrom time-invariant model")+ 
  labs(color="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=10),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(cols=vars(sel_ordered),scale="free")+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10),legend.position=c(0.1,0.8)) 
ggsave("estimated_S_msy_by_method_medians.pdf",g,width=5,height=4,units="in")  

##========================================================================##
##====================## probability of Smsy greater than traditional method
##========================================================================##
S_msy_probs<-dplyr::select(S_msy_diffs,-sel_ordered,-method,-trends)
S_msy_prob<-apply(S_msy_probs,1,function(x) length(x[which(x>0)])/length(x))
S_msy_scen<-dplyr::select(S_msy_diffs,sel_ordered,method,trends)
S_msy_prob<-data.frame(cbind(S_msy_scen,S_msy_prob))
S_msy_prob<-dplyr::filter(S_msy_prob,trends!="continuing trends")

##=============================================================## plot probs
cols<-c("dodgerblue2","goldenrod1") ## when using fct_inorder(method)
p<-S_msy_prob %>% ggplot(aes(x=fct_inorder(trends),y=S_msy_prob,col=fct_inorder(method),group=fct_inorder(method)))+
  geom_hline(yintercept=0.5,linetype="solid",size=0.1)+ 
  geom_line(lwd=0.5) +
  #geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point(size=2.5) + 
  coord_cartesian(ylim=c(0,1))+
  scale_colour_manual(values=cols)+
  theme_classic()+
  labs(x="",y="Probability of estimate greater\nthan from time-invariant model")+ 
  labs(color="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=10),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(cols=vars(sel_ordered),scale="free")+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10),legend.position=c(0.1,0.8))
ggsave("estimated_S_msy_by_method_probability.pdf",g,width=5,height=4,units="in")  

##========================================================================##
##========================================================================##
##=========================================## management performance metrics
##========================================================================##
##========================================================================##
## using simulated or observed data?
##---------------------------------## number of reconstructed observed years
ny_obs<-dim(data.frame(obs.list[[1]][[1]]))[1]
##--------------------------------------------------------## new trend names
trend_names_old<-trend_names
trend_names<-c("no trends","age-length trends","ASL trends stabilized","ASL trends continued" )
replace<-data.frame(from=trend_names_old,to=trend_names)
scenarios<-FindReplace(scenarios,Var="trends",replaceData=replace, from="from",to="to",exact=FALSE)
##------------------------------------------------------## array for results
scen_names<-paste0("scen=",seq(nscen))
iter_names<-paste0("iter=",seq(niter))
dnames<-list(scen_names,iter_names)
myarray<-array(dim=c(nscen,niter),dimnames=dnames)

##========================================================================##
##===========================================## performance trends over time
##========================================================================##
year_index_list<-list()
future_yrs<-seq(10,50,10)
n_index<-length(future_yrs)
##------------------------------------------------------------## each period
#for(i in 1:n_index) { year_index_list[[i]]<-(nyh+1+(i-1)*10):(nyh+i*10) }
##-------------------------------------------------------------## cumulative
for(i in 1:n_index) { year_index_list[[i]]<-(nyh+1):(nyh+i*10) }
##==========================================## loop scenarios and iterations
av_harv_list<-av_ret_list<-av_esc_list<-list()
av_harv_diff_list<-av_ret_diff_list<-av_esc_diff_list<-list()
for(i in 1:n_index) { 
  year_index<-year_index_list[[i]]
  av_harv<-av_ret<-av_esc<-myarray  
  for(j in 1:nscen) {
    for(k in 1:niter) { 
      harv_jk<-data.frame(obs.list[[j]][[k]])$obsHarv[year_index] 
      if(is.null(harv_jk)){next}else{ av_harv[j,k]<-mean(harv_jk,na.rm=T)}
      esc_jk<-data.frame(obs.list[[j]][[k]])$obsEsc[year_index] 
      if(is.null(esc_jk)){next}else{ av_esc[j,k]<-mean(esc_jk,na.rm=T)}
      ret_jk<-data.frame(obs.list[[j]][[k]])$obsRet[year_index] 
      if(is.null(ret_jk)){next}else{ av_ret[j,k]<-mean(ret_jk,na.rm=T)}
    } ## end k-loop
  } ## end j-loop
  av_harv_list[[i]]<-av_harv
  av_esc_list[[i]]<-av_esc
  av_ret_list[[i]]<-av_ret
  ## compute differences relative to traditional method as reference
  #print(paste0("done i=",i," part1"))
  av_harv_diff<-av_ret_diff<-av_esc_diff<-myarray  
  for(j in 1:nscen) {
    for(k in 1:niter) {
    if(j %in% ref_msy_ind) j_ref<-j
    if(j %in% ypr_msy_ind) j_ref<-j-4
    if(j %in% dlm_msy_ind) j_ref<-j-8
    av_harv_diff[j,k]<-(av_harv[j,k]-av_harv[j_ref,k])*100/av_harv[j_ref,k]
    av_esc_diff[j,k]<-(av_esc[j,k]-av_esc[j_ref,k])*100/av_esc[j_ref,k]
    av_ret_diff[j,k]<-(av_ret[j,k]-av_ret[j_ref,k])*100/av_ret[j_ref,k]
    } ## end k-loop
  } ## end j-loop
  av_harv_diff_list[[i]]<-av_harv_diff
  av_esc_diff_list[[i]]<-av_esc_diff
  av_ret_diff_list[[i]]<-av_ret_diff
  #print(paste0("done i=",i," part2"))
} ## end i-loop over year indices

##==============================================================## scenarios
df_diff_scen<-dplyr::select(scenarios,trends,selectivity,mgmt)
##----------------------------------------------------## 'method'
df_diff_scen$method<-"TRM"
df_diff_scen$method[grepl("eq",df_diff_scen$mgmt)]<-"YPR"
df_diff_scen$method[grepl("dlm",df_diff_scen$mgmt)]<-"DLM"

##================================================================## harvest
df_av_harv<-array(dim=c(nscen,n_index))
for(i in 1:n_index) { df_av_harv[,i]<-apply(av_harv_diff_list[[i]],1, function(x) median(x,na.rm=T)) }
colnames(df_av_harv)<-future_yrs
##-------------------------------------------------## combine with scenarios
df_av_harv<-data.frame(cbind(df_diff_scen,df_av_harv))
df_av_harv<-dplyr::filter(df_av_harv,method!="TRM")
##----------------------------------------------------------## select trends
df_av_harv<-dplyr::filter(df_av_harv,trends=="ASL trends stabilized")
df_av_harv<-dplyr::select(df_av_harv,-trends,-mgmt)
##-----------------------------------------------------------## pivot longer
df_av_harv_plot<-df_av_harv %>% pivot_longer(!c(selectivity,method), names_to="period", values_to="median") %>% data.frame() 
df_av_harv_plot$period<-as.numeric(gsub("X","",df_av_harv_plot$period))
##-------------------------------------------------------## ordered factors
#df_av_harv_plot$sel_ordered<-factor(df_av_harv_plot$selectivity,levels=c("6 inch gillnet", "unselective", "8.5 inch gillnet"))
df_av_harv_plot$sel_ordered<-factor(df_av_harv_plot$selectivity,levels=c("6 inch gillnet", "unselective", "8.5 inch gillnet"),labels=c("small-mesh","unselective","large-mesh"))
##-----------------------------------------------------------## metric label
df_av_harv_plot$metric_label<-"Mean harvest"

##=================================================================## return
df_av_ret<-array(dim=c(nscen,n_index))
for(i in 1:n_index) {df_av_ret[,i]<-apply(av_ret_diff_list[[i]],1, function(x) median(x,na.rm=T)) }
colnames(df_av_ret)<-future_yrs
##-------------------------------------------------## combine with scenarios
df_av_ret<-data.frame(cbind(df_diff_scen,df_av_ret))
df_av_ret<-dplyr::filter(df_av_ret,method!="TRM")
##----------------------------------------------------------## select trends
df_av_ret<-dplyr::filter(df_av_ret,trends=="ASL trends stabilized")
df_av_ret<-dplyr::select(df_av_ret,-trends,-mgmt)
##-----------------------------------------------------------## pivot longer
df_av_ret_plot<-df_av_ret %>% pivot_longer(!c(selectivity,method), names_to="period", values_to="median") %>% data.frame() 
df_av_ret_plot$period<-as.numeric(gsub("X","",df_av_ret_plot$period))
##-------------------------------------------------------## ordered factors
#df_av_ret_plot$sel_ordered<-factor(df_av_ret_plot$selectivity,levels=c("6 inch gillnet", "unselective", "8.5 inch gillnet"))
df_av_ret_plot$sel_ordered<-factor(df_av_ret_plot$selectivity,levels=c("6 inch gillnet", "unselective", "8.5 inch gillnet"),labels=c("small-mesh","unselective","large-mesh"))
##-----------------------------------------------------------## metric label
df_av_ret_plot$metric_label<-"Mean return"

##=================================================================## return
df_av_esc<-array(dim=c(nscen,n_index))
for(i in 1:n_index) {df_av_esc[,i]<-apply(av_esc_diff_list[[i]],1, function(x) median(x,na.rm=T)) }
colnames(df_av_esc)<-future_yrs
##-------------------------------------------------## combine with scenarios
df_av_esc<-data.frame(cbind(df_diff_scen,df_av_esc))
df_av_esc<-dplyr::filter(df_av_esc,method!="TRM")
##----------------------------------------------------------## select trends
df_av_esc<-dplyr::filter(df_av_esc,trends=="ASL trends stabilized")
df_av_esc<-dplyr::select(df_av_esc,-trends,-mgmt)
##-----------------------------------------------------------## pivot longer
df_av_esc_plot<-df_av_esc %>% pivot_longer(!c(selectivity,method), names_to="period", values_to="median") %>% data.frame() 
df_av_esc_plot$period<-as.numeric(gsub("X","",df_av_esc_plot$period))
##-------------------------------------------------------## ordered factors
# df_av_esc_plot$sel_ordered<-factor(df_av_esc_plot$selectivity,levels=c("6 inch gillnet", "unselective", "8.5 inch gillnet"))
df_av_esc_plot$sel_ordered<-factor(df_av_esc_plot$selectivity,levels=c("6 inch gillnet", "unselective", "8.5 inch gillnet"),labels=c("small-mesh","unselective","large-mesh"))
##-----------------------------------------------------------## metric label
df_av_esc_plot$metric_label<-"Mean escapement"

##=================================================================## return
df_plot<-df_av_ret_plot
p<-df_plot %>% ggplot(aes(x=period,y=median,col=fct_inorder(method))) +
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point() + 
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  #coord_cartesian(xlim=c(5,55))+
  theme_classic() +
  labs(x="Year",y="Median % difference\nfrom time-invariant model") + 
  labs(fill="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        legend.key.size=unit(0.4,'cm'),
        panel.border=element_rect(fill=NA,size=1)
  )
g<-p+facet_grid(cols=vars(sel_ordered),rows=vars(metric_label),scale="free") +theme(strip.text.x=element_text(size=9),strip.text.y=element_text(size=11), legend.position=c(0.08,0.82),legend.title=element_blank(),legend.text=element_text(size=9),legend.key=element_rect(colour=NA,fill=NA))
ggsave("Future_mean_return.pdf",g,width=6,height=2.5,units="in")

##=================================================## harvest and escapement
df_plot<-data.frame(rbind(df_av_esc_plot,df_av_harv_plot))
p<-df_plot %>% ggplot(aes(x=period,y=median,col=fct_inorder(method))) +
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point() + 
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  #coord_cartesian(xlim=c(5,55))+
  theme_classic() +
  labs(x="Year",y="Median % difference from time-invariant model") + 
  labs(fill="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        legend.key.size=unit(0.4,'cm'),
        panel.border=element_rect(fill=NA,size=1)
  )
g<-p+facet_grid(cols=vars(sel_ordered),rows=vars(metric_label),scale="free") +theme(strip.text.x=element_text(size=9),strip.text.y=element_text(size=11), legend.position=c(0.08,0.9),legend.title=element_blank(),legend.text=element_text(size=9),legend.key=element_rect(colour=NA,fill=NA))
ggsave("Future_mean_harvest_escapement.pdf",g,width=6,height=4,units="in")

##========================================================================##
##==========================================## performance for last 50 years
##========================================================================##
##--------------------------------------------------------## years evaluated
# year_index<-1:ny_obs ## all reconstructed years
year_index<-(nyh+1):ny_obs ## most recent e.g. 50 years

##=============================================## for each estimation method
sum_log_harv<-cv_harv<-av_harv<-av_log_harv<-p_no_harv<-myarray ## harvest 
sum_log_esc<-cv_esc<-av_esc<-av_log_esc<-myarray ## escapement 
sum_log_ret<-cv_ret<-av_ret<-av_log_ret<-myarray ## total return 
p_below_thresh<-p_above_thresh<-myarray ## conservation threshold
thresh<-0.5 ## proportion of maximum recruitment at equilibrium
##------------------------------------------## loop scenarios and iterations
for(j in 1:nscen) {
  for(k in 1:niter) { 
    ##---------------------------------------------------------## escapement
    # esc_jk<-data.frame(data.list[[j]][[k]])$Esc[year_index] ## simulated 
    esc_jk<-data.frame(obs.list[[j]][[k]])$obsEsc[year_index] ## observed 
    if(is.null(esc_jk)) { next } else { 
      sum_log_esc[j,k]<-sum(log(esc_jk+0.0001))
      cv_esc[j,k]<-1/(sd(esc_jk,na.rm=T)/mean(esc_jk,na.rm=T)) ## 1/CV
      av_esc[j,k]<-mean(esc_jk,na.rm=T)
      av_log_esc[j,k]<-mean(log(esc_jk),na.rm=T)
    }
    ##----------------------------------------------------------## harvest
    # harv_jk<-data.frame(data.list[[j]][[k]])$Harv[year_index] ## simulated
    harv_jk<-data.frame(obs.list[[j]][[k]])$obsHarv[year_index] ## observed 
    if(is.null(harv_jk)) { next } else { 
      sum_log_harv[j,k]<-sum(log(harv_jk+0.0001))
      cv_harv[j,k]<-1/(sd(harv_jk,na.rm=T)/mean(harv_jk,na.rm=T)) ## 1/CV
      av_harv[j,k]<-mean(harv_jk,na.rm=T)
      av_log_harv[j,k]<-mean(log(harv_jk),na.rm=T)
      p_no_harv[j,k]<-length(which(harv_jk==0))/length(harv_jk)
    }
    ##----------------------------------------------------------## return
    # ret_jk<-data.frame(data.list[[j]][[k]])$Ret[year_index] ## simulated 
    ret_jk<-data.frame(obs.list[[j]][[k]])$obsRet[year_index] ## observed 
    if(is.null(ret_jk)) { next } else { 
      sum_log_ret[j,k]<-sum(log(ret_jk+0.0001))
      cv_ret[j,k]<-1/(sd(ret_jk,na.rm=T)/mean(ret_jk,na.rm=T)) ## 1/CV
      av_ret[j,k]<-mean(ret_jk,na.rm=T)
      av_log_ret[j,k]<-mean(log(ret_jk),na.rm=T)
      ## maximum recruitment given ricker parameters
      ricker_parms<-data.frame(sr_sim.list[[j]][[k]]) 
      max_rec<-(ricker_parms$alpha/ricker_parms$beta)*exp(-1) ## Hilborn ~
      ## probability above/below threshold
      p_below_jk<-length(which(ret_jk<max_rec*thresh))/length(ret_jk)
      p_below_thresh[j,k]<-p_below_jk 
      p_above_jk<-length(which(ret_jk>max_rec*thresh))/length(ret_jk)
      p_above_thresh[j,k]<-1-p_below_jk
    }
  } ## end k-loop
} ## end j-loop
p_below_thresh[is.na(p_below_thresh)]<-0

##=================================================## ratios and differences
av_harv_diff<-cv_harv_diff<-cv_harv_ratio<-av_log_harv_diff<-myarray
av_esc_diff<-av_log_esc_diff<-myarray
av_ret_diff<-av_log_ret_diff<-p_above_thresh_ratio<-p_no_harv_ratio<-myarray
##-----------------------------------------------------------## get indices
ref_msy_ind<-which(grepl("smsy_goal",scenario_names,fixed=T)) 
ypr_msy_ind<-which(grepl("s_eq_goal",scenario_names,fixed=T)) 
dlm_msy_ind<-which(grepl("smsy_dlm_goal",scenario_names,fixed=T))  
##---------------------------------------------------## compute differences
for(j in 1:nscen) {
  for(k in 1:niter) {
    if(j %in% ref_msy_ind) j_ref<-j
    if(j %in% ypr_msy_ind) j_ref<-j-4
    if(j %in% dlm_msy_ind) j_ref<-j-8
    cv_harv_diff[j,k]<-(cv_harv[j,k]-cv_harv[j_ref,k])*100/cv_harv[j_ref,k]
    cv_harv_ratio[j,k]<-cv_harv[j,k]/cv_harv[j_ref,k]
    av_harv_diff[j,k]<-(av_harv[j,k]-av_harv[j_ref,k])*100/av_harv[j_ref,k]
    av_log_harv_diff[j,k]<-(av_log_harv[j,k]-av_log_harv[j_ref,k])
    av_esc_diff[j,k]<-(av_esc[j,k]-av_esc[j_ref,k])*100/av_esc[j_ref,k]
    av_log_esc_diff[j,k]<-(av_log_esc[j,k]-av_log_esc[j_ref,k])
    av_ret_diff[j,k]<-(av_ret[j,k]-av_ret[j_ref,k])*100/av_ret[j_ref,k]
    av_log_ret_diff[j,k]<-(av_log_ret[j,k]-av_log_ret[j_ref,k])
    p_above_thresh_ratio[j,k]<-p_above_thresh[j,k]/p_above_thresh[j_ref,k]
    p_no_harv_ratio[j,k]<-p_no_harv[j,k]/p_no_harv[j_ref,k]
  } ## end k-loop
} ## end j-loop

##========================================================================##
##======================## plot metrics that compare to time-invariant model
##========================================================================##
df_diff<-dplyr::select(scenarios,trends,selectivity,mgmt)
##----------------------------------------------------## 'method'
df_diff$method<-"TRM"
df_diff$method[grepl("eq",df_diff$mgmt)]<-"YPR"
df_diff$method[grepl("dlm",df_diff$mgmt)]<-"DLM"
##----------------------------------------------------## median all metrics
df_diff$av_harv<-apply(av_harv_diff,1,function(x) median(x,na.rm=T))
df_diff$av_esc<-apply(av_esc_diff,1,function(x) median(x,na.rm=T))
df_diff$av_ret<-apply(av_ret_diff,1,function(x) median(x,na.rm=T))
df_diff$p_harv_above_ref<-apply(av_harv_diff,1,function(x) length(x[which(x>0)])/length(x))
df_diff$p_esc_above_ref<-apply(av_esc_diff,1,function(x) length(x[which(x>0)])/length(x))
df_diff$p_ret_above_ref<-apply(av_ret_diff,1,function(x) length(x[which(x>0)])/length(x))
##----------------------------------------------------------## metric labels
labs<-c("Mean harvest","Mean escapement","Mean return","P(harv>harv_ref)","P(esc>esc_ref)","P(ret>ret_ref)")

##==========================================================## long format
df_diff_plot<-df_diff %>% pivot_longer(!c(trends,selectivity,mgmt,method), names_to="metric", values_to="median") %>% data.frame()
##-----------------------------------------------## add labels for metrics
metrics<-data.frame(name=colnames(df_diff)[-c(1:4)])
metrics$label<-labs
for(l in 1:dim(df_diff_plot)[1]) { df_diff_plot$metric_label[l]<-metrics$label[df_diff_plot$metric[l]==metrics$name] }
##-------------------------------------------------------## ordered factors
sel_levels<-c("6 inch gillnet", "unselective", "8.5 inch gillnet")
sel_labels<-c("small-mesh","unselective","large-mesh")
df_diff_plot$sel_ordered<-factor(df_diff_plot$selectivity,levels=sel_levels, labels=sel_labels)
##-------------------------------------------------------## drop references
## values are differences compared to 'traditional' MSY model 
df_diff_plot<-df_diff_plot[df_diff_plot$mgmt!="smsy_goal",]
df_diff_plot<-df_diff_plot[df_diff_plot$mgmt!="umsy_goal",]
##----------------------------------------------------------## select trends
df_diff_plot<-dplyr::filter(df_diff_plot,trends!="age-length trends")
df_diff_plot$trends<-as.factor(as.character(df_diff_plot$trends))
##-----------------------------------------## split ratios and % differences
df_diff_plot1<-df_diff_plot[df_diff_plot$metric %in% c("av_harv","av_esc"),]
df_diff_plot2<-df_diff_plot[df_diff_plot$metric %in% c("av_ret"),]
df_diff_plot3<-df_diff_plot[df_diff_plot$metric %in% c("p_ret_above_ref"),]

##============================================## plot harvest and escapement
p<-df_diff_plot1 %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point() + 
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="",y="Median % difference from time-invariant model") + 
  labs(color="")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank()
  )
g<-p+facet_grid(rows=vars(metric_label),cols=vars(sel_ordered),scale="free") +theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10),legend.position=c(0.08,0.94))
ggsave("Metrics_differences_harv_esc.pdf",g,width=5,height=5.5,units="in")

##============================================================## plot return
p<-df_diff_plot2 %>%
  filter(method!="DLM") %>%
  ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=1) +
  geom_point(size=2.5) + 
  geom_point(shape=1,size=3,fill=NA,color="black") +
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1),breaks=seq(-20,30,5)) +
  theme_classic() +
  labs(x="",y="Median % difference\nfrom time-invariant model") + 
  labs(color="")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank()
  )
g<-p+facet_grid(cols=vars(sel_ordered),scale="free") +theme(strip.text.x=element_text(size=10),legend.position=c(0.08,0.94))
ggsave("Metrics_differences_return.pdf",g,width=5,height=4,units="in")

##============================================================## plot return
p<-df_diff_plot2 %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),shape=fct_inorder(sel_ordered),linetype=fct_inorder(sel_ordered),group=interaction(method,sel_ordered))) +
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=0.8) +
  geom_point(size=2) + 
  scale_linetype_manual(values=c("solid","dashed","dotted")) +
  scale_shape_manual(values=c(15,17,19)) +
  scale_colour_manual(values=c("dodgerblue2","goldenrod1")) +
  scale_y_continuous(expand=c(0.1,0.1),breaks=seq(-20,30,5)) +
  scale_x_discrete(expand=c(0.1,0.1))+
  theme_classic() +
  labs(x="",y="Median % difference\nfrom time-invariant model") + 
  labs(col="Method",linetype="Selectivity",shape="Selectivity")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        legend.position="right",
        legend.key.size=unit(0.6,'cm'),
        legend.key.width=unit(0.8,'cm'),
        panel.border=element_rect(fill=NA,size=1)
  )
ggsave("Metrics_differences_return_v2.pdf",p,width=4.6,height=4,units="in")

##===========================================================## plot P(>ref)
p<-df_diff_plot3 %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_hline(yintercept=0.5,linetype="solid",size=0.1)+ 
  geom_line(lwd=1) +
  geom_point(size=2.5) + 
  geom_point(shape=1,size=3,fill=NA,color="black") +
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="",y="Probability of mean return larger\nthan for time-invariant model") + 
  labs(color="")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank()
  )
g<-p+facet_grid(cols=vars(sel_ordered),scale="free") +theme(strip.text.x=element_text(size=10),legend.position=c(0.08,0.94))
ggsave("Metrics_P_return_higher_ref.pdf",g,width=5,height=4,units="in")

##========================================================================##
##======================## differences all iterations harvest vs escapement
##========================================================================##
df_scen<-dplyr::select(scenarios,trends,selectivity,mgmt)
df_harv<-data.frame(cbind(df_scen,av_harv_diff))
df_harv_long<-df_harv %>% pivot_longer(!c(selectivity,mgmt,trends),names_to="iter",values_to="harvest") %>% data.frame()
df_harv_long<-dplyr::select(df_harv_long,-iter)
df_esc<-data.frame(cbind(df_scen,av_esc_diff))
df_esc_long<-df_esc %>% pivot_longer(!c(selectivity,mgmt,trends),names_to="iter",values_to="escapement") %>% data.frame()
df_esc_long<-dplyr::select(df_esc_long,-iter)
df_all<-data.frame(cbind(df_harv_long,escapement=df_esc_long$escapement))
##---------------------------------------------------------------## 'method'
df_all$method<-"TRM"
df_all$method[grepl("eq",df_all$mgmt)]<-"YPR" #"Yield-Per-Recruit"
df_all$method[grepl("dlm",df_all$mgmt)]<-"DLM" #"Dynamic Linear Model"
df_all<-dplyr::filter(df_all,method!="TRM")
df_all$method<-as.factor(as.character(df_all$method))
##-------------------------------------------------------## ordered factors
sel_levels<-c("6 inch gillnet", "unselective", "8.5 inch gillnet")
sel_labels<-c("small-mesh","unselective","large-mesh")
df_all$sel_ordered<-factor(df_all$selectivity,levels=sel_levels,labels=sel_labels)
# df_all<-dplyr::filter(df_all,sel_ordered!="unselective")
##----------------------------------------------------------## select trends
trends_used<-trend_names[c(4)]
df_all<-dplyr::filter(df_all,trends %in% trends_used)
##---------------------------------------------------------## select columns
df_all<-dplyr::select(df_all,-mgmt,-selectivity)
##-------------------------------------------------------------------## plot
p<-df_all %>% ggplot(aes(x=harvest,y=escapement,col=fct_inorder(method), group=fct_inorder(method)))+
  geom_hline(yintercept=0,linetype="solid",color="black",size=0.1)+
  geom_vline(xintercept=0,linetype="solid",color="black",size=0.1)+
  #geom_density_2d(size=0.25)+
  geom_point(size=1,shape=16,alpha=0.5)+
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  #coord_cartesian(xlim=c(-65,160),ylim=c(-65,160))+ 
  theme_classic()+
  labs(x="Mean harvest difference (%) relative to time-invariant model",
  y="Mean escapement difference (%)\nrelative to time-invariant model")+
  labs(color="")+
  # labs(color="Estimation method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text.x=element_text(angle=90,size=10,vjust=0.5,hjust=0.5),
        axis.text.y=element_text(size=10,angle=0,vjust=0.5,hjust=0.5),
        axis.title.x=element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        axis.title.y=element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        panel.border=element_rect(fill=NA,size=1),
        legend.background=element_rect(fill='transparent'),
        legend.position="none",
        # legend.position="top",
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(trends~sel_ordered,scales="fixed")+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10))
## ,legend.position=c(0.25,0.9)
ggsave("Performance-difference-escepement-vs-harvest-all.pdf",g,width=7.5, height=3,units="in")

##================================================## including density plots
use_selectivity<-"large-mesh"
#xlim<-c(-65,160);ylim<-c(-65,160)
df_one<-dplyr::filter(df_all,sel_ordered==use_selectivity)
main<-df_one %>% ggplot(aes(x=harvest,y=escapement,col=fct_inorder(method), group=fct_inorder(method)))+
  #annotate("text",x=100,y=200,size=3,label=use_selectivity)+
  geom_hline(yintercept=0,linetype="solid",color="black",size=0.1)+
  geom_vline(xintercept=0,linetype="solid",color="black",size=0.1)+
  #geom_density_2d(size=0.25)+
  geom_point(size=1,shape=16,alpha=0.5)+
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  # coord_cartesian(xlim=xlim,ylim=ylim)+
  ## don't use axes limits or else marginal axes don't align !!!
  theme_classic()+
  labs(x="Mean harvest difference (%)\nrelative to time-invariant model",
  y="Mean escapement difference (%)\nrelative to time-invariant model")+
  labs(color="")+
  # labs(color="Estimation method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text.x=element_text(angle=90,size=10,vjust=0.5,hjust=0.5),
        axis.text.y=element_text(size=10,angle=0,vjust=0.5,hjust=0.5),
        axis.title.x=element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        axis.title.y=element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        panel.border=element_rect(fill=NA,size=1),
        legend.background=element_rect(fill='transparent'),
        legend.position="none",
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        plot.margin=unit(c(1,1,1,1),"lines")
  )
g<-ggMarginal(main,groupColour=TRUE,groupFill=TRUE,margins="both",size=4, type="density")
ggsave(paste0("Performance-difference-and-density-",use_selectivity,".pdf"), g,width=4,height=4, units="in")

##========================================================================##
##=========================================================## plot quantiles
##========================================================================##
use_qunatiles<-c(0.05,0.1,0.5,0.9,0.95)
n_quants<-length(use_qunatiles)
for(q in 1:n_quants){
use_quant<-use_qunatiles[q]
##----------------------------------------------------## 'method'
df<-dplyr::select(scenarios,trends,selectivity,mgmt)
##----------------------------------------------------## 'method'
df$method<-"TRM"
df$method[grepl("eq",df$mgmt)]<-"YPR"
df$method[grepl("dlm",df$mgmt)]<-"DLM"
##----------------------------------------------------## median all metrics
df$av_harv<-apply(av_harv,1,function(x) quantile(x,prob=use_quant,na.rm=T))
df$av_harv<-df$av_harv/1e3
df$cv_harv<-apply(cv_harv,1,function(x) quantile(x,prob=use_quant,na.rm=T))
#df$av_esc<-apply(av_esc,1,function(x) quantile(x,prob=use_quant,na.rm=T))
#df$av_esc<-df$av_esc/1e3
df$av_ret<-apply(av_ret,1,function(x) quantile(x,prob=use_quant,na.rm=T))
df$av_ret<-df$av_ret/1e3
df$p_above_thresh<-apply(p_above_thresh,1,function(x) quantile(x,prob=use_quant,na.rm=T))
#df$p_no_harv<-apply(p_no_harv,1,function(x) quantile(x,prob=use_quant,na.rm=T))
##----------------------------------------------------------## metric labels
labs<-c("Mean harvest","Harvest stability","Mean spawner\nescapement","Probability\nabove threshold") #, "P no harvest")
n_metrics<-length(labs)
##==========================================================## long format
dfp<-df %>% pivot_longer(!c(trends,selectivity,mgmt,method), names_to="metric", values_to="median") %>% data.frame()
##-----------------------------------------------## add labels for metrics
metrics<-data.frame(name=colnames(df)[-c(1:4)])
metrics$label<-labs
for(i in 1:dim(dfp)[1]) {
  dfp$metric_label[i]<-metrics$label[dfp$metric[i]==metrics$name]
}
##-------------------------------------------------------## ordered factors
sel_levels<-c("6 inch gillnet", "unselective", "8.5 inch gillnet")
sel_labels<-c("small-mesh", "unselective", "large-mesh")
dfp$sel_ordered<-factor(dfp$selectivity,levels=sel_levels,labels=sel_labels)

##===================================================================## plot
p<-dfp %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point() + 
  scale_colour_manual(values=c("gray","dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="",y="Median value") + 
  labs(color="Accounting for trends")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(rows=vars(metric_label),cols=vars(sel_ordered),scale="free_y")+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10)) 
# ggsave(paste0("Performance metrics q=",use_quant,".pdf"),g,width=6.5,height=6.5,units="in")
} ## end loop over quantiles

##========================================================================##
##===========================================================## plot medians
##========================================================================##
df<-dplyr::select(scenarios,trends,selectivity,mgmt)
##----------------------------------------------------## 'method'
df$method<-"TRM"
df$method[grepl("eq",df$mgmt)]<-"YPR"
df$method[grepl("dlm",df$mgmt)]<-"DLM"
##----------------------------------------------------## median all metrics
df$av_harv<-apply(av_harv,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_harv<-df$av_harv/1e3
df$cv_harv<-apply(cv_harv,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_ret<-apply(av_ret,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_ret<-df$av_ret/1e3
# df$av_esc<-apply(av_esc,1,function(x) median(x,prob=use_quant,na.rm=T))
# df$av_esc<-df$av_esc/1e3
df$p_above_thresh<-apply(p_above_thresh,1,function(x) median(x,na.rm=T))
# df$p_no_harv<-apply(p_no_harv,1,function(x) median(x,na.rm=T))
##----------------------------------------------------------## metric labels
labs<-c("Mean harvest\n(thousands)","Harvest stability\n(1/CV)","Mean return\n(thousands)","Probability\nabove threshold") #, "P no harvest")
n_metrics<-length(labs)
##==========================================================## long format
dfp<-df %>% pivot_longer(!c(trends,selectivity,mgmt,method), names_to="metric", values_to="median") %>% data.frame()
##-----------------------------------------------## add labels for metrics
metrics<-data.frame(name=colnames(df)[-c(1:4)])
metrics$label<-labs
for(i in 1:dim(dfp)[1]) {
dfp$metric_label[i]<-metrics$label[dfp$metric[i]==metrics$name]
}
##----------------------------------------------------------## select trends
dfp<-dplyr::filter(dfp,trends!="age-length trends")
dfp$trends<-as.factor(as.character(dfp$trends))
##-------------------------------------------------------## ordered factors
sel_levels<-c("6 inch gillnet", "unselective", "8.5 inch gillnet")
sel_labels<-c("small-mesh", "unselective", "large-mesh")
dfp$sel_ordered<-factor(dfp$selectivity,levels=sel_levels,labels=sel_labels)

##=================================================## medians by target type
p<-dfp %>% 
  filter(method!="DLM") %>%
  filter(metric!="cv_harv") %>%
  ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point() + 
  scale_colour_manual(values=c("gray","dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="",y="") + 
  labs(color="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(rows=vars(metric_label),cols=vars(sel_ordered),scale="free_y",switch="y")+theme(strip.text.x=element_text(size=9),strip.text.y=element_text(size=9))+theme(strip.placement="outside")
ggsave(paste0("Performance_metrics.pdf"),g,width=5.2,height=5.5,units="in")

##========================================================================##
##========================================================================##
##============================================================## radar plots
##========================================================================##
##========================================================================##
mylabs<-c("Mean harvest","Harvest\nstability","Mean return"," Prob(>\n threshold)") 
cols_method<-c("dodgerblue2","goldenrod1","gray50") ## to macth order below
use_selectivities<-selectivity_names[c(3,2,1)]
n_used_sel<-length(use_selectivities)
# trends_used<-trend_names[c(1,3,4)]
trends_used<-trend_names[c(4)]
n_trend_used<-length(trends_used)
df_names<-data.frame(expand_grid(use_selectivities,trends_used))
df_names<-apply(df_names,1,function(x) paste0(x,collapse=" "))
##----------------------------------------## make data frames for plotting
cnt<-1
df_pl_all<-NA
df_list<-list()
for(j in 1:n_used_sel) {
  use_selectivity<-use_selectivities[j]
  for(i in 1:n_trend_used) {
    trend_name<-trends_used[i]
    df_pl<-dfp[dfp$trends==trend_name,]
    df_pl<-df_pl[df_pl$selectivity==use_selectivity,]
    df_pl<-df_pl %>% select(mgmt,median,metric_label)
    if(cnt==1) { df_pl_all<-df_pl } 
    else { df_pl_all<-data.frame(rbind(df_pl_all,df_pl)) }
    df_plw<-df_pl %>% pivot_wider(names_from=metric_label,values_from=median)
    df_list[[cnt]]<-df_plw
    cnt<-cnt+1
  } ## end i-loop
} ## end j-loop
names(df_list)<-df_names
##--------------------------## min/max for each metric across all scenarios
maxs<-df_pl_all %>% group_by(metric_label) %>% summarize(max=max(median)) %>% data.frame
maxs[,2]<-signif(1.01*maxs[,2],5)
# maxs$max[maxs$metric_label=="P > threshold"]<-1
mins<-df_pl_all %>% group_by(metric_label) %>% summarize(min=min(median)) %>% data.frame
mins[,2]<-signif(0.99*mins[,2],5)
mins$min<-0 ## center is actually zero
##---------------------------------------------------------## make gg plots
match_radarchart<-TRUE ## match order with radarchart()?
scale<-"within" ## scale radarplots "within" or "across" scenarios
n_used<-length(df_list)
cnt<-1
plts<-list()
for(j in 1:n_used_sel) {
  use_selectivity<-use_selectivities[j]  
  my_gg_labs<-mylabs
  if(match_radarchart) my_gg_labs<-my_gg_labs[c(1,4,3,2)] 
  for(i in 1:n_trend_used) {
    trend_name<-trends_used[i]
    df_use<-df_list[[cnt]]
    ##----------------## scale relative to min/max values for that scenario
    if(scale=="within") {
      df_range<-apply(df_use[,-1],2,function(x) (x-min(x))/max(x-min(x)))
      df_gg<-data.frame(cbind(df_use[,1],df_range)) 
    }
    ##----------------## scale relative to min/max values for all scenarios
    if(scale=="across") {
      df_gg<-df_use
      nms<-mins$metric_label
      for(k in 1:length(nms)){ 
      df_gg[,names(df_gg)==nms[k]]<-(df_gg[,names(df_gg)==nms[k]]- mins$min[mins$metric_label==nms[k]])/(maxs$max[maxs$metric_label==nms[k]]) }
    } ## end if statement
    ##------------------------------------------------## make ggradar plot
    if(match_radarchart) df_gg<-df_gg[,c(1,2,5,4,3)] 
    plt<-df_gg %>%
      ggradar(
        group.point.size=0,
        group.line.width=1,
        group.colours=cols_method,
        grid.label.size=0,
        grid.line.width=0.25,
        gridline.min.linetype="solid",
        gridline.mid.linetype="solid",
        gridline.max.linetype="solid",
        gridline.min.colour="gray",
        gridline.mid.colour="white",
        gridline.max.colour="gray",
        axis.line.colour="gray",
        axis.label.size=3,
        axis.labels=paste0(my_gg_labs),
        background.circle.colour="gray90",
        grid.min=0,  
        grid.mid=0.5,
        grid.max=1,  
        plot.title=paste0(use_selectivities[j]," - ",trends_used[i]),
        plot.extent.x.sf=1.25,
        plot.extent.y.sf=1.25,
        #legend.title="",
        #legend.text.size=3,
        plot.legend=FALSE
      )+
      theme(plot.title=element_text(size=10,face="bold",color="black"))
    plts[[cnt]]<-plt
    cnt<-cnt+1
  } ## end i-loop
} ## end j-loop
##-------------------------------------------------## arrange plots and save
plts_all<-arrangeGrob(grobs=plts,layout_matrix=matrix(1:n_used,ncol=n_used_sel,nrow=n_trend_used))
ggsave("Radarplots.pdf",plts_all,width=4*n_used_sel,height=3*n_trend_used,units="in")

##========================================================================##
##========================================================================##
##========================================================================##
##=====## 2nd set of scenarios on estimation methods and management startegy
##========================================================================##
##========================================================================##
##========================================================================##
scenarios<-scenarios_all[scenarios_all$factorMSY!=1,]
nscen<-dim(scenarios)[1]		
##---------------------------------------------------------## scenario names
scenario_names<-paste0(scenarios$trends," ",scenarios$mgmt," ", scenarios$factorMSY)
##--------------------------------------------------------## names of trends
trend_names<-unique(scenarios$trends)
n_trends<-length(trend_names)
##-------------------------------## names of estimation methods used in mgmt
mgmt_names<-unique(scenarios$mgmt)
n_mgmt<-length(mgmt_names)
##------------------------------------------## factors used in mgmt strategy
factorMSY_names<-unique(scenarios$factorMSY)
n_factors<-length(factorMSY_names)

##========================================================================##
##========================================================================##
##=========================================## differences in S_msy estimates
##========================================================================##
##========================================================================##
review_years<-seq((nyi+20),ny,goalfreq)
nrev<-length(review_years)
myarray<-array(NA,dim=c(nscen,niter,nrev))
S_msy_scen<-myarray
for(j in 1:nscen) { 
  for(k in 1:niter) { 
    s_msy_est<-as.vector(unlist(S_msy.list[[j]][[k]]))
    if(length(s_msy_est)!=0) S_msy_scen[j,k,]<-s_msy_est
  } ## end j-loop
} ## end i-loop
##------------------------------------------------------------## time period
## evaluated after 50 years, i.e. after the historical period
S_msy_hist<-S_msy_scen[,,which(review_years==nyh)]
##------------------------------------------------------## S_msy differences
## DLM and YPR analysis compared to traditional model
ref_msy_ind<-which(grepl("smsy_goal",scenario_names,fixed=T)) 
ypr_msy_ind<-which(grepl("s_eq_goal",scenario_names,fixed=T)) 
dlm_msy_ind<-which(grepl("smsy_dlm_goal",scenario_names,fixed=T))  
n_diffs<-length(ref_msy_ind) ## number of scenarios to compare
S_msy_diff_YPR<-S_msy_diff_DLM<-array(NA,dim=c(n_diffs,niter))
for(i in 1:n_diffs){
  ref_msys<-S_msy_hist[ref_msy_ind[i],]
  ypr_msys<-S_msy_hist[ypr_msy_ind[i],]
  dlm_msys<-S_msy_hist[dlm_msy_ind[i],]
  S_msy_diff_YPR[i,]<-100*(ypr_msys-ref_msys)/ref_msys
  S_msy_diff_DLM[i,]<-100*(dlm_msys-ref_msys)/ref_msys
}
##-----------------------------------------------## add scenario information
use_all<-dplyr::select(scenarios,trends,mgmt,factorMSY)
S_msy_df<-data.frame(cbind(use_all,S_msy_hist))
S_msy_df<-dplyr::filter(S_msy_df,trends!="continuing ASL trends")
S_msy_df<-dplyr::filter(S_msy_df,trends!="continuing trends")
use_ypr<-use_all[ypr_msy_ind,]
S_msy_df_YPR<-data.frame(cbind(use_ypr,S_msy_diff_YPR))
S_msy_df_YPR<-dplyr::filter(S_msy_df_YPR,trends!="continuing ASL trends")
S_msy_df_YPR<-dplyr::filter(S_msy_df_YPR,trends!="continuing trends")
use_dlm<-use_all[dlm_msy_ind,]
S_msy_df_DLM<-data.frame(cbind(use_dlm,S_msy_diff_DLM))
S_msy_df_DLM<-dplyr::filter(S_msy_df_DLM,trends!="continuing ASL trends")
S_msy_df_DLM<-dplyr::filter(S_msy_df_DLM,trends!="continuing trends")
n_diff<-dim(S_msy_df_DLM)[1]
##-------------------------------------------## add which alternative method
S_msy_df_YPR$method<-"YPR"
S_msy_df_DLM$method<-"DLM"
S_msy_diffs<-data.frame(rbind(S_msy_df_YPR,S_msy_df_DLM))
S_msy_diffs<-dplyr::select(S_msy_diffs,-mgmt)

##============================================================## select data
plot_smsy<-S_msy_diffs %>% pivot_longer(!c(trends,factorMSY,method), names_to="iteration", values_to="value") %>% data.frame()
plot_smsy<-plot_smsy[complete.cases(plot_smsy),] ## <1%
plot_smsy<-dplyr::select(plot_smsy,-iteration)

##================================================================## medians
S_msy_scen<-dplyr::select(S_msy_diffs,factorMSY,method,trends)
S_msy_iter<-dplyr::select(S_msy_diffs,-factorMSY,-method,-trends)
S_msy_med<-apply(S_msy_iter,1,function(x) median(x,na.rm=T))
S_msy_diff_med<-data.frame(cbind(S_msy_scen,median=S_msy_med))

##===========================================================## for plotting
method_colors<-c("dodgerblue2", "goldenrod1")
jig<-position_dodge(width=0.75)
##----------------------------------## function to plot median and quantiles
## using 90% CIs (5th-95th quantile)
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) }
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) }

##=========================================================## plot quantiles
method_colors<-c("dodgerblue2", "goldenrod1")
jig<-position_dodge(width=0.25)
colors<-rep(method_colors,n_diff)
ymin<-quantile(plot_smsy$value,probs=0.01,na.rm=T)
ymax<-quantile(plot_smsy$value,probs=0.99,na.rm=T)
p<-plot_smsy %>% ggplot(aes(x=fct_inorder(trends),y=value,fill=fct_inorder(method)))+
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  stat_summary(fun.data=summary_CI90,position=jig,size=0.2,color=colors)+
  stat_summary(fun.data=summary_CI50,position=jig,size=0.5,color=colors)+
  # scale_y_continuous(limits=c(ymin,ymax),breaks=seq(-50,150,25))+
  # ylim(-20,30)+
  ## Setting limits affects medians/CIs as it drops data from stat_summary
  ## Do not use 'ylim' or scale with 'limits' > only use coord_cartesian
  scale_y_continuous(breaks=seq(-100,150,25))+
  coord_cartesian(ylim=c(-55,80))+
  theme_classic()+
  labs(fill="Method")+
  labs(x="",y="Difference in S_MSY (%)\nfrom time-invariant model")+ 
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=10),
        panel.border=element_rect(fill=NA,size=1),
        legend.position="none",
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(cols=vars(factorMSY),scale="free")+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10))
ggsave("estimated_S_msy_by_method_quantiles_factorMSY.pdf",g,width=5,height=4,units="in")  

##===========================================================## plot medians
cols<-c("dodgerblue2","goldenrod1") ## when using fct_inorder(method)
ymin<-min(S_msy_diff_med$median,na.rm=T)*1.1
ymax<-max(S_msy_diff_med$median,na.rm=T)*1.1
p<-S_msy_diff_med %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=fct_inorder(method)))+
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=0.5) +
  #geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point(size=2.5) + 
  scale_y_continuous(breaks=seq(-100,150,10))+
  coord_cartesian(ylim=c(-20,30))+
  scale_colour_manual(values=cols)+
  theme_classic()+
  labs(x="",y="Difference in S_MSY (%)\nfrom time-invariant model")+ 
  labs(color="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=10),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(cols=vars(factorMSY),scale="free")+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10),legend.position=c(0.1,0.8))
ggsave("estimated_S_msy_by_method_medians_factorMSY.pdf",g,width=5,height=4,units="in")  

##========================================================================##
##====================## probability of Smsy greater than traditional method
##========================================================================##
S_msy_probs<-dplyr::select(S_msy_diffs,-factorMSY,-method,-trends)
S_msy_prob<-apply(S_msy_probs,1,function(x) length(x[which(x>0)])/length(x))
S_msy_scen<-dplyr::select(S_msy_diffs,factorMSY,method,trends)
S_msy_prob<-data.frame(cbind(S_msy_scen,S_msy_prob))
S_msy_prob<-dplyr::filter(S_msy_prob,trends!="continuing trends")

##===========================================================## plot probs
cols<-c("dodgerblue2","goldenrod1") ## when using fct_inorder(method)
p<-S_msy_prob %>% ggplot(aes(x=fct_inorder(trends),y=S_msy_prob,col=fct_inorder(method),group=fct_inorder(method)))+
  geom_hline(yintercept=0.5,linetype="solid",size=0.1)+ 
  geom_line(lwd=0.5) +
  #geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point(size=2.5) + 
  coord_cartesian(ylim=c(0,1))+
  scale_colour_manual(values=cols)+
  theme_classic()+
  labs(x="",y="Probability of estimate greater\nthan from time-invariant model")+ 
  labs(color="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=10),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(cols=vars(factorMSY),scale="free")+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10),legend.position=c(0.1,0.8)) 
ggsave("estimated_S_msy_by_method_probability_factorMSY.pdf",g,width=5,height=4,units="in")  

##========================================================================##
##========================================================================##
##=========================================## management performance metrics
##========================================================================##
##========================================================================##
## using simulated or observed data?
##---------------------------------## number of reconstructed observed years
ny_obs<-dim(data.frame(obs.list[[1]][[1]]))[1]
##--------------------------------------------------------## new trend names
trend_names_old<-trend_names
trend_names<-c("no trends","age-length trends","ASL trends stabilized","ASL trends continued" )
replace<-data.frame(from=trend_names_old,to=trend_names)
scenarios<-FindReplace(scenarios,Var="trends",replaceData=replace, from="from",to="to",exact=FALSE)
##------------------------------------------------------## array for results
scen_names<-paste0("scen=",seq(nscen))
iter_names<-paste0("iter=",seq(niter))
dnames<-list(scen_names,iter_names)
myarray<-array(dim=c(nscen,niter),dimnames=dnames)

##========================================================================##
##===========================================## performance trends over time
##========================================================================##
year_index_list<-list()
future_yrs<-seq(10,50,10)
n_index<-length(future_yrs)
thresh<-0.5
##------------------------------------------------------------## each period
#for(i in 1:n_index) { year_index_list[[i]]<-(nyh+1+(i-1)*10):(nyh+i*10) }
##-------------------------------------------------------------## cumulative
for(i in 1:n_index) { year_index_list[[i]]<-(nyh+1):(nyh+i*10) }
##==========================================## loop scenarios and iterations
av_harv_list<-av_ret_list<-av_esc_list<-p_above_thresh_list<-list()
av_harv_diff_list<-av_ret_diff_list<-av_esc_diff_list<-list()
for(i in 1:n_index) { 
  year_index<-year_index_list[[i]]
  av_harv<-av_ret<-av_esc<-p_above_thresh<-myarray  
  for(j in 1:nscen) {
    for(k in 1:niter) { 
      harv_jk<-data.frame(obs.list[[j]][[k]])$obsHarv[year_index] 
      if(is.null(harv_jk)){next}else{ av_harv[j,k]<-mean(harv_jk,na.rm=T)}
      esc_jk<-data.frame(obs.list[[j]][[k]])$obsEsc[year_index] 
      if(is.null(esc_jk)){next}else{ av_esc[j,k]<-mean(esc_jk,na.rm=T)}
      ret_jk<-data.frame(obs.list[[j]][[k]])$obsRet[year_index] 
      if(is.null(ret_jk)){next}else{ av_ret[j,k]<-mean(ret_jk,na.rm=T)}
      ##--------------------------------------## probability above threshold
      ricker_parms<-data.frame(sr_sim.list[[j]][[k]]) 
      max_rec<-(ricker_parms$alpha/ricker_parms$beta)*exp(-1)
      p_below_jk<-length(which(ret_jk<max_rec*thresh))/length(ret_jk)
      p_above_thresh[j,k]<-1-p_below_jk
    } ## end k-loop
  } ## end j-loop
  av_harv_list[[i]]<-av_harv
  av_esc_list[[i]]<-av_esc
  av_ret_list[[i]]<-av_ret
  p_above_thresh_list[[i]]<-p_above_thresh
  ## compute differences relative to traditional method as reference
  #print(paste0("done i=",i," part1"))
  av_harv_diff<-av_ret_diff<-av_esc_diff<-myarray  
  for(j in 1:nscen) {
    for(k in 1:niter) {
      if(j %in% ref_msy_ind) j_ref<-j
      if(j %in% ypr_msy_ind) j_ref<-j-4
      if(j %in% dlm_msy_ind) j_ref<-j-8
      av_harv_diff[j,k]<-(av_harv[j,k]-av_harv[j_ref,k])*100/av_harv[j_ref,k]
      av_esc_diff[j,k]<-(av_esc[j,k]-av_esc[j_ref,k])*100/av_esc[j_ref,k]
      av_ret_diff[j,k]<-(av_ret[j,k]-av_ret[j_ref,k])*100/av_ret[j_ref,k]
    } ## end k-loop
  } ## end j-loop
  av_harv_diff_list[[i]]<-av_harv_diff
  av_esc_diff_list[[i]]<-av_esc_diff
  av_ret_diff_list[[i]]<-av_ret_diff
  #print(paste0("done i=",i," part2"))
} ## end i-loop over year indices

##========================================================================##
##======================================================## conservation risk
##========================================================================##
df_risk_scen<-dplyr::select(scenarios,trends,factorMSY,mgmt)
##----------------------------------------------------## 'method'
df_risk_scen$method<-"TRM"
df_risk_scen$method[grepl("eq",df_risk_scen$mgmt)]<-"YPR"
df_risk_scen$method[grepl("dlm",df_risk_scen$mgmt)]<-"DLM"
##==================================================## prob above threshold
df_p_above_thresh<-array(dim=c(nscen,n_index))
for(i in 1:n_index) { df_p_above_thresh[,i]<-apply(p_above_thresh_list[[i]],1, function(x) median(x,na.rm=T)) }
colnames(df_p_above_thresh)<-future_yrs
##-------------------------------------------------## combine with scenarios
df_p_above_thresh<-data.frame(cbind(df_risk_scen,df_p_above_thresh))
##----------------------------------------------------------## select trends
df_p_above_thresh<-dplyr::filter(df_p_above_thresh,trends=="ASL trends continued")
df_p_above_thresh<-dplyr::select(df_p_above_thresh,-trends,-mgmt)
##-----------------------------------------------------------## pivot longer
df_pl<-df_p_above_thresh %>% pivot_longer(!c(factorMSY,method), names_to="period", values_to="median") %>% data.frame() 
df_pl$period<-as.numeric(gsub("X","",df_pl$period))
##-----------------------------------------------------------## metric label
df_pl$metric_label<-"Probability above threshold"
df_pl$MSYfactor<-factor(df_pl$factorMSY,levels=c(0.75,1.5),labels=c("aggressive","precautionary"))

##===================================================================## plot
p<-df_pl %>% ggplot(aes(x=period,y=median,col=fct_inorder(method))) +
  #geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  #geom_hline(yintercept=1,linetype="solid",size=0.1)+ 
  geom_line(lwd=1,aes(linetype=fct_inorder(MSYfactor))) +
  #geom_point(shape=1,size=2,fill=NA,color="black") +
  #geom_point() + 
  scale_colour_manual(values=c("gray","dodgerblue2","goldenrod1"))+
  #scale_y_continuous(expand=c(0.1,0.1)) +
  #scale_x_continuous(expand=c(0.1,0.1)) +
  #coord_cartesian(ylim=c(0,1))+
  theme_classic() +
  labs(x="Year",y="Median probability of\nabundance above threshold") + 
  labs(col="Method",linetype="Strategy")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        legend.key.size=unit(0.6,'cm'),
        panel.border=element_rect(fill=NA,size=1)
  )
ggsave("Future_prob_above_threshold_factorMSY.pdf",p,width=4.5,height=3,units="in")

##========================================================================##
##==========================## differences compared to time-invariant model
##========================================================================##

##==============================================================## scenarios
df_diff_scen<-dplyr::select(scenarios,trends,factorMSY,mgmt)
##----------------------------------------------------## 'method'
df_diff_scen$method<-"TRM"
df_diff_scen$method[grepl("eq",df_diff_scen$mgmt)]<-"YPR"
df_diff_scen$method[grepl("dlm",df_diff_scen$mgmt)]<-"DLM"

##================================================================## harvest
df_av_harv<-array(dim=c(nscen,n_index))
for(i in 1:n_index) { df_av_harv[,i]<-apply(av_harv_diff_list[[i]],1, function(x) median(x,na.rm=T)) }
colnames(df_av_harv)<-future_yrs
##-------------------------------------------------## combine with scenarios
df_av_harv<-data.frame(cbind(df_diff_scen,df_av_harv))
df_av_harv<-dplyr::filter(df_av_harv,method!="TRM")
##----------------------------------------------------------## select trends
df_av_harv<-dplyr::filter(df_av_harv,trends=="ASL trends continued")
df_av_harv<-dplyr::select(df_av_harv,-trends,-mgmt)
##-----------------------------------------------------------## pivot longer
df_av_harv_plot<-df_av_harv %>% pivot_longer(!c(factorMSY,method), names_to="period", values_to="median") %>% data.frame() 
df_av_harv_plot$period<-as.numeric(gsub("X","",df_av_harv_plot$period))
##-----------------------------------------------------------## metric label
df_av_harv_plot$metric_label<-"Mean harvest"

##=================================================================## return
df_av_ret<-array(dim=c(nscen,n_index))
for(i in 1:n_index) {df_av_ret[,i]<-apply(av_ret_diff_list[[i]],1, function(x) median(x,na.rm=T)) }
colnames(df_av_ret)<-future_yrs
##-------------------------------------------------## combine with scenarios
df_av_ret<-data.frame(cbind(df_diff_scen,df_av_ret))
df_av_ret<-dplyr::filter(df_av_ret,method!="TRM")
##----------------------------------------------------------## select trends
df_av_ret<-dplyr::filter(df_av_ret,trends=="ASL trends continued")
df_av_ret<-dplyr::select(df_av_ret,-trends,-mgmt)
##-----------------------------------------------------------## pivot longer
df_av_ret_plot<-df_av_ret %>% pivot_longer(!c(factorMSY,method), names_to="period", values_to="median") %>% data.frame() 
df_av_ret_plot$period<-as.numeric(gsub("X","",df_av_ret_plot$period))
##-----------------------------------------------------------## metric label
df_av_ret_plot$metric_label<-"Mean return"

##============================================================## escapement
df_av_esc<-array(dim=c(nscen,n_index))
for(i in 1:n_index) {df_av_esc[,i]<-apply(av_esc_diff_list[[i]],1, function(x) median(x,na.rm=T)) }
colnames(df_av_esc)<-future_yrs
##-------------------------------------------------## combine with scenarios
df_av_esc<-data.frame(cbind(df_diff_scen,df_av_esc))
df_av_esc<-dplyr::filter(df_av_esc,method!="None")
##----------------------------------------------------------## select trends
df_av_esc<-dplyr::filter(df_av_esc,trends=="ASL trends continued")
df_av_esc<-dplyr::select(df_av_esc,-trends,-mgmt)
##-----------------------------------------------------------## pivot longer
df_av_esc_plot<-df_av_esc %>% pivot_longer(!c(factorMSY,method), names_to="period", values_to="median") %>% data.frame() 
df_av_esc_plot$period<-as.numeric(gsub("X","",df_av_esc_plot$period))
##-----------------------------------------------------------## metric label
df_av_esc_plot$metric_label<-"Mean escapement"

##=================================================================## return
df_plot<-df_av_ret_plot
p<-df_plot %>% ggplot(aes(x=period,y=median,col=fct_inorder(method))) +
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=2) +
  #geom_point(shape=1,size=2,fill=NA,color="black") +
  #geom_point() + 
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  #coord_cartesian(xlim=c(5,55))+
  theme_classic() +
  labs(x="Year",y="Median % difference\nfrom time-invariant model") + 
  labs(fill="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        legend.key.size=unit(0.4,'cm'),
        panel.border=element_rect(fill=NA,size=1)
  )
g<-p+facet_grid(cols=vars(factorMSY),rows=vars(metric_label),scale="free") +theme(strip.text.x=element_text(size=9),strip.text.y=element_text(size=11), legend.position=c(0.1,0.8),legend.title=element_blank(),legend.text=element_text(size=9),legend.key=element_rect(colour=NA,fill=NA))
ggsave("Future_mean_return_differences_factorMSY.pdf",g,width=6,height=2.5,units="in")

##=================================================## harvest and escapement
df_plot<-data.frame(rbind(df_av_esc_plot,df_av_harv_plot))
p<-df_plot %>% ggplot(aes(x=period,y=median,col=fct_inorder(method))) +
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=2) +
  #geom_point(shape=1,size=2,fill=NA,color="black") +
  #geom_point() + 
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  #coord_cartesian(xlim=c(5,55))+
  theme_classic() +
  labs(x="Year",y="Median % difference from time-invariant model") + 
  labs(fill="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        legend.key.size=unit(0.4,'cm'),
        panel.border=element_rect(fill=NA,size=1)
  )
g<-p+facet_grid(cols=vars(factorMSY),rows=vars(metric_label),scale="free") +theme(strip.text.x=element_text(size=9),strip.text.y=element_text(size=11), legend.position=c(0.08,0.9),legend.title=element_blank(),legend.text=element_text(size=9),legend.key=element_rect(colour=NA,fill=NA))
ggsave("Future_mean_harvest_escapement_differences_factorMSY.pdf",g,width=6,height=4,units="in")

##========================================================================##
##==========================================## performance for last 50 years
##========================================================================##
##--------------------------------------------------------## years evaluated
# year_index<-1:ny_obs ## all reconstructed years
year_index<-(nyh+1):ny_obs ## most recent e.g. 50 years

##=============================================## for each estimation method
sum_log_harv<-cv_harv<-av_harv<-av_log_harv<-p_no_harv<-myarray ## harvest 
sum_log_esc<-cv_esc<-av_esc<-av_log_esc<-myarray ## escapement 
sum_log_ret<-cv_ret<-av_ret<-av_log_ret<-myarray ## total return 
p_below_thresh<-p_above_thresh<-myarray ## conservation threshold
thresh<-0.5 ## proportion of maximum recruitment at equilibrium

##------------------------------------------## loop scenarios and iterations
for(j in 1:nscen) {
  for(k in 1:niter) { 
    ##---------------------------------------------------------## escapement
    # esc_jk<-data.frame(data.list[[j]][[k]])$Esc[year_index] ## simulated 
    esc_jk<-data.frame(obs.list[[j]][[k]])$obsEsc[year_index] ## observed 
    if(is.null(esc_jk)) { next } else { 
      sum_log_esc[j,k]<-sum(log(esc_jk+0.0001))
      cv_esc[j,k]<-1/(sd(esc_jk,na.rm=T)/mean(esc_jk,na.rm=T)) ## 1/CV
      av_esc[j,k]<-mean(esc_jk,na.rm=T)
      av_log_esc[j,k]<-mean(log(esc_jk),na.rm=T)
    }
    ##----------------------------------------------------------## harvest
    # harv_jk<-data.frame(data.list[[j]][[k]])$Harv[year_index] ## simulated
    harv_jk<-data.frame(obs.list[[j]][[k]])$obsHarv[year_index] ## observed 
    if(is.null(harv_jk)) { next } else { 
      sum_log_harv[j,k]<-sum(log(harv_jk+0.0001))
      cv_harv[j,k]<-1/(sd(harv_jk,na.rm=T)/mean(harv_jk,na.rm=T)) ## 1/CV
      av_harv[j,k]<-mean(harv_jk,na.rm=T)
      av_log_harv[j,k]<-mean(log(harv_jk),na.rm=T)
      p_no_harv[j,k]<-length(which(harv_jk==0))/length(harv_jk)
    }
    ##----------------------------------------------------------## return
    # ret_jk<-data.frame(data.list[[j]][[k]])$Ret[year_index] ## simulated 
    ret_jk<-data.frame(obs.list[[j]][[k]])$obsRet[year_index] ## observed 
    if(is.null(ret_jk)) { next } else { 
      sum_log_ret[j,k]<-sum(log(ret_jk+0.0001))
      cv_ret[j,k]<-1/(sd(ret_jk,na.rm=T)/mean(ret_jk,na.rm=T)) ## 1/CV
      av_ret[j,k]<-mean(ret_jk,na.rm=T)
      av_log_ret[j,k]<-mean(log(ret_jk),na.rm=T)
      ## maximum recruitment given ricker parameters
      ricker_parms<-data.frame(sr_sim.list[[j]][[k]]) 
      max_rec<-(ricker_parms$alpha/ricker_parms$beta)*exp(-1) ## Hilborn ~
      ## probability above/below threshold
      p_below_jk<-length(which(ret_jk<max_rec*thresh))/length(ret_jk)
      p_below_thresh[j,k]<-p_below_jk 
      p_above_thresh[j,k]<-1-p_below_jk
    }
  } ## end k-loop
} ## end j-loop
p_below_thresh[is.na(p_below_thresh)]<-0

##=================================================## ratios and differences
av_harv_diff<-cv_harv_diff<-cv_harv_ratio<-av_log_harv_diff<-myarray
av_esc_diff<-av_log_esc_diff<-myarray
av_ret_diff<-av_log_ret_diff<-p_above_thresh_ratio<-p_no_harv_ratio<-myarray
##-----------------------------------------------------------## get indices
ref_msy_ind<-which(grepl("smsy_goal",scenario_names,fixed=T)) 
ypr_msy_ind<-which(grepl("s_eq_goal",scenario_names,fixed=T)) 
dlm_msy_ind<-which(grepl("smsy_dlm_goal",scenario_names,fixed=T))  
##---------------------------------------------------## compute differences
for(j in 1:nscen) {
  for(k in 1:niter) {
    if(j %in% ref_msy_ind) j_ref<-j
    if(j %in% ypr_msy_ind) j_ref<-j-4
    if(j %in% dlm_msy_ind) j_ref<-j-8
    cv_harv_diff[j,k]<-(cv_harv[j,k]-cv_harv[j_ref,k])*100/cv_harv[j_ref,k]
    cv_harv_ratio[j,k]<-cv_harv[j,k]/cv_harv[j_ref,k]
    av_harv_diff[j,k]<-(av_harv[j,k]-av_harv[j_ref,k])*100/av_harv[j_ref,k]
    av_log_harv_diff[j,k]<-(av_log_harv[j,k]-av_log_harv[j_ref,k])
    av_esc_diff[j,k]<-(av_esc[j,k]-av_esc[j_ref,k])*100/av_esc[j_ref,k]
    av_log_esc_diff[j,k]<-(av_log_esc[j,k]-av_log_esc[j_ref,k])
    av_ret_diff[j,k]<-(av_ret[j,k]-av_ret[j_ref,k])*100/av_ret[j_ref,k]
    av_log_ret_diff[j,k]<-(av_log_ret[j,k]-av_log_ret[j_ref,k])
    p_above_thresh_ratio[j,k]<-p_above_thresh[j,k]/p_above_thresh[j_ref,k]
    p_no_harv_ratio[j,k]<-p_no_harv[j,k]/p_no_harv[j_ref,k]
  } ## end k-loop
} ## end j-loop

##========================================================================##
##=================================## plot conservation risk for all methods
##========================================================================##
df_risk_scen<-dplyr::select(scenarios,trends,factorMSY,mgmt)
##----------------------------------------------------## 'method'
df_risk_scen$method<-"TRM"
df_risk_scen$method[grepl("eq",df_risk_scen$mgmt)]<-"YPR"
df_risk_scen$method[grepl("dlm",df_risk_scen$mgmt)]<-"DLM"
##----------------------------------------------------------------## metric
df_risk_scen$p_above_thresh<-apply(p_above_thresh,1,function(x) median(x,na.rm=T))
##-----------------------------------------------------------## metric label
#df_risk_scen$metric_label<-"Probability above threshold"
df_risk_scen$MSYfactor<-factor(df_risk_scen$factorMSY,levels=c(0.75,1.5),labels=c("aggressive","conservative"))
##----------------------------------------------------------## long format
df_risk_plot<-df_risk_scen %>% pivot_longer(!c(trends,factorMSY,mgmt,method,MSYfactor),names_to="metric",values_to="median") %>% data.frame()
##----------------------------------------------------------## select trends
df_risk_plot<-dplyr::filter(df_risk_plot,trends!="age-length trends")
df_risk_plot$trends<-as.factor(as.character(df_risk_plot$trends))
# df_risk_plot$MSYfactor<-factor(df_risk_plot$factorMSY,levels=c(0.75,1.5), labels=c("aggressive","conservative"))
df_risk_plot$factorMSY<-as.factor(df_risk_plot$factorMSY)

##===================================================================## plot
p<-df_risk_plot %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=fct_inorder(method)))+
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=0.5,aes(linetype=factorMSY)) +
  geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point() + 
  scale_colour_manual(values=c("gray","dodgerblue2","goldenrod1"))+
  theme_classic() +
  labs(x="",y="Median probability of\nabundance above threshold") + 
  labs(col="Method",linetype="Strategy")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        legend.key.size=unit(0.6,'cm'),
        panel.border=element_rect(fill=NA,size=1)
  )
g<-p+facet_grid(rows=vars(factorMSY),scale="free") +theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10),legend.position=c(0.08,0.94))
ggsave("Metrics_prob_above_threshold_factorMSY.pdf",g,width=3,height=6,units="in")

##========================================================================##
##======================## plot metrics that compare to time-invariant model
##========================================================================##
df_diff<-dplyr::select(scenarios,trends,factorMSY,mgmt)
##----------------------------------------------------## 'method'
df_diff$method<-"TRM"
df_diff$method[grepl("eq",df_diff$mgmt)]<-"YPR"
df_diff$method[grepl("dlm",df_diff$mgmt)]<-"DLM"
##----------------------------------------------------## median all metrics
df_diff$av_harv<-apply(av_harv_diff,1,function(x) median(x,na.rm=T))
df_diff$av_esc<-apply(av_esc_diff,1,function(x) median(x,na.rm=T))
df_diff$av_ret<-apply(av_ret_diff,1,function(x) median(x,na.rm=T))
df_diff$p_harv_above_ref<-apply(av_harv_diff,1,function(x) length(x[which(x>0)])/length(x))
df_diff$p_esc_above_ref<-apply(av_esc_diff,1,function(x) length(x[which(x>0)])/length(x))
df_diff$p_ret_above_ref<-apply(av_ret_diff,1,function(x) length(x[which(x>0)])/length(x))
##----------------------------------------------------------## metric labels
labs<-c("Mean harvest","Mean escapement","Mean return","P(harv>harv_ref)","P(esc>esc_ref)","P(ret>ret_ref)")
##==========================================================## long format
df_diff_plot<-df_diff %>% pivot_longer(!c(trends,factorMSY,mgmt,method), names_to="metric", values_to="median") %>% data.frame()
##-----------------------------------------------## add labels for metrics
metrics<-data.frame(name=colnames(df_diff)[-c(1:4)])
metrics$label<-labs
for(l in 1:dim(df_diff_plot)[1]) { df_diff_plot$metric_label[l]<-metrics$label[df_diff_plot$metric[l]==metrics$name] }
##-------------------------------------------------------## drop references
## values are differences compared to 'traditional' MSY model 
df_diff_plot<-df_diff_plot[df_diff_plot$mgmt!="smsy_goal",]
df_diff_plot<-df_diff_plot[df_diff_plot$mgmt!="umsy_goal",]
##----------------------------------------------------------## select trends
df_diff_plot<-dplyr::filter(df_diff_plot,trends!="age-length trends")
df_diff_plot$trends<-as.factor(as.character(df_diff_plot$trends))
##-----------------------------------------## split ratios and % differences
df_diff_plot1<-df_diff_plot[df_diff_plot$metric %in% c("av_harv","av_esc"),]
df_diff_plot2<-df_diff_plot[df_diff_plot$metric %in% c("av_ret"),]
df_diff_plot3<-df_diff_plot[df_diff_plot$metric %in% c("p_ret_above_ref"),]

##=================================================## harvest and escapement
p<-df_diff_plot1 %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point() + 
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="",y="Median % difference from time-invariant model") + 
  labs(color="")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank()
  )
g<-p+facet_grid(rows=vars(metric_label),cols=vars(factorMSY),scale="free") +theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10),legend.position=c(0.08,0.94))
ggsave("Metrics_differences_harv_esc_factorMSY.pdf",g,width=5,height=5.5,units="in")

##=================================================================## return
p<-df_diff_plot2 %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_hline(yintercept=1,linetype="solid",size=0.1)+ 
  geom_line(lwd=1) +
  geom_point(size=2.5) + 
  geom_point(shape=1,size=3,fill=NA,color="black") +
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="",y="Median % difference\nfrom time-invariant model") + 
  labs(color="")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank()
  )
g<-p+facet_grid(cols=vars(factorMSY),scale="free") +theme(strip.text.x=element_text(size=10),legend.position=c(0.08,0.94))
ggsave("Metrics_differences_return_factorMSY.pdf",g,width=5,height=4,units="in")

##================================================================## P(>ref)
p<-df_diff_plot3 %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_hline(yintercept=0.5,linetype="solid",size=0.1)+ 
  geom_line(lwd=1) +
  geom_point(size=2.5) + 
  geom_point(shape=1,size=3,fill=NA,color="black") +
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="",y="Probability of mean return larger\nthan for time-invariant model") + 
  labs(color="")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank()
  )
g<-p+facet_grid(cols=vars(factorMSY),scale="free") +theme(strip.text.x=element_text(size=10),legend.position=c(0.08,0.94))
ggsave("Metrics_P_return_higher_ref_factorMSY.pdf",g,width=5,height=4,units="in")

##========================================================================##
##======================## differences all iterations harvest vs escapement
##========================================================================##
df_scen<-dplyr::select(scenarios,trends,factorMSY,mgmt)
df_harv<-data.frame(cbind(df_scen,av_harv_diff))
df_harv_long<-df_harv %>% pivot_longer(!c(factorMSY,mgmt,trends),names_to="iter",values_to="harvest") %>% data.frame()
df_harv_long<-dplyr::select(df_harv_long,-iter)
df_esc<-data.frame(cbind(df_scen,av_esc_diff))
df_esc_long<-df_esc %>% pivot_longer(!c(factorMSY,mgmt,trends),names_to="iter",values_to="escapement") %>% data.frame()
df_esc_long<-dplyr::select(df_esc_long,-iter)
df_all<-data.frame(cbind(df_harv_long,escapement=df_esc_long$escapement))
##---------------------------------------------------------------## 'method'
df_all$method<-"Time-invariant model"
df_all$method[grepl("eq",df_all$mgmt)]<-"YPR" #"Yield-Per-Recruit"
df_all$method[grepl("dlm",df_all$mgmt)]<-"DLM" #"Dynamic Linear Model"
df_all<-dplyr::filter(df_all,method!="Time-invariant model")
df_all$method<-as.factor(as.character(df_all$method))
##----------------------------------------------------------## select trends
trends_used<-trend_names[c(4)]
df_all<-dplyr::filter(df_all,trends %in% trends_used)
##---------------------------------------------------------## select columns
df_all<-dplyr::select(df_all,-mgmt)
##-------------------------------------------------------------------## plot
p<-df_all %>% ggplot(aes(x=harvest,y=escapement,col=fct_inorder(method), group=fct_inorder(method)))+
  geom_hline(yintercept=0,linetype="solid",color="black",size=0.1)+
  geom_vline(xintercept=0,linetype="solid",color="black",size=0.1)+
  #geom_density_2d(size=0.25)+
  geom_point(size=1,shape=16,alpha=0.5)+
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  coord_cartesian(xlim=c(-110,160),ylim=c(-110,160))+
  theme_classic()+
  labs(x="Mean harvest difference (%) relative to time-invariant model",
       y="Mean escapement difference (%)\nrelative to time-invariant model")+
  labs(color="")+
  # labs(color="Estimation method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text.x=element_text(angle=90,size=10,vjust=0.5,hjust=0.5),
        axis.text.y=element_text(size=10,angle=0,vjust=0.5,hjust=0.5),
        axis.title.x=element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        axis.title.y=element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        panel.border=element_rect(fill=NA,size=1),
        legend.background=element_rect(fill='transparent'),
        legend.position="none",
        # legend.position="top",
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(trends~factorMSY,scales="fixed")+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10))
## ,legend.position=c(0.25,0.9)
ggsave("Performance-difference-escepement-vs-harvest-all_factorMSY.pdf",g,width=7.5, height=3,units="in")

##================================================## including density plots
use_factorMSY<-1.5
#xlim<-c(-110,160);ylim<-c(-110,160)
xmax<-max(df_one$harvest,na.rm=T);ymax<-max(df_one$escapement,na.rm=T)
df_one<-dplyr::filter(df_all,factorMSY==use_factorMSY)
main<-df_one %>% ggplot(aes(x=harvest,y=escapement,col=fct_inorder(method), group=fct_inorder(method)))+
  # annotate("text",x=xlim[2]-35,y=ylim[2]-5,size=3,label=use_factorMSY)+
  annotate("text",x=xmax-10,y=ymax-10,size=3,label=use_factorMSY)+
  geom_hline(yintercept=0,linetype="solid",color="black",size=0.1)+
  geom_vline(xintercept=0,linetype="solid",color="black",size=0.1)+
  #geom_density_2d(size=0.25)+
  geom_point(size=1,shape=16,alpha=0.5)+
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  #coord_cartesian(xlim=xlim,ylim=ylim)+
  theme_classic()+
  labs(x="Mean harvest difference (%)\nrelative to time-invariant model",
       y="Mean escapement difference (%)\nrelative to time-invariant model")+
  labs(color="")+
  # labs(color="Estimation method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text.x=element_text(angle=90,size=10,vjust=0.5,hjust=0.5),
        axis.text.y=element_text(size=10,angle=0,vjust=0.5,hjust=0.5),
        axis.title.x=element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        axis.title.y=element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        panel.border=element_rect(fill=NA,size=1),
        legend.background=element_rect(fill='transparent'),
        legend.position="none",
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        plot.margin=unit(c(1,1,1,1),"lines")
  )
g<-ggMarginal(main,groupColour=TRUE,groupFill=TRUE,margins="both",size=4, type="density")
ggsave(paste0("Performance-difference-and-density-",use_factorMSY,"_factorMSY.pdf"), g,width=4,height=4, units="in")

##========================================================================##
##===========================================================## plot medians
##========================================================================##
df<-dplyr::select(scenarios,trends,factorMSY,mgmt)
##----------------------------------------------------## 'method'
df$method<-"TRM"
df$method[grepl("eq",df$mgmt)]<-"YPR"
df$method[grepl("dlm",df$mgmt)]<-"DLM"
##----------------------------------------------------## median all metrics
df$av_harv<-apply(av_harv,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_harv<-df$av_harv/1e3
df$cv_harv<-apply(cv_harv,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_ret<-apply(av_ret,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_ret<-df$av_ret/1e3
#df$av_esc<-apply(av_esc,1,function(x) median(x,prob=use_quant,na.rm=T))
df$p_above_thresh<-apply(p_above_thresh,1,function(x) median(x,na.rm=T))
# df$p_no_harv<-apply(p_no_harv,1,function(x) median(x,na.rm=T))
##----------------------------------------------------------## metric labels
labs<-c("Mean harvest\n(thousands)","Harvest stability\n(1/CV)","Mean return\n(thousands)","Probability\nabove threshold")
n_metrics<-length(labs)
##==========================================================## long format
#if(n_select>1) dfp<-df %>% pivot_longer(!c(trends,selectivity,mgmt,method), names_to="metric", values_to="median") %>% data.frame()
if(n_factors>1) dfp<-df %>% pivot_longer(!c(trends,factorMSY,mgmt,method), names_to="metric", values_to="median") %>% data.frame()
##-----------------------------------------------## add labels for metrics
metrics<-data.frame(name=colnames(df)[-c(1:4)])
metrics$label<-labs
for(i in 1:dim(dfp)[1]) {
  dfp$metric_label[i]<-metrics$label[dfp$metric[i]==metrics$name]
}
##----------------------------------------------------------## select trends
dfp<-dplyr::filter(dfp,trends!="age-length trends")
dfp$trends<-as.factor(as.character(dfp$trends))

##-------------------------------------------------------## ordered factors
#fac_levels<-c("0.75","1.5")
#fac_labels<-c("0.75*SMSY","1.5*SMSY")
#fac_labels<-c(expression(0.75*S[MSY]),expression(1.5*S[MSY]))
#dfp$factorMSY<-factor(dfp$factorMSY,levels=fac_levels,labels=fac_labels)

##=================================================## medians by target type
p<-dfp %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point() + 
  scale_colour_manual(values=c("gray","dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="",y="") + 
  labs(color="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8)
  )
# col.labs<-c(expression(0.75*S[MSY]),expression(1.5*S[MSY]))
#col.labs<-c(bquote(0.75*S[MSY]),bquote(1.5*S[MSY]))
#names(col.labs)<-c("0.75","1.5")
g<-p+facet_grid(rows=vars(metric_label),cols=vars(factorMSY),scale="free_y", switch="y")+theme(strip.text.x=element_text(size=9),strip.text.y=element_text(size=9))+theme(strip.placement="outside")
## ,labeller=labeller(factorMSY=col.labs)
ggsave(paste0("Performance_metrics_factorMSY.pdf"),g,width=4.2,height=6.5,units="in")

##=================================================## medians by target type
dfp_few<-dplyr::filter(dfp,metric %in% c("av_ret","p_above_thresh"))
p<-dfp_few %>% ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2,fill=NA,color="black") +
  geom_point() + 
  scale_colour_manual(values=c("gray","dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="",y="Median value") + 
  labs(color="Accounting\nfor trends")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(rows=vars(metric_label),cols=vars(factorMSY),scale="free_y")+theme(strip.text.x=element_text(size=9),strip.text.y=element_text(size=9))
ggsave(paste0("Performance_metrics_few_factorMSY.pdf"),g,width=4,height=4,units="in")

##========================================================================##
##========================================================================##
##============================================================## radar plots
##========================================================================##
##========================================================================##
mylabs<-c("Mean harvest","Harvest\nstability","Mean return"," Prob(>\n threshold)") 
cols_method<-c("dodgerblue2","goldenrod1","gray50") ## to macth order below
use_factorMSYs<-factorMSY_names
n_used_fac<-length(use_factorMSYs)
# trends_used<-trend_names[c(1,3,4)]
trends_used<-trend_names[c(4)]
n_trend_used<-length(trends_used)
df_names<-data.frame(expand_grid(use_factorMSYs,trends_used))
df_names<-apply(df_names,1,function(x) paste0(x,collapse=" "))
##----------------------------------------## make data frames for plotting
cnt<-1
df_pl_all<-NA
df_list<-list()
for(j in 1:n_used_fac) {
  use_factorMSY<-use_factorMSYs[j]
  for(i in 1:n_trend_used) {
    trend_name<-trends_used[i]
    df_pl<-dfp[dfp$trends==trend_name,]
    df_pl<-df_pl[df_pl$factorMSY==use_factorMSY,]
    df_pl<-df_pl %>% select(mgmt,median,metric_label)
    if(cnt==1) { df_pl_all<-df_pl } 
    else { df_pl_all<-data.frame(rbind(df_pl_all,df_pl)) }
    df_plw<-df_pl %>% pivot_wider(names_from=metric_label,values_from=median)
    df_list[[cnt]]<-df_plw
    cnt<-cnt+1
  } ## end i-loop
} ## end j-loop
names(df_list)<-df_names
##--------------------------## min/max for each metric across all scenarios
maxs<-df_pl_all %>% group_by(metric_label) %>% summarize(max=max(median)) %>% data.frame
maxs[,2]<-signif(1.01*maxs[,2],5)
# maxs$max[maxs$metric_label=="P > threshold"]<-1
mins<-df_pl_all %>% group_by(metric_label) %>% summarize(min=min(median)) %>% data.frame
mins[,2]<-signif(0.99*mins[,2],5)
mins$min<-0 ## center is actually zero
##---------------------------------------------------------## make gg plots
match_radarchart<-TRUE ## match order with radarchart()?
scale<-"across" ## scale radarplots "within" or "across" scenarios
n_used<-length(df_list)
cnt<-1
plts<-list()
for(j in 1:n_used_fac) {
  use_factorMSY<-use_factorMSYs[j]  
  my_gg_labs<-mylabs
  if(match_radarchart) my_gg_labs<-my_gg_labs[c(1,4,3,2)] 
  for(i in 1:n_trend_used) {
    trend_name<-trends_used[i]
    df_use<-df_list[[cnt]]
    ##----------------## scale relative to min/max values for that scenario
    if(scale=="within") {
      df_range<-apply(df_use[,-1],2,function(x) (x-min(x))/max(x-min(x)))
      df_gg<-data.frame(cbind(df_use[,1],df_range)) 
    }
    ##----------------## scale relative to min/max values for all scenarios
    if(scale=="across") {
      df_gg<-df_use
      nms<-mins$metric_label
      for(k in 1:length(nms)){ 
        df_gg[,names(df_gg)==nms[k]]<-(df_gg[,names(df_gg)==nms[k]]- mins$min[mins$metric_label==nms[k]])/(maxs$max[maxs$metric_label==nms[k]]) }
    } ## end if statement
    ##------------------------------------------------## make ggradar plot
    if(match_radarchart) df_gg<-df_gg[,c(1,2,5,4,3)] 
    plt<-df_gg %>%
      ggradar(
        group.point.size=0,
        group.line.width=1,
        group.colours=cols_method,
        grid.label.size=0,
        grid.line.width=0.25,
        gridline.min.linetype="solid",
        gridline.mid.linetype="solid",
        gridline.max.linetype="solid",
        gridline.min.colour="gray",
        gridline.mid.colour="white",
        gridline.max.colour="gray",
        axis.line.colour="gray",
        axis.label.size=3,
        axis.labels=paste0(my_gg_labs),
        background.circle.colour="gray90",
        grid.min=0,  
        grid.mid=0.5,
        grid.max=1,  
        plot.title=paste0(use_factorMSYs[j]," - ",trends_used[i]),
        plot.extent.x.sf=1.25,
        plot.extent.y.sf=1.25,
        #legend.title="",
        #legend.text.size=3,
        plot.legend=FALSE
      )+
      theme(plot.title=element_text(size=10,face="bold",color="black"))
    plts[[cnt]]<-plt
    cnt<-cnt+1
  } ## end i-loop
} ## end j-loop
##-------------------------------------------------## arrange plots and save
plts_all<-arrangeGrob(grobs=plts,layout_matrix=matrix(1:n_used,ncol=n_used_fac,nrow=n_trend_used))
ggsave("Radarplots_factorMSY.pdf",plts_all,width=4*n_used_fac,height=3*n_trend_used,units="in")

##========================================================================##
##========================================================================##
##============================================================## end of code
##========================================================================##
##========================================================================##