##========================================================================##
##                                                                        ##
##                 Plot simulation model results for MSE                  ##
##                                                                        ##
##========================================================================##
##===============================================================## packages
pkgs<-c("here","tidyverse","dplyr","scales","RColorBrewer","readxl","gsl","vioplot", "caroline","ggplot2","viridis","forcats","gridExtra","lemon","fmsb","MetBrewer","ggradar","showtext","DataCombine","ggridges","ggExtra","remotes","ggpubr")
if(length(setdiff(pkgs,rownames(installed.packages())))>0) {install.packages(setdiff(pkgs,rownames(installed.packages())),dependencies=TRUE)}
invisible(lapply(pkgs,library,character.only=T))

##========================================================================##
##========================================================================##
##===========================================================## load results
##========================================================================##
##========================================================================##
path<-paste0(here::here(),"/R/out/")
##=============================================================## simulation
myfiles<-list.files(path,pattern="\\.Rdata$")
myfile<-myfiles[grep("parms",myfiles,invert=TRUE)]
load(paste0(here::here(),"/R/out/",myfile))
timestamp<-substr(myfile,5,nchar(myfile)-6)
##=============================================================## parameters
parameters<-readRDS(paste0(path,"run_",timestamp,"_parms.Rdata"))
npara<-length(parameters)
for(i in 1:npara) { assign(names(parameters)[i],unlist(parameters[i])) }
##==============================================================## scenarios
scenarios<-read_excel(paste0(path,"run_",timestamp,"_scen.xlsx"))[,-1]
scenarios_all<-data.frame(scenarios)
##==============================================================## directory
dir.create(file.path(paste0(path,timestamp)),showWarnings=F)
setwd(file.path(paste0(path,timestamp)))

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
##========================## Figure 2: trends in mean fecundity and egg mass
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
s_labels<-c("small-mesh","unselective","large-mesh")
plot_egg_mass$sel_ordered<-factor(plot_egg_mass$selectivity,levels=s_levels,labels=s_labels)
plot_egg_mass<-dplyr::select(plot_egg_mass,-iteration,-mgmt,-selectivity)
##----------------------------------## function to plot median and quantiles
## using 90% CIs (5th-95th quantile)
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) }
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) }
##-----------------------------------------------------------------## ggplot
jig<-position_dodge(width=0.6)
p<-plot_egg_mass %>% 
  ggplot(aes(x=fct_inorder(trends),y=value,fill=sel_ordered))+
  geom_violin(position=jig,trim=F,lwd=0.1,scale="width",width=0.6)+ 
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
        panel.border=element_rect(fill=NA,linewidth=1),
        legend.key.size=unit(0.4,'cm'),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        legend.position=c(0.82,0.86)
  )
p<-p+guides(fill=guide_legend(override.aes=list(size=0.3)))
ggsave("Figure2.pdf",p,width=4.2,height=5.5,units="in")

##========================================================================##
##========================================================================##
##===============================## Figure 3: differences in S_msy estimates
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
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) }
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) }

##===========================================================## get legend
method_colors<-c("dodgerblue2", "goldenrod1")
ymin<-min(S_msy_diff_med$median,na.rm=T)*1.1
ymax<-max(S_msy_diff_med$median,na.rm=T)*1.1
p<-S_msy_diff_med %>% 
  filter(method %in% c("YPR","DLM")) %>%
  ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method)))+
  scale_colour_manual(values=method_colors)+
  geom_point(size=2.5) + 
  theme_classic()+
  labs(color="Method")+
  theme(legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        legend.position=c(0.3,0.65))

lgnd<-get_legend(p,position=NULL)
as_ggplot(lgnd)

##=========================================================## plot quantiles
jig<-position_dodge(width=0.5)
colors<-rep(method_colors,n_diff)
ymin<-quantile(plot_smsy$value,probs=0.01,na.rm=T)
ymax<-quantile(plot_smsy$value,probs=0.99,na.rm=T)
p<-plot_smsy %>% ggplot(aes(x=fct_inorder(trends),y=value,fill=fct_inorder(method)))+
  geom_hline(yintercept=0,linetype="solid",linewidth=0.1)+ 
  stat_summary(fun.data=summary_CI90,position=jig,linewidth=0.2,color=colors)+
  stat_summary(fun.data=summary_CI50,position=jig,linewidth=0.6,color=colors)+
  scale_y_continuous(breaks=seq(-100,150,25))+
  coord_cartesian(ylim=c(-55,80))+
  theme_classic()+
  labs(fill="Method")+ 
  labs(x="",y=expression("Difference in "*S[MSY]*" (%)"))+ 
  theme(strip.background=element_blank(),
        axis.line=element_line(linewidth=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=10),
        panel.border=element_rect(fill=NA,linewidth=1),
        legend.position="none")
col.labs<-c("large-mesh","unselective","small-mesh")
names(col.labs)<-c("8.5 inch gillnet","unselective","6 inch gillnet")
g<-p+
  facet_grid(cols=vars(sel_ordered),labeller=labeller(sel_ordered=col.labs))+
  theme(strip.text.x=element_text(size=10))
g<-cowplot::plot_grid(g,lgnd,ncol=2,rel_widths=c(0.8,0.2)) ## add legend
ggsave("Figure3.pdf",g,width=6,height=4,units="in")  

##========================================================================##
##========================================================================##
##===============================## Figure 4: management performance metrics
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
  # filter(method!="DLM") %>%
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
        axis.line=element_line(linewidth=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,size=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(rows=vars(metric_label),cols=vars(sel_ordered),scales="free_y",switch="y")+theme(strip.text.x=element_text(size=9),strip.text.y=element_text(size=9))+theme(strip.placement="outside")
ggsave("Figure4.pdf",g,width=5.2,height=5.5,units="in")




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

##============================================================## plot return
p<-df_diff_plot2 %>%
  # filter(method!="DLM") %>%
  ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(linewidth=1) +
  geom_point(size=2.5) + 
  geom_point(shape=1,size=3,fill=NA,color="black") +
  scale_colour_manual(values=c("dodgerblue2","goldenrod1"))+
  scale_y_continuous(expand=c(0.1,0.1),breaks=seq(-20,30,5)) +
  theme_classic() +
  labs(x="",y="Median difference in returns (%)") + 
  labs(color="")+
  theme(strip.background=element_blank(),
        axis.line=element_line(linewidth=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,linewidth=1),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank()
  )
g<-p+facet_grid(cols=vars(sel_ordered)) +theme(strip.text.x=element_text(size=10),legend.position=c(0.08,0.94))
ggsave("Figure5.pdf",g,width=5,height=4,units="in")

##========================================================================##
##========================================================================##
##============================================================## end of code
##========================================================================##
##========================================================================##