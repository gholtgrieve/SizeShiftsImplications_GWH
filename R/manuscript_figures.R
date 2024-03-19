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
##===========================================================## load results
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

##===========================================## trends for old scenario file
scenarios_all<-scenarios_all %>%
  mutate(trends=case_when(
    trends=="age-size trends" ~ "age-length trends",
    trends=="age-size-sex trends" ~ "age-sex-length trends",
    TRUE ~ as.character(trends)
  ))

##========================================================================##
##==========================================================## plot settings
##========================================================================##
colors<-c("deepskyblue3","orange")
##--------------------------------------------------------------## functions
## 90% and 50% ranges (5th-95th and 25th-75th quantiles)
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) } 
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) } 

##========================================================================##
##===========================================================## factorMSY==1
##========================================================================##
scenarios<-scenarios_all[scenarios_all$factorMSY==1,]
nscen<-nscen_base_cases<-dim(scenarios)[1]		
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
##===============================================================## Figure 2
##========================================================================##
## trends in mean fecundity and egg mass over 50-year period 
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
##----------------------------------------------------------------## make df
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
  coord_cartesian(ylim=c(-73,40))+
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
##===============================================================## Figure 3
##========================================================================##
## plot differences in S_msy estimates compared to time-invariant model
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

##===========================================================## get legend
ymin<-min(S_msy_diff_med$median,na.rm=T)*1.1
ymax<-max(S_msy_diff_med$median,na.rm=T)*1.1
p<-S_msy_diff_med %>% 
  filter(method %in% c("YPR","DLM")) %>%
  ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method)))+
  scale_colour_manual(values=colors)+
  geom_point(size=2.5) + 
  theme_classic()+
  labs(color="Method")+
  theme(legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.position=c(0.5,0.65))

lgnd<-get_legend(p,position=NULL)
as_ggplot(lgnd)

##=========================================================## plot quantiles
jig<-position_dodge(width=0.5)
cols<-rep(colors,n_diff)
ymin<-quantile(plot_smsy$value,probs=0.01,na.rm=T)
ymax<-quantile(plot_smsy$value,probs=0.99,na.rm=T)
p<-plot_smsy %>% ggplot(aes(x=fct_inorder(trends),y=value,fill=fct_inorder(method)))+
  geom_hline(yintercept=0,linetype="solid",linewidth=0.1)+ 
  stat_summary(fun.data=summary_CI90,position=jig,linewidth=0.2,color=cols)+
  stat_summary(fun.data=summary_CI50,position=jig,linewidth=0.6,color=cols)+
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
g<-cowplot::plot_grid(g,lgnd,ncol=2,rel_widths=c(0.85,0.15)) ## add legend
ggsave("Figure3.pdf",g,width=6,height=4,units="in")  

##========================================================================##
##===============================================================## Figure 4
##========================================================================##
## plot management performance metrics for each estimation method
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
##==============================================## performance last 50 years
year_index<-(nyh+1):ny_obs ## most recent years after historical period
##--------------------------------------------------------## for each method
sum_log_harv<-cv_harv<-av_harv<-av_log_harv<-p_no_harv<-myarray ## harvest 
sum_log_esc<-cv_esc<-av_esc<-av_log_esc<-myarray ## escapement 
sum_log_ret<-cv_ret<-av_ret<-av_log_ret<-myarray ## total return 
p_below_thresh<-p_above_thresh<-myarray ## conservation threshold
thresh<-0.5 ## proportion of maximum recruitment at equilibrium
##------------------------------------------## loop scenarios and iterations
for(j in 1:nscen) {
  for(k in 1:niter) { 
    ##---------------------------------------------------------## escapement
    esc_jk<-data.frame(obs.list[[j]][[k]])$obsEsc[year_index] ## observed 
    if(is.null(esc_jk)) { next } else { 
      sum_log_esc[j,k]<-sum(log(esc_jk+0.0001))
      cv_esc[j,k]<-1/(sd(esc_jk,na.rm=T)/mean(esc_jk,na.rm=T)) ## 1/CV
      av_esc[j,k]<-mean(esc_jk,na.rm=T)
      av_log_esc[j,k]<-mean(log(esc_jk),na.rm=T)
    }
    ##----------------------------------------------------------## harvest
    harv_jk<-data.frame(obs.list[[j]][[k]])$obsHarv[year_index] ## observed 
    if(is.null(harv_jk)) { next } else { 
      sum_log_harv[j,k]<-sum(log(harv_jk+0.0001))
      cv_harv[j,k]<-1/(sd(harv_jk,na.rm=T)/mean(harv_jk,na.rm=T)) ## 1/CV
      av_harv[j,k]<-mean(harv_jk,na.rm=T)
      av_log_harv[j,k]<-mean(log(harv_jk),na.rm=T)
      p_no_harv[j,k]<-length(which(harv_jk==0))/length(harv_jk)
    }
    ##----------------------------------------------------------## return
    ret_jk<-data.frame(obs.list[[j]][[k]])$obsRet[year_index] ## observed 
    if(is.null(ret_jk)) { next } else { 
      sum_log_ret[j,k]<-sum(log(ret_jk+0.0001))
      cv_ret[j,k]<-1/(sd(ret_jk,na.rm=T)/mean(ret_jk,na.rm=T)) ## 1/CV
      av_ret[j,k]<-mean(ret_jk,na.rm=T)
      av_log_ret[j,k]<-mean(log(ret_jk),na.rm=T)
      ## maximum recruitment given Ricker parameters
      ricker_parms<-data.frame(sr_sim.list[[j]][[k]]) 
      max_rec<-(ricker_parms$alpha/ricker_parms$beta)*exp(-1) 
      ## probability above/below threshold
      p_below_jk<-length(which(ret_jk<max_rec*thresh))/length(ret_jk)
      p_below_thresh[j,k]<-p_below_jk 
      p_above_jk<-length(which(ret_jk>max_rec*thresh))/length(ret_jk)
      p_above_thresh[j,k]<-1-p_below_jk
    }
  } ## end k-loop
} ## end j-loop
p_below_thresh[is.na(p_below_thresh)]<-0

##===========================================================## plot medians
df<-dplyr::select(scenarios,trends,selectivity,mgmt)
##---------------------------------------------------------------## 'method'
df$method<-"TRM"
df$method[grepl("eq",df$mgmt)]<-"YPR"
df$method[grepl("dlm",df$mgmt)]<-"DLM"
##----------------------------------------------------## median all metrics
df$av_harv<-apply(av_harv,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_harv<-df$av_harv/1e3
df$cv_harv<-apply(cv_harv,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_ret<-apply(av_ret,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_ret<-df$av_ret/1e3
df$p_above_thresh<-apply(p_above_thresh,1,function(x) median(x,na.rm=T))
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
  filter(metric!="cv_harv") %>%
  ggplot(aes(x=fct_inorder(trends),y=median,col=fct_inorder(method),group=method)) +
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2.5,fill=NA,color="black") +
  geom_point(size=2.2) + 
  scale_colour_manual(values=c("darkgray",colors))+
  scale_y_continuous(expand=c(0.08,0.08)) +
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
g<-p+
  facet_grid(rows=vars(metric_label),cols=vars(sel_ordered),scales="free_y",switch="y")+
  theme(strip.text.x=element_text(size=9),strip.text.y=element_text(size=9))+
  theme(strip.placement="outside")
ggsave("Figure4.pdf",g,width=5.2,height=5.5,units="in")

##========================================================================##
##===============================================================## Figure 5
##========================================================================##
## plot median percent difference in return compared to time-invariant model

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

##=======================================================## plot differences
df_diff<-dplyr::select(scenarios,trends,selectivity,mgmt)
##---------------------------------------------------------------## 'method'
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
  ggplot(aes(x=fct_inorder(trends),y=median,color=fct_inorder(method),group=method)) +
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(linewidth=0.6) +
  geom_point(size=2.5) + 
  geom_point(shape=1,size=3,fill=NA,color="black") +
  scale_colour_manual(values=colors)+
  scale_y_continuous(expand=c(0.1,0.1),breaks=seq(-20,30,5)) +
  theme_classic() +
  labs(x="",y="Median difference in returns (%)") + 
  labs(color="Method")+
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
g<-p+facet_grid(cols=vars(sel_ordered)) +
  theme(strip.text.x=element_text(size=10))
ggsave("Figure5.pdf",g,width=6,height=4,units="in")

##========================================================================##
##===============================================================## Figure 6
##========================================================================##
## plot median percent difference from time-invariant model over time 

##===========================================## performance trends over time
year_index_list<-list()
future_yrs<-seq(10,50,10)
n_index<-length(future_yrs)
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
} ## end i-loop over year indices

##==============================================================## scenarios
df_diff_scen<-dplyr::select(scenarios,trends,selectivity,mgmt)
##----------------------------------------------------## 'method'
df_diff_scen$method<-"TRM"
df_diff_scen$method[grepl("eq",df_diff_scen$mgmt)]<-"YPR"
df_diff_scen$method[grepl("dlm",df_diff_scen$mgmt)]<-"DLM"

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
df_av_ret_plot$sel_ordered<-factor(df_av_ret_plot$selectivity,levels=c("6 inch gillnet", "unselective", "8.5 inch gillnet"),labels=c("small-mesh","unselective","large-mesh"))
##-----------------------------------------------------------## metric label
df_av_ret_plot$metric_label<-"Mean return"

##=================================================================## return
df_plot<-df_av_ret_plot
p<-df_plot %>% 
  ggplot(aes(x=period,y=median,color=fct_inorder(method)))+
  geom_hline(yintercept=0,linetype="solid",size=0.1)+ 
  geom_line(linewidth=0.5) +
  geom_point(shape=1,size=2.5,fill=NA,color="black") +
  geom_point(size=2.2) + 
  scale_colour_manual(values=colors)+
  scale_y_continuous(expand=c(0.1,0.1)) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="Year",y="Median difference in returns (%)") + 
  labs(color="Method")+
  theme(strip.background=element_blank(),
        axis.line=element_line(size=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        panel.border=element_rect(fill=NA,size=1)
  )
g<-p+
  facet_grid(cols=vars(sel_ordered),rows=vars(metric_label),switch="y")+
  theme(strip.text.x=element_text(size=10),strip.text.y=element_blank())
ggsave("Figure6.pdf",g,width=6.5,height=3,units="in")

##========================================================================##
##===============================================================## Figure 7
##========================================================================##
## harvest vs escapement % difference to time-invariant model all iterations 
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

##=================================================## plot including density
p<-df_all %>% 
  dplyr::filter(sel_ordered=="large-mesh") %>% ## large-mesh selectivity
  ggplot(aes(x=harvest,y=escapement,col=fct_inorder(method),group=fct_inorder(method)))+
  geom_hline(yintercept=0,linetype="solid",color="black",size=0.1)+
  geom_vline(xintercept=0,linetype="solid",color="black",size=0.1)+
  geom_point(size=1,shape=16,alpha=0.5)+
  scale_colour_manual(values=colors)+
  scale_x_continuous(limits=c(-80,190),breaks=seq(-50,150,50))+
  scale_y_continuous(limits=c(-80,190),breaks=seq(-50,150,50))+
  theme_classic()+
  labs(x="Difference in mean harvest (%)",
       y="Difference in mean escapement (%)")+
  labs(color="")+
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
g<-ggMarginal(p,groupColour=TRUE,groupFill=TRUE,margins="both",size=5,type="density")
ggsave("Figure7.pdf",g,width=4,height=4,units="in")

##========================================================================##
##===============================================================## Figure 8
##========================================================================##

##===========================================================## factorMSY!=1
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

##=========================================## management performance metrics
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
##------------------------------------------## performance for last 50 years
year_index<-(nyh+1):ny_obs ## recent 50 years
##--------------------------------------------------------## for each method
sum_log_harv<-cv_harv<-av_harv<-av_log_harv<-p_no_harv<-myarray ## harvest 
sum_log_esc<-cv_esc<-av_esc<-av_log_esc<-myarray ## escapement 
sum_log_ret<-cv_ret<-av_ret<-av_log_ret<-myarray ## total return 
p_below_thresh<-p_above_thresh<-myarray ## conservation threshold
thresh<-0.5 ## proportion of maximum recruitment at equilibrium
##------------------------------------------## loop scenarios and iterations
nbc<-nscen_base_cases ## factorMSY==1 cases
for(j in 1:nscen) {
  for(k in 1:niter) { 
    ##---------------------------------------------------------## escapement
    esc_jk<-data.frame(obs.list[[nbc+j]][[k]])$obsEsc[year_index]  
    if(is.null(esc_jk)) { next } else { 
      sum_log_esc[j,k]<-sum(log(esc_jk+0.0001))
      cv_esc[j,k]<-1/(sd(esc_jk,na.rm=T)/mean(esc_jk,na.rm=T)) ## 1/CV
      av_esc[j,k]<-mean(esc_jk,na.rm=T)
      av_log_esc[j,k]<-mean(log(esc_jk),na.rm=T)
    }
    ##----------------------------------------------------------## harvest
    harv_jk<-data.frame(obs.list[[nbc+j]][[k]])$obsHarv[year_index]  
    if(is.null(harv_jk)) { next } else { 
      sum_log_harv[j,k]<-sum(log(harv_jk+0.0001))
      cv_harv[j,k]<-1/(sd(harv_jk,na.rm=T)/mean(harv_jk,na.rm=T)) ## 1/CV
      av_harv[j,k]<-mean(harv_jk,na.rm=T)
      av_log_harv[j,k]<-mean(log(harv_jk),na.rm=T)
      p_no_harv[j,k]<-length(which(harv_jk==0))/length(harv_jk)
    }
    ##----------------------------------------------------------## return
    ret_jk<-data.frame(obs.list[[nbc+j]][[k]])$obsRet[year_index]  
    if(is.null(ret_jk)) { next } else { 
      sum_log_ret[j,k]<-sum(log(ret_jk+0.0001))
      cv_ret[j,k]<-1/(sd(ret_jk,na.rm=T)/mean(ret_jk,na.rm=T)) ## 1/CV
      av_ret[j,k]<-mean(ret_jk,na.rm=T)
      av_log_ret[j,k]<-mean(log(ret_jk),na.rm=T)
      ## maximum recruitment given Ricker parameters
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

##===========================================================## plot medians
df<-dplyr::select(scenarios,trends,factorMSY,mgmt,selectivity)
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
df$p_above_thresh<-apply(p_above_thresh,1,function(x) median(x,na.rm=T))
##----------------------------------------------------------## metric labels
labs<-c("Mean harvest\n(thousands)","Harvest stability\n(1/CV)","Mean return\n(thousands)","Probability\nabove threshold")
n_metrics<-length(labs)
df_factorMSY<-df

##===========================================================## selectivity
## make sure to use only 'unselective' in case of other scenarios were run
df<-df_factorMSY %>% 
  dplyr::filter(selectivity=="unselective") %>%
  dplyr::select(-selectivity)

##==========================================================## long format
dfp<-df %>% pivot_longer(!c(trends,factorMSY,mgmt,method), names_to="metric", values_to="median") %>% data.frame()
##-----------------------------------------------## add labels for metrics
metrics<-data.frame(name=colnames(df)[-c(1:4)])
metrics$label<-labs
for(i in 1:dim(dfp)[1]) {
  dfp$metric_label[i]<-metrics$label[dfp$metric[i]==metrics$name]
}
##----------------------------------------------------------## select trends
dfp<-dplyr::filter(dfp,trends!="age-length trends")
dfp$trends<-as.factor(as.character(dfp$trends))
dfp_factorMSY<-dfp

##=================================================## medians by target type
p<-dfp_factorMSY %>% 
  dplyr::filter(metric!="cv_harv") %>%
  ggplot(aes(x=fct_inorder(trends),y=median,color=fct_inorder(method),group=method)) +
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2.5,fill=NA,color="black") +
  geom_point(size=2.2) + 
  scale_colour_manual(values=c("darkgray",colors))+
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
g<-p+facet_grid(rows=vars(metric_label),cols=vars(factorMSY),scales="free_y", switch="y")+
  theme(strip.text.x=element_text(size=9),strip.text.y=element_text(size=9))+
  theme(strip.placement="outside")
ggsave("Figure8.pdf",g,width=4.2,height=5.5,units="in")

##========================================================================##
##========================================================================##
##========================================================================##