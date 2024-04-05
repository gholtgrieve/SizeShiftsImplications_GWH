##==================================================================##
##                                                                  ##
##                Plot simulation model results for MSE             ##
##                                                                  ##
##==================================================================##
pkgs<-c("here","tidyverse","readxl","ggExtra")
if(length(setdiff(pkgs,rownames(installed.packages())))>0) {install.packages(setdiff(pkgs,rownames(installed.packages())),dependencies=TRUE)}
invisible(lapply(pkgs,library,character.only=T))

##==================================================================##
##=====================================================## load results
##==================================================================##
path<-paste0(here::here(),"/R/out/")

##=======================================================## simulation
# myfiles<-list.files(path,pattern="\\.Rdata$")
# myfile<-myfiles[grep("parms",myfiles,invert=TRUE)]
# load(paste0(here::here(),"/R/out/",myfile))
# timestamp<-substr(myfile,5,nchar(myfile)-6)

##=======================================================## saved file
timestamp<-"2024-03-22 13-19-52.529517"
load(paste0(path,"run_",timestamp,".Rdata"))

##=======================================================## parameters
parameters<-readRDS(paste0(path,"run_",timestamp,"_parms.Rdata"))
npara<-length(parameters)
for(i in 1:npara) { assign(names(parameters)[i],unlist(parameters[i])) }

##========================================================## scenarios
scenarios_all<-read_excel(paste0(path,"run_",timestamp,"_scen.xlsx"))[,-1]
scenarios_all<-data.frame(scenarios_all)
##-----------------------------------------## rename and order factors
scenarios_all$selectivity<-factor(
  scenarios_all$selectivity,
  levels=c("6 inch gillnet","unselective","8.5 inch gillnet"),
  labels=c("small-mesh","unselective","large-mesh"))
##-----------------------------------------## rename management method
scenarios_all<-scenarios_all %>%
  mutate(mgmt=case_when(
    mgmt=="smsy_goal" ~ "TRM",
    mgmt=="s_eq_goal" ~ "YPR",
    mgmt=="smsy_dlm_goal" ~ "DLM"))
##---------------------------------------------## rename S_MSY factors
scenarios_all$factorMSY<-factor(
  scenarios_all$factorMSY,
  levels=c("0.75","1","1.5"),
  labels=c("aggressive","1","precautionary"))

##==================================================================##
##====================================================## plot settings
##==================================================================##
colors<-c("deepskyblue3","orange")
##--------------------------------------------------------## functions
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) } 
summary_CI80<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.1),na.rm=T),ymax=quantile(x,prob=c(0.9),na.rm=T))) } 
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) } 
##--------------------------------------------------------## directory
dir.create(file.path(paste0(path,timestamp)),showWarnings=F)
setwd(file.path(paste0(path,timestamp)))

##==================================================================##
##=====================================================## factorMSY==1
##==================================================================##
scenarios<-scenarios_all[scenarios_all$factorMSY==1,]
nscen<-nscen_base_cases<-dim(scenarios)[1]		
##---------------------------------------------------## scenario names
scenario_names<-paste0(scenarios$trends," ",scenarios$mgmt," ", scenarios$selectivity)
##--------------------------------------------------## names of trends
trend_names<-unique(scenarios$trends)
n_trends<-length(trend_names)
##-------------------------## names of estimation methods used in mgmt
mgmt_names<-unique(scenarios$mgmt)
n_mgmt<-length(mgmt_names)
##-------------------------------------------## names on selectivities
selectivity_names<-unique(scenarios$selectivity)
n_select<-length(selectivity_names)

##==================================================================##
##=========================================================## Figure 2
##==================================================================##
review_years<-seq((nyi+20),ny,goalfreq)
nrev<-length(review_years)
myarray<-array(NA,dim=c(nscen,niter,nrev))
S_msy_scen<-myarray
for(j in 1:nscen) { 
  for(k in 1:niter) { 
    s_msy_est<-as.vector(unlist(S_msy.list[[j]][[k]]))
    if(length(s_msy_est)!=0) S_msy_scen[j,k,]<-s_msy_est
  } 
} 
##------------------------------------------------------## time period
## evaluated after 50 years, i.e. after the historical period
S_msy_hist<-S_msy_scen[,,which(review_years==nyh)]
##---------------------------------------------------## add scenarios
use_all<-dplyr::select(scenarios,trends,selectivity,mgmt)
S_msy_df<-data.frame(cbind(use_all,S_msy_hist)) %>%
  dplyr::filter(trends!="continuing trends")
##-----------------------------------------------------## long format
plot_smsy_df<-S_msy_df %>% 
  pivot_longer(!c(trends,selectivity,mgmt),
               names_to="iteration",
               values_to="value") %>% 
  dplyr::select(-iteration) %>%
  data.frame()
plot_smsy_df<-plot_smsy_df[complete.cases(plot_smsy_df),] ## <1%
##------------------------------------------------------------## plot
## 50% and 80% quantiles!
p<-plot_smsy_df %>% 
  ggplot(aes(x=fct_inorder(trends),y=value,color=fct_inorder(mgmt)))+
  scale_color_manual(values=c("darkgray",colors))+
  #stat_summary(fun.data=summary_CI90,position=jig,linewidth=0.1)+
  stat_summary(fun.data=summary_CI80,position=jig,linewidth=0.2)+
  stat_summary(fun.data=summary_CI50,position=jig,linewidth=0.6)+
  # coord_cartesian(ylim=c(4800,19500))+
  scale_y_continuous(breaks=seq(0,20000,4000),expand=c(0.02,0.02))+
  theme_classic()+
  labs(color="Method")+ 
  labs(x="",y=expression(""*S[MSY]*""))+ 
  facet_grid(cols=vars(selectivity),scales="free_y")+
  # facet_wrap(vars(selectivity),ncol=n_select,scales="free_y")+
  theme(strip.background=element_blank(),
        strip.text.x=element_text(size=10),
        axis.line=element_line(linewidth=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=10),
        panel.border=element_rect(fill=NA,linewidth=0.5),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8))
ggsave("Figure2.pdf",p,width=6,height=4,units="in")  

##====================================================## rename trends
scenarios_all<-scenarios_all %>%
  mutate(trends=case_when(
    trends=="age-sex-length trends" ~ "ASL trends stabilized",
    trends=="continuing trends" ~ "ASL trends continued",
    TRUE ~ as.character(trends)))

##==================================================================##
##=========================================================## Figure 3
##==================================================================##
## plot management performance metrics for each estimation method
##---------------------------## number of reconstructed observed years
ny_obs<-dim(data.frame(obs.list[[1]][[1]]))[1]
##------------------------------------------------## array for results
scen_names<-paste0("scen=",seq(nscen))
iter_names<-paste0("iter=",seq(niter))
dnames<-list(scen_names,iter_names)
myarray<-array(dim=c(nscen,niter),dimnames=dnames)
###---------------------------------------## performance last 50 years
year_index<-(nyh+1):ny_obs ## most recent years after historical period
##----------------------------------------------------------## metrics
av_harv<-av_ret<-av_esc<-p_over_Rmax50<-p_over_Seq50<-S_ratio<-myarray 
##------------------------------------## loop scenarios and iterations
for(j in 1:nscen) {
  for(k in 1:niter) {
    ## recruitment
    rec_jk<-data.frame(obs.list[[j]][[k]])$recRec[year_index]  
    if(is.null(rec_jk)) { next } else { 
      ## probability above/below threshold
      ricker_parms<-data.frame(sr_sim.list[[j]][[k]])
      R_max<-(ricker_parms$alpha/ricker_parms$beta)*exp(-1)
      p_over_Rmax50[j,k]<-length(which(rec_jk>R_max*0.5))/length(rec_jk)
    }
    ## stock recruit parameters and reference points
    ricker_parms<-data.frame(sr_sim.list[[j]][[k]])
    alpha<-ricker_parms$alpha
    beta<-ricker_parms$beta
    ## biological reference points
    S_eq<-log(alpha)/beta
    S_msy<-(1-lambert_W0(exp(1-log(alpha))))/beta
    S_max<-1/beta
    S_ratio[j,k]<-S_max/S_msy
    ## escapement
    esc_jk<-data.frame(obs.list[[j]][[k]])$obsEsc[year_index]  
    if(is.null(esc_jk)) { next } else {
      ## probability of escapement above 50% of S_zero
      p_over_Seq50[j,k]<-length(which(esc_jk>0.5*S_eq))/length(esc_jk)
    }
    ## harvest
    harv_jk<-data.frame(obs.list[[j]][[k]])$obsHarv[year_index]  
    if(is.null(harv_jk)) { next } else { 
      ## long-term average harvest
      av_harv[j,k]<-mean(harv_jk,na.rm=T)
    }
    ## return
    ret_jk<-data.frame(obs.list[[j]][[k]])$obsRet[year_index]  
    if(is.null(ret_jk)) { next } else { 
      ## long-term average return
      av_ret[j,k]<-mean(ret_jk,na.rm=T)
    }
  } ## end k-loop
} ## end j-loop
##--------------------------------------------------------## scenarios
df<-dplyr::select(scenarios,trends,selectivity,mgmt)
##-----------------------------------------------## median all metrics
df$av_harv<-apply(av_harv,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_harv<-df$av_harv/1e3
df$av_ret<-apply(av_ret,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_ret<-df$av_ret/1e3
df$p_over_Rmax50<-apply(p_over_Rmax50,1,function(x) median(x,na.rm=T))
df$p_over_Seq50<-apply(p_over_Seq50,1,function(x) median(x,na.rm=T))
##------------------------------------------------------## long format
dfp<-df %>% 
  dplyr::select(-S_ratio) %>%
  pivot_longer(!c(trends,selectivity,mgmt),names_to="metric",values_to="median") %>% 
  data.frame()
##-------------------------------------------## add labels for metrics
dfp<-dfp %>%
  mutate(metric_label=case_when(
    metric=="av_harv" ~ "Mean harvest\n(thousands)",
    metric=="av_ret" ~ "Mean return\n(thousands)",
    metric=="p_over_Rmax50" ~ "Probability\nabove 50% Rmax",
    metric=="p_over_Seq50" ~ "Probability\nabove 50% S0"))
##----------------------------------------------------## select trends
dfp<-dplyr::filter(dfp,trends!="age-length trends")
dfp$trends<-as.factor(as.character(dfp$trends))
##-------------------------------------------------------------## plot
p<-dfp %>% 
  ggplot(aes(x=fct_inorder(trends),y=median,
             col=fct_inorder(mgmt),group=mgmt)) +
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2.5,fill=NA,color="black") +
  geom_point(size=2.2) + 
  scale_colour_manual(values=c("darkgray",colors))+
  scale_y_continuous(expand=c(0.07,0.07),limits=c(0,NA)) +
  theme_classic() +
  labs(x="",y="") + 
  labs(color="Method")+
  theme(strip.background=element_blank(),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        axis.line=element_line(linewidth=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,linewidth=0.5),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8))
g<-p+facet_grid(rows=vars(metric_label),cols=vars(selectivity),scales="free_y",switch="y")+theme(strip.placement="outside")
ggsave("Figure3.pdf",g,width=5.5,height=7,units="in")

##=========================================================## S_ratios
## ratio of S_max over S_msy across stochastic simulations
S_ratios<-dplyr::select(scenarios,trends,selectivity,mgmt) %>%
  add_column(S_ratio=apply(S_ratio,1,function(x) median(x,na.rm=T)))%>%
  data.frame()
##---------------------------------------------------------## by trend
S_ratios_by_trend<-S_ratios %>%
  group_by(trends) %>%
  summarize(S_ratio=median(S_ratio))  
##---------------------------------------------------## by selectivity
S_ratios_by_selectivity<-S_ratios %>%
  group_by(selectivity) %>%
  summarize(S_ratio=median(S_ratio))  
##----------------------------------## by management/estimation method
S_ratios_by_mgmt<-S_ratios %>%
  group_by(mgmt) %>%
  summarize(S_ratio=median(S_ratio))  
##---------------------------------------------------## overall median
S_ratios_median<-median(S_ratios$S_ratio) ## ~1.5
## use for precautionary scenarios to emulate S_max as reference point

##==================================================================##
##=========================================================## Figure 4
##==================================================================##
## plot median % difference in return compared to time-invariant model
myarray<-array(dim=c(nscen,niter),dimnames=dnames)
av_ret_diff<-av_esc_diff<-myarray
##-----------------------------------------------------## get indices
ref_msy_ind<-which(grepl("TRM",scenario_names,fixed=T)) 
ypr_msy_ind<-which(grepl("YPR",scenario_names,fixed=T)) 
dlm_msy_ind<-which(grepl("DLM",scenario_names,fixed=T))  
##---------------------------------------------## compute differences
for(j in 1:nscen) {
  for(k in 1:niter) {
    if(j %in% ref_msy_ind) j_ref<-j
    if(j %in% ypr_msy_ind) j_ref<-j-4
    if(j %in% dlm_msy_ind) j_ref<-j-8
    av_esc_diff[j,k]<-(av_esc[j,k]-av_esc[j_ref,k])*100/av_esc[j_ref,k]
    av_ret_diff[j,k]<-(av_ret[j,k]-av_ret[j_ref,k])*100/av_ret[j_ref,k]
  } 
} 
##--------------------------------------------------------## scenarios
df_diff<-dplyr::select(scenarios,trends,selectivity,mgmt)
##---------------------------------------------------------## medians
df_diff$av_esc<-apply(av_esc_diff,1,function(x) median(x,na.rm=T))
df_diff$av_ret<-apply(av_ret_diff,1,function(x) median(x,na.rm=T))
##-----------------------------------------------------## long format
df_diff_plot<-df_diff %>% 
  pivot_longer(!c(trends,selectivity,mgmt),names_to="metric",values_to="median") %>% 
  data.frame()
##-----------------------------------------## add labels for metrics
df_diff_plot<-df_diff_plot %>%
  mutate(metric_label=case_when(
    metric=="av_esc" ~ "Mean escapement\n(thousands)",
    metric=="av_ret" ~ "Mean return\n(thousands)"))
##----------------------------------------------------## select trends
df_diff_plot<-dplyr::filter(df_diff_plot,trends!="age-length trends")
df_diff_plot$trends<-as.factor(as.character(df_diff_plot$trends))
##-------------------------------------------------## drop references
df_diff_plot<-df_diff_plot[df_diff_plot$mgmt!="TRM",]
##-----------------------------------------------------## plot return
p<-df_diff_plot %>%
  filter(metric=="av_ret") %>%
  ggplot(aes(x=fct_inorder(trends),y=median,
             color=fct_inorder(mgmt),group=mgmt)) +
  geom_hline(yintercept=0,linetype="solid",linewidth=0.1)+ 
  geom_line(linewidth=0.6) +
  geom_point(size=2.5) + 
  geom_point(shape=1,size=3,fill=NA,color="black") +
  scale_colour_manual(values=colors)+
  scale_y_continuous(expand=c(0.1,0.1),breaks=seq(-20,30,5)) +
  theme_classic() +
  labs(x="",y="Median difference in returns (%)") + 
  labs(color="Method")+
  facet_grid(cols=vars(selectivity)) +
  theme(strip.background=element_blank(),
        strip.text.x=element_text(size=10),
        axis.line=element_line(linewidth=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,linewidth=0.5),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank())
ggsave("Figure4.pdf",p,width=6,height=4,units="in")

##==================================================================##
##=========================================================## Figure 5
##==================================================================##
## plot % difference from time-invariant model over time 
##-------------------------------------------------## trends over time
year_index_list<-list()
future_yrs<-seq(10,50,10)
n_index<-length(future_yrs)
##-------------------------------------------------------## cumulative
for(i in 1:n_index) { year_index_list[[i]]<-(nyh+1):(nyh+i*10) }
##--------------------------## loop scenarios, iterations, and periods
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
} 
##--------------------------------------------------------## scenarios
df_diff_scen<-dplyr::select(scenarios,trends,selectivity,mgmt)
##----------------------------------------------------------## return
df_av_ret<-array(dim=c(nscen,n_index))
for(i in 1:n_index) {
  df_av_ret[,i]<-apply(av_ret_diff_list[[i]],1,function(x)median(x,na.rm=T))
  }
colnames(df_av_ret)<-future_yrs
##-------------------------------------------## combine with scenarios
df_av_ret<-data.frame(cbind(df_diff_scen,df_av_ret)) %>% 
  filter(mgmt!="TRM") %>%
  filter(trends=="ASL trends stabilized") %>%
  select(-trends)
##-----------------------------------------------------## pivot longer
df_av_ret_plot<-df_av_ret %>% 
  pivot_longer(!c(selectivity,mgmt),names_to="period",values_to="median")
df_av_ret_plot$period<-as.numeric(gsub("X","",df_av_ret_plot$period))
##------------------------------------------------------------## plot
p<-df_av_ret_plot %>% 
  ggplot(aes(x=period,y=median,color=fct_inorder(mgmt)))+
  geom_hline(yintercept=0,linetype="solid",linewidth=0.1)+ 
  geom_line(linewidth=0.5) +
  geom_point(shape=1,size=2.5,fill=NA,color="black") +
  geom_point(size=2.2) + 
  scale_colour_manual(values=colors)+
  scale_y_continuous(expand=c(0.1,0.1)) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  theme_classic() +
  labs(x="Year",y="Median difference in returns (%)") + 
  labs(color="Method")+
  facet_grid(cols=vars(selectivity),switch="y")+
  theme(strip.background=element_blank(),
        strip.text.x=element_text(size=10),
        strip.text.y=element_blank(),
        axis.line=element_line(linewidth=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        panel.border=element_rect(fill=NA,linewidth=0.5))
ggsave("Figure5.pdf",p,width=6.5,height=3,units="in")

##==================================================================##
##=========================================================## Figure 6
##==================================================================##
## harvest vs escapement % difference to time-invariant model all iterations
df_scen<-dplyr::select(scenarios,trends,selectivity,mgmt)
df_harv_long<-data.frame(cbind(df_scen,av_harv_diff)) %>% 
  pivot_longer(!c(selectivity,mgmt,trends),names_to="iter",values_to="harvest") %>% 
  dplyr::select(-iter)
df_esc_long<-data.frame(cbind(df_scen,av_esc_diff)) %>% 
  pivot_longer(!c(selectivity,mgmt,trends),names_to="iter",values_to="escapement") %>% 
  dplyr::select(-iter)
df_all<-data.frame(cbind(df_harv_long,escapement=df_esc_long$escapement))
##-------------------------------------------------------------## plot
p<-df_all %>% 
  filter(mgmt!="TRM") %>% ## alternatives relative to traditional
  filter(selectivity=="large-mesh") %>% ## large-mesh selectivity
  filter(trends=="ASL trends continued") %>% ## continued trends
  dplyr::select(-c(selectivity,trends)) %>%
  ggplot(aes(x=harvest,y=escapement,
             color=fct_inorder(mgmt),group=fct_inorder(mgmt)))+
  geom_hline(yintercept=0,linetype="solid",color="black",linewidth=0.1)+
  geom_vline(xintercept=0,linetype="solid",color="black",linewidth=0.1)+
  # geom_point(size=1,shape=16,alpha=0.1)+
  geom_point(size=0.75,shape=16)+
  scale_colour_manual(values=colors)+
  scale_x_continuous(limits=c(-80,190),breaks=seq(-50,150,50))+
  scale_y_continuous(limits=c(-80,190),breaks=seq(-50,150,50))+
  theme_classic()+
  labs(x="Difference in mean harvest (%)",
       y="Difference in mean escapement (%)")+
  labs(color="")+
  theme(strip.background=element_blank(),
        axis.line=element_line(linewidth=0.1),
        axis.text.x=element_text(angle=90,size=10,vjust=0.5,hjust=0.5),
        axis.text.y=element_text(size=10,angle=0,vjust=0.5,hjust=0.5),
        axis.title.x=element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        axis.title.y=element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        panel.border=element_rect(fill=NA,linewidth=0.5),
        legend.background=element_rect(fill='transparent'),
        # legend.box.background=element_rect(linewidth=0.1),
        # legend.position=c(0.9,0.1),
        legend.position=c(0.1,0.95),
        legend.key.size=unit(0.25,'cm'),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        plot.margin=unit(c(1,1,1,1),"lines"))
g<-ggMarginal(p,groupColour=TRUE,groupFill=TRUE,margins="both",size=5,type="density")
ggsave("Figure6.pdf",g,width=4,height=4,units="in")

##==================================================================##
##=========================================================## Figure 7
##==================================================================##
## scenarios assuming 0.75*S_MSY or 1.5*S_MSY management goals
##-----------------------------------------------------## factorMSY!=1
scenarios<-scenarios_all[scenarios_all$factorMSY!=1,]
scenarios$factorMSY<-as.factor(as.character(scenarios$factorMSY))
nscen<-dim(scenarios)[1]		
##---------------------------------------------------## scenario names
scenario_names<-paste0(scenarios$trends," ",scenarios$mgmt," ", scenarios$factorMSY)
##--------------------------------------------------## names of trends
trend_names<-unique(scenarios$trends)
n_trends<-length(trend_names)
##-------------------------## names of estimation methods used in mgmt
mgmt_names<-unique(scenarios$mgmt)
n_mgmt<-length(mgmt_names)
##------------------------------------## factors used in mgmt strategy
factorMSY_names<-unique(scenarios$factorMSY)
n_factors<-length(factorMSY_names)
##---------------------------## number of reconstructed observed years
ny_obs<-dim(data.frame(obs.list[[1]][[1]]))[1]
##------------------------------------------------## array for results
scen_names<-paste0("scen=",seq(nscen))
iter_names<-paste0("iter=",seq(niter))
dnames<-list(scen_names,iter_names)
myarray<-array(dim=c(nscen,niter),dimnames=dnames)
##------------------------------------## performance for last 50 years
year_index<-(nyh+1):ny_obs ## recent 50 years
##----------------------------------------------------------## metrics
av_harv<-av_ret<-av_esc<-p_over_Rmax50<-p_over_Seq50<-myarray 
##------------------------------------## loop scenarios and iterations
nbc<-nscen_base_cases ## factorMSY==1 cases
for(j in 1:nscen) {
  for(k in 1:niter) { 
    ## recruitment
    rec_jk<-data.frame(obs.list[[nbc+j]][[k]])$recRec[year_index]  
    if(is.null(rec_jk)) { next } else { 
      ## probability above/below threshold
      ricker_parms<-data.frame(sr_sim.list[[nbc+j]][[k]])
      R_max<-(ricker_parms$alpha/ricker_parms$beta)*exp(-1)
      p_over_Rmax50[j,k]<-length(which(rec_jk>R_max*0.5))/length(ret_jk)
    }
    ## escapement
    esc_jk<-data.frame(obs.list[[nbc+j]][[k]])$obsEsc[year_index]  
    if(is.null(esc_jk)) { next } else { 
      ## probability of escapement above 50% of S_zero
      ricker_parms<-data.frame(sr_sim.list[[nbc+j]][[k]])
      if(is.null(ricker_parms$alpha)) { next } else { 
        S_eq<-log(ricker_parms$alpha)/ricker_parms$beta
        p_over_Seq50[j,k]<-length(which(esc_jk>0.5*S_eq))/length(esc_jk)
      }
    }
    ## harvest
    harv_jk<-data.frame(obs.list[[nbc+j]][[k]])$obsHarv[year_index]  
    if(is.null(harv_jk)) { next } else { 
      ## long-term average harvest
      av_harv[j,k]<-mean(harv_jk,na.rm=T)
    }
    ## return
    ret_jk<-data.frame(obs.list[[nbc+j]][[k]])$obsRet[year_index]  
    if(is.null(ret_jk)) { next } else { 
      ## long-term average return
      av_ret[j,k]<-mean(ret_jk,na.rm=T)
    }
  } ## end k-loop
} ## end j-loop
##--------------------------------------------------------## scenarios
df<-dplyr::select(scenarios,trends,factorMSY,mgmt,selectivity)
##----------------------------------------------## median all metrics
df$av_harv<-apply(av_harv,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_harv<-df$av_harv/1e3
df$av_ret<-apply(av_ret,1,function(x) median(x,prob=use_quant,na.rm=T))
df$av_ret<-df$av_ret/1e3
df$p_over_Rmax50<-apply(p_over_Rmax50,1,function(x) median(x,na.rm=T))
df$p_over_Seq50<-apply(p_over_Seq50,1,function(x) median(x,na.rm=T))
##---------------------------------------------------------## filters
df<-df %>% 
  filter(selectivity=="unselective") %>% ## not needed for new scenarios
  dplyr::select(-selectivity) %>% 
  filter(trends!="age-length trends")
##-----------------------------------------------------## long format
dfp<-df %>% 
  pivot_longer(!c(trends,factorMSY,mgmt),names_to="metric",values_to="median")
##------------------------------------------## add labels for metrics
dfp<-dfp %>%
  mutate(metric_label=case_when(
    metric=="av_harv" ~ "Mean harvest\n(thousands)",
    metric=="av_ret" ~ "Mean return\n(thousands)",
    metric=="p_over_Rmax50" ~ "Probability\nabove 50% Rmax",
    metric=="p_over_Seq50" ~ "Probability\nabove 50% S0"))
##-------------------------------------------------------------## plot
p<-dfp %>% 
  ggplot(aes(x=fct_inorder(trends),y=median,
             color=fct_inorder(mgmt),group=mgmt)) +
  geom_line(lwd=0.5) +
  geom_point(shape=1,size=2.5,fill=NA,color="black") +
  geom_point(size=2.2) + 
  scale_colour_manual(values=c("darkgray",colors))+
  scale_y_continuous(expand=c(0.07,0.07),limits=c(0,NA)) +
  theme_classic() +
  labs(x="",y="") + 
  labs(color="Method")+
  theme(strip.background=element_blank(),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        axis.line=element_line(linewidth=0.1),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title=element_text(size=12),
        panel.border=element_rect(fill=NA,linewidth=0.5),
        legend.key.size=unit(0.5,'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8)
  )
g<-p+facet_grid(rows=vars(metric_label),cols=vars(factorMSY),scales="free_y", switch="y")+theme(strip.placement="outside")
ggsave("Figure7.pdf",g,width=4.5,height=7,units="in")

##==================================================================##
##==================================================================##
##==================================================================##