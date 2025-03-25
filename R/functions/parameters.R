## List of input parameters for simulation model
parms.list<-list( 

##################### Model Operating Parameters ####################
seednum=seednum, ## seed
j=j, ## scenario
k=k, ## iteration
print=F, ## print progress for diagnostics 
plot=F, ## plot data during model run (only first few iterations)
sim_recruits="eggmass", ##  from 'spawners' 'fecundity' 'eggmass'
est_method="GLS", ## SR fit method ('GLS' - previously 'JAGS','STAN')
ricker_type="const_beta", ## 'const_beta' or 'const_rmax'
var="both", ## vary 'alpha' or 'beta' or 'both' in DLM?
optimize="logFmax", ## optimize F as 'Fmax' or logFmax'
nyi=10, ## initial years not analyzed (to produce return)
nyh=50, ## number of 'historical' years with trends
ny=110, ## number of total years with trends (ny-10 reconstructed)
futureT=scen$futureT[j], ## continue trends beyond historical period?
#####################################################################


################### Harvest Management Parameters ##################
obserr=0.2, ## observation error of escapement (default=0.2)
hobserr=0.1, ## observation error of harvest (default=0.1)
harverr=0.15, ## harvest implementation error (default=0.3) 
reglength=20, ## inverse probability of regime shift
regstr=1, ## low/high productivity factor (default=1)
harvmgmt=paste0(scen$mgmt)[j], ## 'smsy', 'umsy' or 'fix_harv_rate'
factorMSY=scen$factorMSY[j], ## aggressive or precautionary
harvrate=0.5, ## harvest rate (only used for 'fix_harv_rate')
goalfreq=10, ## review freq in yrs ('F'=no review)
#####################################################################


################### Recruitment Parameters ##########################
###  Updated for Kuskokwim based on Staton et al. 2021, Supplemental 
###  Material, N0 model, section 10, Table 4a.      #################
procerr=0.38, ## recruitment process error (median sigma_R in Table 4a)
rho=0.593, ## recruitment autocorrelation (median phi in Table 4a)
alpha_mean=6.463, ## mean productivity at low abundance (median alpha in Table 4a)
beta_mean=1.039e-5, ## mean density effect (median beta_e5 in Table 4a but replace positive exponent with negative, i.e., e5 becomes e-5
sr_corr=0.0, ## strength of alpha-beta correlation (default=0)
sr_parms_sd=0.2, ## stochastic variation SR parameters
#####################################################################


#################### Selectivity Parameters ########################
maxsel=scen$maxsel[j], ## size of maximum selectivity (default=0.2)
sdsel=scen$sdsel[j], ## selectivity (see Bromaghin 2005)
#####################################################################


#################### Age-Sex-Length Parameters ######################
ages=seq(1,9,1), ## all possible ages including FW and marine life stages.
meanageini=5.5, ## initial mean age, Confirmed as 5.5 based on Staton, et al. 2021, Supplemental Material, EM-ASL scenario, section 6, Fig. 7a.
# Calculated weighted average of ages in 1970s as weighted average across year classes and sexes.
agetrend=-0.4*scen$ageT[j], ## mean age trend.  This is based on an observed decline in average age of population presented in (default=-0.4). 
# Staton, et al. 2021, Supplemental Material, EM-ASL scenario, section 6, Fig. 7a.
# Calculated as difference in weighted average of ages in 1970 and 2015.
sdage=0.6, ## variation in cohort age composition (default=0.6)
vonB_Linf=1200, ## von B average asymptotic length (default=1200). #might need to change 
vonB_k=0.325, ## von B growth rate coefficient (default=0.325). #might need to change 
ocean0s=100, ## size in mm of ocean-0 fish.  Originally 150 in model.  Decreased to 100 assuming that Kusko Chinook are smaller than average at all age classes.  This is a guess.
sdSaA=0.01, ## sd of size-at-age anomalies (default=0.01)
sizetrends=c(0,0,0,0,10,-40,rep(-100,3))*scen$sizeT[j], ## change in size by age class over time period in mm. 
# Based on visually estimated change 1970-2015 shown in Olhberger et al. 2020, Fig 5d. No data for ages 1-3, so set to 0. Age 4 shows no change, age 5 
# shows increase of ~10 mm, age 6 shows decline of ~40 mm, and age 7 shows decline of ~100 mm, which is assumed the same for ages 8 and 9 (no data).
propF=0.41, ## proportion female in return, set to 0.41 based on Fig 5e in Ohlberger et al. 2020.  Also estimated at 0.41 in Staton et al. 2021,
# Supplemental Material, EM-ASL scenario, section 6, Fig. 7b.
propFtrend=-0.275*scen$sexT[j], ## trend in proportion female on logit scale. This is based on an observed decline in propotion female 
# in Staton, et al. 2021, Supplemental Material, EM-ASL scenario, section 6, Fig. 7b. Visually estimated to change from 0.41 to 0.3 
# between 1970 and 2015 (~27.5% decline).  This is a more conservative change than suggested in Fig 5e in Ohlberger et al. 2020.
agediff=1, ## age difference between sexes (default= 1)
allometry=c(size_ref=800,
            fec_ref=6600,
            b_fec=2.4,  #from Yukon, keep for now
            egg_ref=916,
            b_eggs=4.8), #from Yukon, keep for now
alt_sr_param=c(alpha_N0=5.07,
               alpha_EASL=0.00177,
               alpha_EMASL=0.01036,
               beta_N0=8.6e-6,
               beta_EASL=2.973e-09,
               beta_EMASL=1.6967e-08) ## see Staton et al. 2021 
#####################################################################
### Ohlberger, J., D.E. Schindler, R.J. Brown, J.M.S. Harding, M.D. Adkison, and A.R. Munro.
### 2020. Analysis of Changes in Quality of Chinook Salmon Escapement in the AYK Region.
### Arctic-Yukon-Kuskokwim Sustainable Salmon Initiative. Anchorage, AK. 47 p. + appendix.
###
### Staton, J., J. Ohlberger, M. Adkison, and A. Munro. 2021. A Bayesian state-space model for
### estimating productivity and recruitment dynamics of Kuskokwim River Chinook salmon.
### Fisheries Research 239:105961.
)
