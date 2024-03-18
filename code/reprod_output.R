##========================================================================##
##                                                                        ##
##     Function to calulate fecundity or egg mass based on female size    ##
##                                                                        ##
##========================================================================##
## uses exponent and reference values to compute log-log intercept/constant
## allows for easy manipulation of exponent without adjusting the constant
##========================================================================##
## size: female length in mm
## allometry: parameters of the allometirc relationships
## a_fec/a_eggs: intercepts on log-log scale
## b_fec/b_eggs: slopes on log-log scale
##========================================================================##
reprod_output<-function(size,allometry){
size_ref<-allometry[names(allometry)=="size_ref"] ## reference size
##----------------------------------------## allometry for fecundity
fec_ref<-allometry[names(allometry)=="fec_ref"] ## fecundity at reference size
b_fec<-allometry[names(allometry)=="b_fec"] ## log-log slope for fecundity
a_fec<-exp(log(fec_ref)-b_fec*log(size_ref)) ## log-log intercept
##----------------------------------------## calculate fecundity
fecundity<-log(a_fec)+b_fec*log(size)
##----------------------------------------## allometry for egg mass
egg_ref<-allometry[names(allometry)=="egg_ref"] ## egg mass at reference size
b_eggs<-allometry[names(allometry)=="b_eggs"] ## log-log slope for egg mass
a_eggs<-exp(log(egg_ref)-b_eggs*log(size_ref)) ## log-log intercept
##----------------------------------------## calculate egg mass
eggmass<-log(a_eggs)+b_eggs*log(size) 
##----------------------------------------## return both as list
return(list(fecundity,eggmass)) 
}
##========================================================================##