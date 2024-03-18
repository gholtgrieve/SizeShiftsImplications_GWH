##========================================================================##
##                                                                        ##
##       Function to generate age composition of brood year recruits      ##
##                                                                        ##
##========================================================================##
## ages: all simulated age classes
## recruits: number of recruits, by sex, from that brood year
## meanage: mean age of recruits, by sex, from that brood year
## sdage: standard deviation of age distribution (age diversity)
##========================================================================##
agecomp<-function(ages,recruits,meanage,sdage){
##------------------------------------## probabilities by age given mean age
probs_a<-dnorm(ages,meanage,sdage)
probs_a<-probs_a/sum(probs_a)
##---------------------------## draw dirichlet based on probabilities by age
d<-50 ## scales variability 
probs<-as.vector(rdirichlet(1,probs_a*d))
##-------------------------------------------------------------## age counts
agecounts<-round(probs*recruits)
return(agecounts) # return recruits by age class
}
##========================================================================##