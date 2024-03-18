##========================================================================##
##                                                                        ##
##                     Harvest Selectivity Function                       ##
##                                                                        ##
##========================================================================##
## size: fish size in mm
## meshsize: size of stretch mesh in inches
## s: standard deviation of selectivity (Bromaghin 2005: s=0.204)
##========================================================================##
selectivity<-function(size,meshsize,s) { 
t<-1.920
h<-0.622
l<--0.547
# s<-0.204 ## input such that it can be set to large value for 'unselective'
perimeter<-25.4*meshsize*2
x<-size/perimeter
select<-(1+l^2/(4*h^2))^h*(1+(x-s*l/2*h-t)^2/(s^2))^-h*exp(-l*(atan((x-s*l/2*h-t)/s)+atan(l/(2*h))))
return(select)	
}