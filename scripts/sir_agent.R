##################################################################################
##################################################################################
# An R script to perform a stochastic epidemic simulation using an
# Agent Based Model and homogeneous mixing in the population.
#
# Author:  Sherry Towers
#          smtowers@asu.edu
# Created: Feb 12, 2016
#
# Copyright Sherry Towers, 2016
#
# This script is not guaranteed to be free of bugs and/or errors
#
# This script can be freely used and shared as long as the author and 
# copyright information in this header remain intact.
##################################################################################
set.seed(578194)

require("deSolve")
source("scripts/sir_agent_func.R")
nrealisations = 25 # Number of Runs

##################################################################################
##################################################################################
# this is a function which, given a value of S,I and R at time t
# calculates the time derivatives of S I and R
# vparameters contains the parameters of the model, like the
# recovery period, gamma, and the transmission rate, beta
# this function gets passed to the deSolve package
##################################################################################
SIRfunc=function(t, x, vparameters){
   S = x[1]  
   I = x[2]  
   R = x[3]  
   if (I<0) I=0 

   with(as.list(vparameters),{
      npop = S+I+R   
      dS = -beta*S*I/npop            
      dI = +beta*S*I/npop - gamma*I  
      dR = +gamma*I                  
      out = c(dS,dI,dR)
      list(out)
   })
}



##################################################################################
##################################################################################
# Set up initial conditions
##################################################################################
N = 1000       # population size
I_0 = 5        # number intially infected people in the population
S_0 = N - I_0
E_0 = 0 * N    # assume no one is exposed at first
R_0 = 0 * N    # assume no one has recovered at first
De  = 1/5      # incubation period of COVID-19 days^{-1}
delta_t = .5          # nominal time step
tbeg  = 0             # begin day
tend  = 80            # end day
gamma = 1/10          # recovery period of influenza in days^{-1}
R0    = 1.34          # R0 of a hypothetical strain of pandemic influenza
beta = R0 * gamma     # "reverse engineer" beta from R0 and gamma

##################################################################################
# first simulate the model with deterministic ODE's, so that we have something
# to compare our stochastic simulation to.
##################################################################################
vt = seq(tbeg,tend,delta_t)
vparameters = c(gamma=gamma,beta=beta)
inits = c(S=S_0,I=I_0,R=R_0)

sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc, vparameters))

##################################################################################
# now plot the results of the deterministic model
##################################################################################
par(mfrow=c(3,1))  # divides the page into two plotting areas 
plot(sirmodel$time,sirmodel$I/N,ylim=c(0,1.5*max(sirmodel$I/N)),type="l",col=1,lwd=5,xlab="Time, in days",ylab="Fraction infected (prevalence)",main=paste("COVID-19 pandemic in population of ",N,sep=""))
cat("The final size of epidemic from the deterministic model is ",max(sirmodel$R/N),"\n")

##################################################################################
# now do several simulations using the agent based model, and overlay the
# results on those from the deterministic model
##################################################################################
vfinal = numeric(0) # we will fill this vector with the epidemic final size estimates from the simulations
all_agents = factor("list", )
for (iter in 1:nrealisations){
   myagent = SIR_agent(N,
                       I_0,
                       S_0,
                       E_0,
                       R_0,
                       De,
                       gamma,
                       R0,
                       tbeg,
                       tend,
                       delta_t)
   lines(myagent$time,myagent$I/N,lwd=2,col=(iter+1),lty=3)
   cat(iter,nrealisations,"The final size of the epidemic from the agent based stochastic model is ",myagent$final_size,"\n")
   vfinal = append(vfinal,myagent$final_size) 
}
lines(sirmodel$time,sirmodel$I/N,lwd=5)
cat("The final size of epidemic from the deterministic model is ",max(sirmodel$R/N),"\n")
legend("topright",legend=c("Deterministic","Agent Based simulation"),col=c(1,2),lwd=3,lty=c(1,3),bty="n")

hist(vfinal,xlab="Distribution of epidemic final size",main="") # histogram the final size
lines(c(max(myagent$I/N),max(sirmodel$I/N)),c(-1000,1000),col=4,lwd=2,lty=2)
legend("topleft",legend=c("Agent Based simulation final size","Deterministic model final size"),col=c(1,4),lwd=2,lty=c(1,2),bty="n")

