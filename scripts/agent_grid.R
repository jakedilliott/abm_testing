##################################################################################
# An R script to run an Agent Based Model with connectedness among the population
# on a square grid.  Infected people can only infect their nearest neighbors on
# the grid.
#
# Author: Sherry Towers
#         admin@sherrytowers.com
#
# Created: Dec 11th, 2012
# Copyright Sherry Towers, 2012
#
# This script is not guaranteed to be free of bugs and/or errors.
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
##################################################################################

##################################################################################
# set up the parameters of our simulation
##################################################################################
nxgrid = 20       # the number of rows and columns in the grid
npop = (nxgrid)^2  
ninitial = 3      # the initial number of people infected in the population

R0 = 3         # Average number of people each infected person infects during the course of their infection
gamma = 1/3    # Recovery rate, in days^{-1}

time_begin = 0    # begin time (in days)
time_end   = 200  # end time
delta_t = 0.10    # The time step, in days

##################################################################################
# set up the positions of the points on the grid, for plotting purposes
# here I add and extra set of points all along the edge of the population, such that
# all population members that I will consider will have 8 neighbors.
#
##################################################################################
vx = seq(1,(nxgrid+2),1) 
vy = seq(1,(nxgrid+2),1)

##################################################################################
# if state is 0 then susceptible, if 1 then infected, if 2 then recovered, if 3 then edge
##################################################################################
wstate = matrix(0,nrow=(nxgrid+2),ncol=(nxgrid+2))
wx = matrix(0,nrow=(nxgrid+2),ncol=(nxgrid+2))
wy = matrix(0,nrow=(nxgrid+2),ncol=(nxgrid+2))
for (i in 1:length(vx)){
  for (j in 1:length(vy)){
     wx[i,j] = vx[i]
     wy[i,j] = vy[j]
     if (i==1|i==(nxgrid+2)) wstate[i,j] = 3
     if (j==1|j==(nxgrid+2)) wstate[i,j] = 3
  } # end loop over vy
} # end loop over vx



##################################################################################
# pick ninitial of the people to be initially infected
##################################################################################
iindex = sample(npop,ninitial) # this randomly chooses the indices of ninitial people out of the sample npop
jindex = which(wstate==3)      # wstate=3 
iindex = iindex[!(iindex%in%jindex)]

while (length(iindex)<ninitial){
   iindex = sample(npop,ninitial)
   iindex = iindex[!iindex%in%jindex]
}
wstate[iindex] = 1


##################################################################################
# let's plot our grid... 
# susceptibles will be lightgray
# infected will be red
# recovered (and immune) will be in blue
#
# scol is a vector of colors, and we will index it using the state vector
##################################################################################
par(pch=20)
scol=c("gray","red","blue",0)
plot(wx,wy,col=scol[(wstate+1)],cex=0.5,axes=F,xlab="",ylab="",main="S=grey   I=red   R=blue")

##################################################################################
# now step through the simulation, in time steps delta_t
# until either we reach the time=time_end, or there are no more infected people
##################################################################################
time = time_begin                 # begin the simulation at time=time_begin
ninf = length(wstate[wstate==1])  # count the number in the population with state=1 (ie; the number of infected)

while(ninf>0&time<time_end){ # continue the simulation until end time, or no more infectives (which ever comes first)
  ################################################################################
  # now loop over the infectives in the population.
  # for each infective, get the indices of their neighbors in the grid
  ################################################################################
  indices_of_infectives = which(wstate==1) # get the indices in the state vector of the infectives
  for (i in 1:length(indices_of_infectives)){  # loop over the infectives

    points(wx,wy,col=scol[(wstate+1)],cex=0.5) # plot the population, with the colors that correspond to their state

    ##############################################################################
    # now get the position of each infective, and identify the neighbors to that
    # infective
    ##############################################################################
    ax = wx[indices_of_infectives[i]] # the x position of the infective
    ay = wy[indices_of_infectives[i]] # the y position of the infective

    indices_of_neighbors = 
             which(
                    (wx==(ax-1)&(wy==(ay-1)|wy==(ay+1))) # kitty corner neighbors on the left
                   |(wx==(ax+1)&(wy==(ay-1)|wy==(ay+1))) # kitty corner neighbors on the right
                   |(wx==ax&(wy==(ay-1)|wy==(ay+1)))     # neighbors above and below
                   |(wy==ay&(wx==(ax-1)|wx==(ax+1)))     # neighbors to left and right
                  )

    state_vector_of_neighbors = wstate[indices_of_neighbors] # the state vector of the neighbors 

    ##############################################################################
    # count the number of susceptibles
    ##############################################################################
    nsus = length(state_vector_of_neighbors[state_vector_of_neighbors==0]) # the number of susceptibles in the neighbors
    ntouch = min(rpois(1,R0*delta_t*gamma),8) # sample a random number of neighbors for the infectives to touch in this
                                              # time step with the Poisson distribution. This number cannot exceed 8.
    if (ntouch>0&nsus>0){
      vtouch = sample(8,ntouch)  # out of the 8 neighbors pick ntouch of them randomly
      indices_of_touched = indices_of_neighbors[vtouch]
      wstate[indices_of_touched][wstate[indices_of_touched]==0] = 1 # if susceptible neighbor was touched, change them to infected
    }

    ##############################################################################
    # now check if the infected person recovers in this time step
    ##############################################################################
    prandom = runif(1)
    if (prandom<(delta_t*gamma)){ # check to see if the infected person recovers
        wstate[indices_of_infectives[i]] = 2 
    }
  } # end loop over infectives

  ##############################################################################
  # now plot the population, colored by their state
  ##############################################################################
  points(wx,wy,col=scol[(wstate+1)],cex=0.5)

  ##############################################################################
  # now increment the time, and count the number of infectives
  ##############################################################################
  time = time + delta_t
  ninf = length(wstate[wstate==1])
} # end while loop that checks if there is still at least one infected person



