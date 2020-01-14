###############################################################################
#
#      _____     _    ____________  __  ___
#     / ___/    | |  / / ____/ __ \/  |/  /
#     \__ \_____| | / / /   / / / / /|_/ /
#    ___/ /_____/ |/ / /___/ /_/ / /  / /
#   /____/      |___/\____/\____/_/  /_/
#
#   Run mosquito model
#   Sean Wu (slwu89@berkeley.edu)
#   January 2020
#
###############################################################################

rm(list=ls());gc()

# to compile the C++ version
library(Rcpp)

# plotting packages
library(here)
library(ggplot2)
library(reshape2)

# for running large numbers of stochastic runs in parallel
library(doSNOW)
library(parallel)
library(abind) # for abind

source(here::here("sim-src/mosquito-equilibrium.R"))
Rcpp::sourceCpp(here::here("sim-src/mosquito-both.cpp"))


###############################################################################
# Parameters
###############################################################################

## Model parameters:
theta <- list(
  ## Mosquito life cycle parameters:
  beta = 21.19, # Number of eggs laid per day by female mosquito
  muEL = 0.034, # Early larval instar daily mortality
  muLL = 0.035, # Late larval instar daily mortality
  muPL = 0.25, # Pupal daily mortality
  durEL = 6.64, # Duration of early instar stage
  durLL = 3.72, # Duration of late instar stage
  durPL = 0.64, # Duration of pupal stage
  durEV = 10, # Duration of latent period in mosquito (days)
  gamma = 13.25, # Effect of density-dependence on late instars relative to early instars
  tau1 = 0.68, # Time spent foraging for a blood meal at 0% ITN coverage
  tau2 = 2.32, # Time spent resting and ovipositing by a mosquito

  ## Intervention parameters (variable):
  ITNcov = 0.5, # ITN coverage
  IRScov = 0.25, # IRS coverave
  time_ITN_on = 250, # When ITNs are applied (days)
  time_IRS_on = 500, # When IRS is applied (days)

  ## Species-specific parameters:
  ## An. gambiae:
  muV = 1/7.6, # Adult mosquito daily mortality
  Q0 = 0.92, # Human blood index
  phiB = 0.89, # Proportion of bites on a person while they are in bed
  phiI = 0.97, # Proportion of bites on a person while they are indoors
  rITN = 0.56, # Probability of mosquito repeating a feeding attempt due to IRS
  sITN = 0.03, # Probability of mosquito feeding and surviving in presence of ITNs
  rIRS = 0.60, # Probability of mosquito repeating a feeding attempt due to IRS
  sIRS = 0, # Probability of mosquito feeding and surviving in presence of IRS

  ## An. arabiensis:
  # muV = 1/7.6, # Adult mosquito daily mortality
  # Q0 = 0.71, # Human blood index
  # phiB = 0.90, # Proportion of bites on a person while they are in bed
  # phiI = 0.96, # Proportion of bites on a person while they are indoors
  # rITN = 0.48, # Probability of mosquito repeating a feeding attempt due to IRS
  # sITN = 0.39, # Probability of mosquito feeding and surviving in presence of ITNs
  # rIRS = 0.60, # Probability of mosquito repeating a feeding attempt due to IRS
  # sIRS = 0, # Probability of mosquito feeding and surviving in presence of IRS

  ## An. funestus:
  # muV = 1/8.9, # Adult mosquito daily mortality
  # Q0 = 0.94, # Human blood index
  # phiB = 0.90, # Proportion of bites on a person while they are in bed
  # phiI = 0.98, # Proportion of bites on a person while they are indoors
  # rITN = 0.56, # Probability of mosquito repeating a feeding attempt due to IRS
  # sITN = 0.03, # Probability of mosquito feeding and surviving in presence of ITNs
  # rIRS = 0.63, # Probability of mosquito repeating a feeding attempt due to IRS
  # sIRS = 0, # Probability of mosquito feeding and surviving in presence of IRS

  ## Additional transmission parameters:
  f0 = 1/3, # Daily biting rate by mosquitoes on animals and humans
  epsilon0 = 10/365, # Daily entomological inolculation rate
  iH_eq = 0.45, # Equilibrium malaria prevalence in humans
  NH_eq = 200, # Equilibrium human population size
  bV = 0.15 # Probability of transmission from human to vector per infectious bite
)


###############################################################################
# deterministic model in R
###############################################################################

source(here::here("sim-src/mosquito-deterministic.R"))

# this is our time step
dt <- 0.5

# assume 500 infectious mosquitoes at equilibrium
IV_eq <- 500

# force of infection on mosquitoes at equilibrium
lambdaV <- (theta$Q0 * theta$f0) * theta$iH_eq * theta$bV
theta$lambdaV <- lambdaV

# use discrete-time equilibria
eq_dt <- calc_eq(theta = theta,dt = dt,IV = IV_eq,lambdaV = lambdaV)

# make the node
node <- make_node()
node$EL <- eq_dt$EL_eq
node$LL <- eq_dt$LL_eq
node$PL <- eq_dt$PL_eq
node$NM <- sum(eq_dt$SV_eq,eq_dt$EV_eq,IV_eq)
node$SV <- eq_dt$SV_eq
node$EV <- eq_dt$EV_eq
node$IV <- IV_eq

# add K to the list of parameters
theta$K <- eq_dt$K_eq

# run the model for tmax days
tmax <- 1e3
time <- seq(from=1,to=tmax,by=dt)

# model output is stored in the matrix sample_pop
sample_grid <- tsamp <- c(0,seq(from=10,to = tmax,by = 1))
sample_pop <- matrix(0,nrow=6,ncol=length(sample_grid),dimnames=list(c("EL","LL","PL","SV","EV","IV"),paste0(sample_grid)))

sample_pop["EL",1] <- node$EL
sample_pop["LL",1] <- node$LL
sample_pop["PL",1] <- node$PL
sample_pop["SV",1] <- node$SV
sample_pop["EV",1] <- node$EV
sample_pop["IV",1] <- node$IV
sample_grid <- sample_grid[-1]

# run simulation
pb <- txtProgressBar(min = 1,max = length(time))
for(t in 1:length(time)){

  # euler step
  euler_step(node = node,pars = theta,tnow = time[t],dt = dt)

  # sample the population (done at the very end of the time-step, because its not part of the dynamics)
  if(time[t] == sample_grid[1]){

    sample_pop["EL",as.character(sample_grid[1])] <- node$EL
    sample_pop["LL",as.character(sample_grid[1])] <- node$LL
    sample_pop["PL",as.character(sample_grid[1])] <- node$PL
    sample_pop["SV",as.character(sample_grid[1])] <- node$SV
    sample_pop["EV",as.character(sample_grid[1])] <- node$EV
    sample_pop["IV",as.character(sample_grid[1])] <- node$IV
    sample_grid <- sample_grid[-1]

  }
  setTxtProgressBar(pb = pb,value = t)
}

data_det_R <- data.frame(time=tsamp,SV=sample_pop["SV",],EV=sample_pop["EV",],IV=sample_pop["IV",])
data_det_R <- melt(data_det_R,"time")

ggplot(data=data_det_R) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()


###############################################################################
# deterministic model in C++
###############################################################################

# the C++ code needs the intervention parameters in a seperate list
int_pars <- list(
  ITNcov = 0.5,
  IRScov = 0.25,
  time_ITN_on = 250,
  time_IRS_on = 500
)

data_det_C <- cpp_deterministic(
  time = time,
  dt = dt,
  EL = eq_dt$EL_eq,
  LL = eq_dt$LL_eq,
  PL = eq_dt$PL_eq,
  SV = eq_dt$SV_eq,
  EV = eq_dt$EV_eq,
  IV = IV_eq,
  K_ = eq_dt$K_eq,
  pars_ = theta,
  int_pars_ = int_pars
)

data_det_C <- as.data.frame(data_det_C[,4:6])
data_det_C <- cbind(time=time,data_det_C)
data_det_C <- melt(data_det_C,"time")

ggplot(data=data_det_C) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()


###############################################################################
# stochastic model in R
###############################################################################

source(here::here("sim-src/mosquito-stochastic.R"))

# because we are using the stochastic model, we will run 100 simulations and plot averages/quantiles

# ensemble run parameters
nruns <- 100

# make a parallel cluster with 4 cores
cl <- makeSOCKcluster(4)
registerDoSNOW(cl)

# combine each matrix (sweep over cells) into a slice of a 3d array
acomb <- function(...) abind(..., along=3)

# progress bar
pb <- txtProgressBar(max = nruns, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

# we want output every day
tsamp <- c(0,seq(from=1,to = tmax,by = 1))

# array's 3rd dimension is over runs
sample_pop_dt <- foreach(i = 1:nruns, .combine = "acomb", .options.snow=opts,.packages=c("foreach")) %dopar% {

  # output
  sample_grid <- tsamp
  sample_pop <- matrix(0,nrow=6,ncol=length(sample_grid),dimnames=list(c("EL","LL","PL","SV","EV","IV"),paste0(sample_grid)))

  # make the node
  node <- make_node()
  node$EL <- as.integer(eq_dt$EL_eq)
  node$LL <- as.integer(eq_dt$LL_eq)
  node$PL <- as.integer(eq_dt$PL_eq)
  node$SV <- as.integer(eq_dt$SV_eq)
  node$EV <- as.integer(eq_dt$EV_eq)
  node$IV <- as.integer(IV_eq)

  # record output
  sample_pop["EL",1] <- node$EL
  sample_pop["LL",1] <- node$LL
  sample_pop["PL",1] <- node$PL
  sample_pop["SV",1] <- node$SV
  sample_pop["EV",1] <- node$EV
  sample_pop["IV",1] <- node$IV
  sample_grid <- sample_grid[-1]

  # run simulation
  for(t in 1:length(time)){

    # euler step
    euler_step(node = node,pars = theta,tnow = time[t],dt = dt)

    # sample the population (done at the very end of the time-step, because its not part of the dynamics)
    if(time[t] == sample_grid[1]){

      sample_pop["EL",as.character(sample_grid[1])] <- node$EL
      sample_pop["LL",as.character(sample_grid[1])] <- node$LL
      sample_pop["PL",as.character(sample_grid[1])] <- node$PL
      sample_pop["SV",as.character(sample_grid[1])] <- node$SV
      sample_pop["EV",as.character(sample_grid[1])] <- node$EV
      sample_pop["IV",as.character(sample_grid[1])] <- node$IV
      sample_grid <- sample_grid[-1]

    }
  }

  sample_pop
}

close(pb)
stopCluster(cl);rm(cl);gc()

# we want to plot the mean and 95% quantiles

# plot the output (correct equilibrium)
mean_SV_dt <- rowMeans(sample_pop_dt["SV",,])
mean_EV_dt <- rowMeans(sample_pop_dt["EV",,])
mean_IV_dt <- rowMeans(sample_pop_dt["IV",,])

traj_SV_dt <- melt(sample_pop_dt["SV",,])
colnames(traj_SV_dt) <- c("time","run","count")
traj_EV_dt <- melt(sample_pop_dt["EV",,])
colnames(traj_EV_dt) <- c("time","run","count")
traj_IV_dt <- melt(sample_pop_dt["IV",,])
colnames(traj_IV_dt) <- c("time","run","count")

quant_95 <- c(0.025,0.975)
quant_SV_dt <- apply(X = sample_pop_dt["SV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = quant_95)
})
quant_EV_dt <- apply(X = sample_pop_dt["EV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = quant_95)
})
quant_IV_dt <- apply(X = sample_pop_dt["IV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = quant_95)
})

plot_datSV_dt <- data.frame(time=tsamp,SV=mean_SV_dt,SV_l=quant_SV_dt[1,],SV_h=quant_SV_dt[2,])
plot_datEV_dt <- data.frame(time=tsamp,EV=mean_EV_dt,EV_l=quant_EV_dt[1,],EV_h=quant_EV_dt[2,])
plot_datIV_dt <- data.frame(time=tsamp,IV=mean_IV_dt,IV_l=quant_IV_dt[1,],IV_h=quant_IV_dt[2,])

# plot mean and 95% quantiles
ggplot() +
  geom_line(data=plot_datSV_dt,aes(x=time,y=SV),color="darkorchid2") +
  geom_ribbon(data=plot_datSV_dt,aes(x=time,ymin=SV_l,ymax=SV_h),alpha=0.35,fill="darkorchid2") +
  geom_line(data=plot_datEV_dt,aes(x=time,y=EV),color="dodgerblue2") +
  geom_ribbon(data=plot_datEV_dt,aes(x=time,ymin=EV_l,ymax=EV_h),alpha=0.35,fill="dodgerblue2") +
  geom_line(data=plot_datIV_dt,aes(x=time,y=IV),color="firebrick2") +
  geom_ribbon(data=plot_datIV_dt,aes(x=time,ymin=IV_l,ymax=IV_h),alpha=0.35,fill="firebrick2") +
  ylab("Counts") +
  xlab("Time") +
  theme_bw()

# plot individual simulation trajectores (takes a while)
ggplot() +
  geom_line(data=traj_SV_dt,aes(x=time,y=count),color="darkorchid2",alpha=0.25) +
  geom_line(data=traj_EV_dt,aes(x=time,y=count),color="dodgerblue2",alpha=0.25) +
  geom_line(data=traj_IV_dt,aes(x=time,y=count),color="firebrick2",alpha=0.25) +
  theme_bw()




###############################################################################
# stochastic model in C++
###############################################################################











###############################################################################
# stochastic approximation
###############################################################################

source(here::here("sim-src/mosquito-stochastic.R"))

# run model against continuous-time equilibria
with(theta,{

  a0 <<- Q0*f0 # Human biting rate at equilibrium
  lambdaV <<- a0*iH_eq*bV # Force of infection in mosquitoes at equilibrium

  IV_eq <<- 500
  EV_eq <<- durEV*IV_eq*muV
  SV_eq <<- ((durEV*IV_eq*muV) + ((durEV^2)*IV_eq*(muV^2))) / (durEV*lambdaV)
  NV_eq <<- IV_eq + EV_eq + SV_eq

  omega_1 <<- -0.5 * ( (gamma*(muLL/muEL)) - (durEL/durLL) + ((gamma-1)*muLL*durEL) )
  omega_2 <<- sqrt(
    (0.25 * (( (gamma*(muLL/muEL)) - (durEL/durLL) + ((gamma-1)*muLL*durEL) )^2)) + (gamma * ((beta*muLL*durEL) / (2*muEL*muV*durLL*(1 + (durPL*muPL)))) )
  )
  omega <<- omega_1 + omega_2

  PL_eq <<- 2*durPL*muV*NV_eq
  LL_eq <<- 2*muV*durLL*(1 + (durPL*muPL))*NV_eq
  EL_eq <<- 2*omega*muV*durLL*(1 + (durPL*muPL))*NV_eq

  K <<- (NV_eq*2*durLL*muV*(1 + (durPL*muPL))*gamma*(omega+1)) / ((omega/(muLL*durEL)) - (1/(muLL*durLL)) - 1)
})

theta_eq <- c(theta,K=K,lambdaV=lambdaV)

# ensemble run parameters
nruns <- 100
tmax <- 1e3
dt <- 1
time <- seq(from=1,to=tmax,by=dt)

# sampling grid
tsamp <- c(0,seq(from=1,to = tmax,by = 1))

# run ensemble of stochastic simulations
cl <- makeSOCKcluster(4)
registerDoSNOW(cl)

# combine each matrix (sweep over cells) into a slice of a 3d array
acomb <- function(...) abind(..., along=3)

# progress bar
pb <- txtProgressBar(max = nruns, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

# array's 3rd dimension is over runs
sample_pop_ct <- foreach(i = 1:nruns, .combine = "acomb", .options.snow=opts,.packages=c("foreach")) %dopar% {

  # output
  sample_grid <- tsamp
  sample_pop <- matrix(0,nrow=6,ncol=length(sample_grid),dimnames=list(c("EL","LL","PL","SV","EV","IV"),paste0(sample_grid)))

  # make the node
  node <- make_node()
  node$EL <- as.integer(EL_eq)
  node$LL <- as.integer(LL_eq)
  node$PL <- as.integer(PL_eq)
  node$SV <- as.integer(SV_eq)
  node$EV <- as.integer(EV_eq)
  node$IV <- as.integer(IV_eq)

  # record output
  sample_pop["EL",1] <- node$EL
  sample_pop["LL",1] <- node$LL
  sample_pop["PL",1] <- node$PL
  sample_pop["SV",1] <- node$SV
  sample_pop["EV",1] <- node$EV
  sample_pop["IV",1] <- node$IV
  sample_grid <- sample_grid[-1]

  # run simulation
  for(t in 1:length(time)){

    # euler step
    euler_step(node = node,pars = theta_eq,tnow = time[t],dt = dt)

    # sample the population (done at the very end of the time-step, because its not part of the dynamics)
    if(time[t] == sample_grid[1]){

      sample_pop["EL",as.character(sample_grid[1])] <- node$EL
      sample_pop["LL",as.character(sample_grid[1])] <- node$LL
      sample_pop["PL",as.character(sample_grid[1])] <- node$PL
      sample_pop["SV",as.character(sample_grid[1])] <- node$SV
      sample_pop["EV",as.character(sample_grid[1])] <- node$EV
      sample_pop["IV",as.character(sample_grid[1])] <- node$IV
      sample_grid <- sample_grid[-1]

    }
  }

  sample_pop
}

close(pb)
stopCluster(cl);rm(cl);gc()

# plot the output (incorrect equilibrium)
mean_SV_ct <- rowMeans(sample_pop_ct["SV",,])
mean_EV_ct <- rowMeans(sample_pop_ct["EV",,])
mean_IV_ct <- rowMeans(sample_pop_ct["IV",,])

quant_SV_ct <- apply(X = sample_pop_ct["SV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})
quant_EV_ct <- apply(X = sample_pop_ct["EV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})
quant_IV_ct <- apply(X = sample_pop_ct["IV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})

plot_datEV_ct <- data.frame(time=tsamp,EV=mean_EV_ct,EV_l=quant_EV_ct[1,],EV_h=quant_EV_ct[2,])
plot_datIV_ct <- data.frame(time=tsamp,IV=mean_IV_ct,IV_l=quant_IV_ct[1,],IV_h=quant_IV_ct[2,])

# ggplot() +
#   geom_line(data=plot_datEV_ct,aes(x=time,y=EV),color="steelblue") +
#   geom_ribbon(data=plot_datEV_ct,aes(x=time,ymin=EV_l,ymax=EV_h),alpha=0.35,fill="steelblue") +
#   geom_line(data=plot_datIV_ct,aes(x=time,y=IV),color="firebrick3") +
#   geom_ribbon(data=plot_datIV_ct,aes(x=time,ymin=IV_l,ymax=IV_h),alpha=0.35,fill="firebrick3") +
#   theme_bw()


# use discrete-time equilibria
eq_dt <- calc_eq(theta = theta_eq,dt = dt,IV = IV_eq,lambdaV = lambdaV)

theta_eq$K <- eq_dt$K_eq

# run ensemble of stochastic simulations
cl <- makeSOCKcluster(4)
registerDoSNOW(cl)

# combine each matrix (sweep over cells) into a slice of a 3d array
acomb <- function(...) abind(..., along=3)

# progress bar
pb <- txtProgressBar(max = nruns, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

# array's 3rd dimension is over runs
sample_pop_dt <- foreach(i = 1:nruns, .combine = "acomb", .options.snow=opts,.packages=c("foreach")) %dopar% {

  # output
  sample_grid <- tsamp
  sample_pop <- matrix(0,nrow=6,ncol=length(sample_grid),dimnames=list(c("EL","LL","PL","SV","EV","IV"),paste0(sample_grid)))

  # make the node
  node <- make_node()
  node$EL <- as.integer(eq_dt$EL_eq)
  node$LL <- as.integer(eq_dt$LL_eq)
  node$PL <- as.integer(eq_dt$PL_eq)
  node$SV <- as.integer(eq_dt$SV_eq)
  node$EV <- as.integer(eq_dt$EV_eq)
  node$IV <- as.integer(IV_eq)

  # record output
  sample_pop["EL",1] <- node$EL
  sample_pop["LL",1] <- node$LL
  sample_pop["PL",1] <- node$PL
  sample_pop["SV",1] <- node$SV
  sample_pop["EV",1] <- node$EV
  sample_pop["IV",1] <- node$IV
  sample_grid <- sample_grid[-1]

  # run simulation
  for(t in 1:length(time)){

    # euler step
    euler_step(node = node,pars = theta_eq,tnow = time[t],dt = dt)

    # sample the population (done at the very end of the time-step, because its not part of the dynamics)
    if(time[t] == sample_grid[1]){

      sample_pop["EL",as.character(sample_grid[1])] <- node$EL
      sample_pop["LL",as.character(sample_grid[1])] <- node$LL
      sample_pop["PL",as.character(sample_grid[1])] <- node$PL
      sample_pop["SV",as.character(sample_grid[1])] <- node$SV
      sample_pop["EV",as.character(sample_grid[1])] <- node$EV
      sample_pop["IV",as.character(sample_grid[1])] <- node$IV
      sample_grid <- sample_grid[-1]

    }
  }

  sample_pop
}

close(pb)
stopCluster(cl);rm(cl);gc()

# plot the output (correct equilibrium)
mean_SV_dt <- rowMeans(sample_pop_dt["SV",,])
mean_EV_dt <- rowMeans(sample_pop_dt["EV",,])
mean_IV_dt <- rowMeans(sample_pop_dt["IV",,])

traj_EV_dt <- melt(sample_pop_dt["EV",,])
colnames(traj_EV_dt) <- c("time","run","count")
traj_IV_dt <- melt(sample_pop_dt["IV",,])
colnames(traj_IV_dt) <- c("time","run","count")

quant_SV_dt <- apply(X = sample_pop_dt["SV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})
quant_EV_dt <- apply(X = sample_pop_dt["EV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})
quant_IV_dt <- apply(X = sample_pop_dt["IV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})

plot_datEV_dt <- data.frame(time=tsamp,EV=mean_EV_dt,EV_l=quant_EV_dt[1,],EV_h=quant_EV_dt[2,])
plot_datIV_dt <- data.frame(time=tsamp,IV=mean_IV_dt,IV_l=quant_IV_dt[1,],IV_h=quant_IV_dt[2,])

ggplot() +
  # incorrect
  geom_line(data=plot_datEV_ct,aes(x=time,y=EV),color="dodgerblue4") +
  geom_ribbon(data=plot_datEV_ct,aes(x=time,ymin=EV_l,ymax=EV_h),alpha=0.35,fill="dodgerblue4") +
  geom_line(data=plot_datIV_ct,aes(x=time,y=IV),color="firebrick4") +
  geom_ribbon(data=plot_datIV_ct,aes(x=time,ymin=IV_l,ymax=IV_h),alpha=0.35,fill="firebrick4") +
  # correct
  geom_line(data=plot_datEV_dt,aes(x=time,y=EV),color="dodgerblue2") +
  geom_ribbon(data=plot_datEV_dt,aes(x=time,ymin=EV_l,ymax=EV_h),alpha=0.35,fill="dodgerblue2") +
  geom_line(data=plot_datIV_dt,aes(x=time,y=IV),color="firebrick2") +
  geom_ribbon(data=plot_datIV_dt,aes(x=time,ymin=IV_l,ymax=IV_h),alpha=0.35,fill="firebrick2") +
  ylab("Counts") +
  xlab("Time") +
  theme_bw()

# show individual trajectories
ggplot() +
  geom_line(data=plot_datEV_dt,aes(x=time,y=EV),color="dodgerblue2") +
  geom_line(data=plot_datIV_dt,aes(x=time,y=IV),color="firebrick2") +
  geom_line(data=traj_EV_dt,aes(x=time,y=count),color="dodgerblue2",alpha=0.25) +
  geom_line(data=traj_IV_dt,aes(x=time,y=count),color="firebrick2",alpha=0.25) +
  theme_bw()


# test C++ version

int_pars <- list(
  ITNcov = 0.5, # ITN coverage
  IRScov = 0.25, # IRS coverave
  time_ITN_on = 1e4, # When ITNs are applied (days)
  time_IRS_on = 1e4 # When IRS is applied (days)
)


sample_pop_cpp <- test_stochastic(time = time,dt = dt,
                              EL_ = as.integer(eq_dt$EL_eq),
                              LL_ = as.integer(eq_dt$LL_eq),
                              PL_ = as.integer(eq_dt$PL_eq),
                              SV_ = as.integer(eq_dt$SV_eq),
                              EV_ = as.integer(eq_dt$EV_eq),
                              IV_ = as.integer(IV_eq),
                              K_ = eq_dt$K_eq,
                              pars_ = theta_eq,int_pars_ = int_pars)


sample_pop_cpp <- as.data.frame(sample_pop_cpp[,4:6])
sample_pop_cpp <- cbind(time=time,sample_pop_cpp)
sample_pop_cpp <- melt(sample_pop_cpp,"time")

ggplot(data=sample_pop_cpp) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()


# run ensemble of stochastic simulations
# combine each matrix (sweep over cells) into a slice of a 3d array
acomb <- function(...) abind(..., along=3)

# progress bar
pb <- txtProgressBar(max = nruns, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

# array's 3rd dimension is over runs
sample_pop_cpp <- foreach(i = 1:nruns, .combine = "acomb", .options.snow=opts,.packages=c("foreach")) %do% {

  sample_pop <- test_stochastic(time = time,dt = dt,
                                EL_ = as.integer(eq_dt$EL_eq),
                                LL_ = as.integer(eq_dt$LL_eq),
                                PL_ = as.integer(eq_dt$PL_eq),
                                SV_ = as.integer(eq_dt$SV_eq),
                                EV_ = as.integer(eq_dt$EV_eq),
                                IV_ = as.integer(IV_eq),
                                K_ = eq_dt$K_eq,
                                pars_ = theta_eq,int_pars_ = int_pars)

  sample_pop
}

# plot the output (correct equilibrium)
mean_SV_cpp <- rowMeans(sample_pop_cpp[,"SV",])
mean_EV_cpp <- rowMeans(sample_pop_cpp[,"EV",])
mean_IV_cpp <- rowMeans(sample_pop_cpp[,"IV",])

traj_SV_cpp <- melt(sample_pop_cpp[,"SV",])
colnames(traj_SV_cpp) <- c("time","run","count")
traj_EV_cpp <- melt(sample_pop_cpp[,"EV",])
colnames(traj_EV_cpp) <- c("time","run","count")
traj_IV_cpp <- melt(sample_pop_cpp[,"IV",])
colnames(traj_IV_cpp) <- c("time","run","count")

quant_SV_cpp <- apply(X = sample_pop_cpp[,"SV",],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})
quant_EV_cpp <- apply(X = sample_pop_cpp[,"EV",],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})
quant_IV_cpp <- apply(X = sample_pop_cpp[,"IV",],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})

plot_datEV_cpp <- data.frame(time=time,EV=mean_EV_cpp,EV_l=quant_EV_cpp[1,],EV_h=quant_EV_cpp[2,])
plot_datIV_cpp <- data.frame(time=time,IV=mean_IV_cpp,IV_l=quant_IV_cpp[1,],IV_h=quant_IV_cpp[2,])

ggplot() +
  # C++
  geom_line(data=plot_datEV_cpp,aes(x=time,y=EV),color="darkorchid4") +
  geom_ribbon(data=plot_datEV_cpp,aes(x=time,ymin=EV_l,ymax=EV_h),alpha=0.35,fill="darkorchid4") +
  geom_line(data=plot_datIV_cpp,aes(x=time,y=IV),color="firebrick4") +
  geom_ribbon(data=plot_datIV_cpp,aes(x=time,ymin=IV_l,ymax=IV_h),alpha=0.35,fill="firebrick4") +
  # R
  geom_line(data=plot_datEV_dt,aes(x=time,y=EV),color="darkorchid2") +
  geom_ribbon(data=plot_datEV_dt,aes(x=time,ymin=EV_l,ymax=EV_h),alpha=0.35,fill="darkorchid2") +
  geom_line(data=plot_datIV_dt,aes(x=time,y=IV),color="firebrick2") +
  geom_ribbon(data=plot_datIV_dt,aes(x=time,ymin=IV_l,ymax=IV_h),alpha=0.35,fill="firebrick2") +
  ylab("Counts") +
  xlab("Time") +
  theme_bw()

ggplot() +
  # geom_line(data=traj_EV_dt,aes(x=time,y=count),color="dodgerblue2",alpha=0.25) +
  # geom_line(data=traj_IV_dt,aes(x=time,y=count),color="firebrick2",alpha=0.25) +
  geom_line(data=traj_SV_cpp,aes(x=time,y=count,group=run),linetype=1,color="mediumpurple3",alpha=0.1) +
  geom_line(data=traj_EV_cpp,aes(x=time,y=count,group=run),linetype=1,color="dodgerblue4",alpha=0.1) +
  geom_line(data=traj_IV_cpp,aes(x=time,y=count,group=run),linetype=1,color="firebrick4",alpha=0.1) +
  theme_bw()

library(tikzDevice)

tikz(file = "svcom.tex", width = 12, height = 8)

ggplot() +
  # C++
  geom_line(data=plot_datEV_cpp,aes(x=time,y=EV),color="darkorchid4") +
  geom_ribbon(data=plot_datEV_cpp,aes(x=time,ymin=EV_l,ymax=EV_h),alpha=0.35,fill="darkorchid4") +
  geom_line(data=plot_datIV_cpp,aes(x=time,y=IV),color="firebrick4") +
  geom_ribbon(data=plot_datIV_cpp,aes(x=time,ymin=IV_l,ymax=IV_h),alpha=0.35,fill="firebrick4") +
  # R
  geom_line(data=plot_datEV_dt,aes(x=time,y=EV),color="darkorchid2") +
  geom_ribbon(data=plot_datEV_dt,aes(x=time,ymin=EV_l,ymax=EV_h),alpha=0.35,fill="darkorchid2") +
  geom_line(data=plot_datIV_dt,aes(x=time,y=IV),color="firebrick2") +
  geom_ribbon(data=plot_datIV_dt,aes(x=time,ymin=IV_l,ymax=IV_h),alpha=0.35,fill="firebrick2") +
  ylab("Counts") +
  xlab("Time") +
  theme_bw()

dev.off()
