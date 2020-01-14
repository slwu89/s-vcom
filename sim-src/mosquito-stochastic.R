################################################################################
#   Exponential distribution of EIP
#   Stochastic Model
################################################################################

# make a single node
make_node <- function(){
  
  node <- list()
  
  # new eggs are generated from a conditionally independent Poisson process
  node$EL_new <- 0
  storage.mode(node$EL_new) <- "integer"
  
  # probabilities & transitions for early-stage instars (names are destinations of jumps)
  node$EL_probs <- rep(0,3)
  node$EL_probs <- setNames(node$EL_probs,c("EL","D","LL"))
  node$EL_transitions <- rep(0,3)
  node$EL_transitions <- setNames(node$EL_transitions,c("EL","D","LL"))
  storage.mode(node$EL_transitions) <- "integer"
  
  # probabilities & transitions for late-stage instars
  node$LL_probs <- rep(0,3)
  node$LL_probs <- setNames(node$LL_probs,c("LL","D","PL"))
  node$LL_transitions <- rep(0,3)
  node$LL_transitions <- setNames(node$LL_transitions,c("LL","D","PL"))
  storage.mode(node$LL_transitions) <- "integer"
  
  # probabilities & transitions for pupae
  node$PL_probs <- rep(0,4)
  node$PL_probs <- setNames(node$PL_probs,c("PL","D","SV_F","SV_M"))
  node$PL_transitions <- rep(0,4)
  node$PL_transitions <- setNames(node$PL_transitions,c("PL","D","SV_F","SV_M"))
  storage.mode(node$PL_transitions) <- "integer"
  
  # probabilities & transitions for susceptible vectors
  node$SV_probs <- rep(0,3)
  node$SV_probs <- setNames(node$SV_probs,c("SV","D","EV"))
  node$SV_transitions <- rep(0,3)
  node$SV_transitions <- setNames(node$SV_transitions,c("SV","D","EV"))
  storage.mode(node$SV_transitions) <- "integer"
  
  # probabilities & transitions for incubating vectors
  node$EV_probs <- rep(0,3)
  node$EV_probs <- setNames(node$EV_probs,c("EV","D","IV"))
  node$EV_transitions <- rep(0,3)
  node$EV_transitions <- setNames(node$EV_transitions,c("EV","D","IV"))
  storage.mode(node$EV_transitions) <- "integer"
  
  # probabilities & transitions for infectious vectors
  node$IV_probs <- rep(0,2)
  node$IV_probs <- setNames(node$IV_probs,c("IV","D"))
  node$IV_transitions <- rep(0,2)
  node$IV_transitions <- setNames(node$IV_transitions,c("IV","D"))
  storage.mode(node$IV_transitions) <- "integer"
  
  # state space
  node$EL <- 0L
  node$LL <- 0L
  node$PL <- 0L
  node$SV <- 0L
  node$EV <- 0L
  node$IV <- 0L
  
  list2env(node,hash=TRUE)
}

# execute a stochastic Euler-step for the model
# rates -> probabilities are defined by discretisation of the underlying
# CTMC system (competing hazards)
euler_step <- function(node,pars,tnow,dt){
  with(pars,{
    
    ########################################
    # INTERVENTION-DEPENDENT PARAMETERS
    ########################################
    
    delta <- 1/(tau1+tau2) # Inverse of gonotrophic cycle without ITNs/IRS
    e_ov <- beta*(exp(muV/delta)-1)/muV # Number of eggs per oviposition per mosquito
    
    ## Derived parameters which depend on intervention status:
    if(tnow > pars$time_ITN_on){
      ITNcov_t <- pars$ITNcov
    } else {
      ITNcov_t <- 0
    }
    if(tnow > pars$time_IRS_on){
      IRScov_t <- pars$IRScov
    } else {
      IRScov_t <- 0
    }
    
    # zCom: Probability of a mosquito being repelled from an ITN or IRS-treated house:
    c0 <- 1 - ITNcov_t - IRScov_t + ITNcov_t*IRScov_t
    cITN <- ITNcov_t - ITNcov_t*IRScov_t
    cIRS <- IRScov_t - ITNcov_t*IRScov_t
    cCom <- ITNcov_t*IRScov_t
    rCom <- rIRS + (1-rIRS)*rITN
    sCom  <- (1-rIRS)*sITN*sIRS
    zCom <- Q0*cITN*phiB*rITN + Q0*cIRS*phiI*rIRS + Q0*cCom*(phiI-phiB)*rIRS + Q0*cCom*phiB*rCom
    
    # deltaCom: Inverse of gonotrophic cycle length with ITNs & IRS:
    deltaCom <- 1/(tau1/(1-zCom) + tau2)
    
    # wCom: Probability that a surviving mosquito succeeds in feeding during a single attempt:
    wCom <- 1 - Q0 + Q0*c0 + Q0*cITN*(1-phiB+phiB*sITN) + Q0*cIRS*(1-phiI+phiI*sIRS) + Q0*cCom*((phiI-phiB)*sIRS + 1-phiI + phiB*sCom)
    
    # muVCom: Female mosquito death rate in presence of ITNs & IRS:
    p10 <- exp(-muV*tau1)
    p1Com <- p10*wCom/(1 - zCom*p10)
    p2 <- exp(-muV*tau2)
    pCom <- (p1Com*p2)^deltaCom
    muVCom <- -log(pCom)
    
    # betaCom: Eggs laid per day by female mosquitoes in presence of ITNs & IRS:
    betaCom <- e_ov*muVCom/(exp(muVCom/deltaCom) - 1)
    
    ########################################
    # EARLY-STAGE LARVAL INSTARS (EL)
    ########################################
    
    # inbound oviposition to EL
    NV <- sum(node$SV,node$EV,node$IV)
    node$EL_new <- rpois(n = 1,lambda = betaCom*NV*dt)
    
    # instantaneous hazards for EL
    haz_EL_mort <- muEL*(1 + ((node$EL+node$LL)/K))
    haz_EL_2LL <- 1/durEL
    
    # jump probabilities
    node$EL_probs[["EL"]] <- exp(-(haz_EL_mort + haz_EL_2LL)*dt)
    node$EL_probs[["D"]] <- (1 - node$EL_probs[["EL"]])*(haz_EL_mort / (haz_EL_mort + haz_EL_2LL)) # death
    node$EL_probs[["LL"]] <- (1 - node$EL_probs[["EL"]])*(haz_EL_2LL / (haz_EL_mort + haz_EL_2LL)) # to late-instar
    
    # jump sizes
    node$EL_transitions[c("EL","D","LL")] <- as.vector(rmultinom(n=1,size=node$EL,prob=node$EL_probs))
    
    ########################################
    # LATE-STAGE LARVAL INSTARS (LL)
    ########################################
    
    # instantaneous hazards for LL
    haz_LL_mort <- muLL*(1 + gamma*((node$EL+node$LL)/K))
    haz_LL_2PL <- 1/durLL
    
    # jump probabilities
    node$LL_probs[["LL"]] <- exp(-(haz_LL_mort + haz_LL_2PL)*dt)
    node$LL_probs[["D"]] <- (1 - node$LL_probs[["LL"]])*(haz_LL_mort / (haz_LL_mort + haz_LL_2PL)) # death
    node$LL_probs[["PL"]] <- (1 - node$LL_probs[["LL"]])*(haz_LL_2PL / (haz_LL_mort + haz_LL_2PL)) # to pupae
    
    # jump sizes
    node$LL_transitions[c("LL","D","PL")] <- as.vector(rmultinom(n=1,size=node$LL,prob=node$LL_probs))
    
    ########################################
    # PUPAE (PL)
    ########################################
    
    # instantaneous hazards for PL
    haz_PL_mort <- muPL
    haz_PL_2SV_F <- (1/durPL)*0.5
    haz_PL_2SV_M <- (1/durPL)*0.5
    
    # jump probabilities
    node$PL_probs[["PL"]] <- exp(-(haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)*dt)
    node$PL_probs[["D"]] <- (1 - node$PL_probs[["PL"]])*(haz_PL_mort / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)) # death
    node$PL_probs[["SV_F"]] <- (1 - node$PL_probs[["PL"]])*(haz_PL_2SV_F / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)) # to susceptible female
    node$PL_probs[["SV_M"]] <- (1 - node$PL_probs[["PL"]])*(haz_PL_2SV_M / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)) # to susceptible males
    
    # jump sizes
    node$PL_transitions[c("PL","D","SV_F","SV_M")] <- as.vector(rmultinom(n=1,size=node$PL,prob=node$PL_probs))
    
    ########################################
    # SUSCEPTIBLE VECTORS (SV)
    ########################################
    
    # instantaneous hazards for SV
    haz_SV_mort <-  muVCom
    haz_SV_inf <- lambdaV
    
    # jump probabilities
    node$SV_probs[["SV"]] <- exp(-(haz_SV_mort + haz_SV_inf)*dt)
    node$SV_probs[["D"]] <- (1 - node$SV_probs[["SV"]])*(haz_SV_mort / (haz_SV_mort + haz_SV_inf)) # death
    node$SV_probs[["EV"]] <- (1 - node$SV_probs[["SV"]])*(haz_SV_inf / (haz_SV_mort + haz_SV_inf)) # to incubating
    
    # jump sizes
    node$SV_transitions[c("SV","D","EV")] <- as.vector(rmultinom(n=1,size=node$SV,prob=node$SV_probs))
    
    ########################################
    # INCUBATING VECTORS (EV)
    ########################################
    
    # instantaneous hazards for EV
    haz_EV_mort <-  muVCom
    haz_EV_inc <- 1/durEV
    
    # jump probabilities
    node$EV_probs[["EV"]] <- exp(-(haz_EV_mort + haz_EV_inc)*dt)
    node$EV_probs[["D"]] <- (1 - node$EV_probs[["EV"]])*(haz_EV_mort / (haz_EV_mort + haz_EV_inc)) # death
    node$EV_probs[["IV"]] <- (1 - node$EV_probs[["EV"]])*(haz_EV_inc / (haz_EV_mort + haz_EV_inc)) # to infectious
    
    # jump sizes
    node$EV_transitions[c("EV","D","IV")] <- as.vector(rmultinom(n=1,size=node$EV,prob=node$EV_probs))
    
    ########################################
    # INFECTIOUS VECTORS (IV)
    ########################################
    
    # instantaneous hazards for IV
    haz_IV_mort <- muVCom
    
    # jump probabilities
    node$IV_probs[["IV"]] <- exp(-haz_IV_mort*dt)
    node$IV_probs[["D"]] <- (1 - node$IV_probs[["IV"]])
    
    # jump sizes
    node$IV_transitions[c("IV","D")] <- as.vector(rmultinom(n=1,size=node$IV,prob=node$IV_probs))
    
    ########################################
    # UPDATE POPULATION
    ########################################
    
    node$EL <- node$EL_transitions[["EL"]] + node$EL_new
    node$LL <- node$LL_transitions[["LL"]] + node$EL_transitions[["LL"]]
    node$PL <- node$PL_transitions[["PL"]] + node$LL_transitions[["PL"]]
    node$SV <- node$SV_transitions[["SV"]] + node$PL_transitions[["SV_F"]]
    node$EV <- node$EV_transitions[["EV"]] + node$SV_transitions[["EV"]]
    node$IV <- node$IV_transitions[["IV"]] + node$EV_transitions[["IV"]]
    
  })
}