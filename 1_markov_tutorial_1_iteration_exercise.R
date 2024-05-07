# install.packages("truncnorm")

rm(list = ls())
gc()

###### sample model parameters ######

# age at the start
n_age_init <- 25
# maximum age
n_age_max <- 110
# number of cycles (one year cycle)
n_t <- n_age_max - n_age_init
# the 4 health states of the model:
v_n <- c("H", "S1", "S2", "D") 
# number of health states 
n_states <- length(v_n)
# discount rate
d_r <- 0.035
#number of MCMC iterations 
n_sim <- 1000
#cost of treatment
c_Trt <- 50

###### Function for creating PSA inputs ######

#data frame with parameter values (n_sim random values for each parameter)
df_psa <- data.frame(
  
  # Transition probabilities (per cycle)
  
  # prob Healthy -> Sick ~ Beta(30, 170)
  p_HS1   = rbeta(n = n_sim,
                  shape1 =  30,
                  shape2 =  170),
  
  # prob Sick    -> Healthy ~ Beta(60, 60)
  # p_S1H
  
  # prob Sick    -> Sicker ~ Beta(84, 716)
  # p_S1S2       
  
  # prob Healthy -> Dead ~ Beta(10, 1990)
  # p_HD
  
  # rate ratio death S1 vs healthy
  hr_S1   = rlnorm(n = n_sim,
                   meanlog =  log(3),
                   sdlog =   0.01),
  
  # rate ratio death S2 vs healthy
  hr_S2   = rlnorm(n = n_sim,
                   meanlog =  log(10),
                   sdlog =  0.02),   
  
  # Cost vectors with length n_sim
  # cost p/cycle in state H
  c_H   = rgamma(n = n_sim, 
                 shape = 100, 
                 scale = 20),
  
  # cost p/cycle in state S1 ~ Gamma(shape = 177.8, scale = 22.5)
  #c_S1
  
  # cost p/cycle in state S2 ~ Gamma(shape = 225, scale = 66.7)
  #c_S2
  
  # cost p/cycle in state D
  c_D   = 0,
  
  # cost p/cycle of treatment
  c_Trt = c_Trt,
  
  # Utility vectors with length n_sim
  # utility when healthy (~N(1, 0.01) -Inf->1)
  u_H   = rtruncnorm(n = n_sim, 
                     mean = 1, 
                     sd = 0.01, 
                     b = 1),
  
  # utility when sick (~N(0.75, 0.02) -Inf->1)
  # u_S1
  
  # utility when sicker (~N(0.50, 0.03) -Inf->1)
  # u_S2
  
  # utility when dead
  u_D   = 0,
  
  # utility when being treated (~N(0.95, 0.02) -Inf->1)
  u_Trt = rtruncnorm(n = n_sim,
                     mean = 0.95,
                     sd = 0.02, 
                     b = 1)  
)

###### Step 2. Markov trace for ONE random iteration ######

df_psa_1 <- df_psa[1,]

# rate of death in healthy rate = - log(1 - probability of death in healthy)
# r_HD
# rate of death in sick = rate of death in healthy * hazard rate if sick
# r_S1D 	  
# rate of death in sicker = rate of death in healthy * hazard rate if sicker
# r_S2D
# probability of death in sick = 1 - exp(-rate of death in sick)
# p_S1D
# probability of death in sicker = 1 - exp(-rate of death in sicker)
# p_S2D
# calculate discount weight for each cycle
v_dwe <- v_dwc <- 1 / (1 + d_r) ^ (0:n_t) 

#transition probability matrix for NO treatment
m_P <- matrix(0,
              nrow = n_states, 
              ncol = n_states,
              dimnames = list(v_n, v_n))
# fill in the transition probability array

### From Healthy
m_P["H", "H"]  <- 1 - (df_psa_1$p_HS1 + df_psa_1$p_HD)
m_P["H", "S1"] <- df_psa_1$p_HS1
m_P["H", "D"]  <- df_psa_1$p_HD
### From Sick
m_P["S1", "H"]
m_P["S1", "S1"]
m_P["S1", "S2"]
m_P["S1", "D"]
### From Sicker
m_P["S2", "S2"]
m_P["S2", "D"]
### From Dead
m_P["D", "D"]

# create empty Markov trace 
m_TR <- matrix(data = NA, 
               nrow = n_t + 1, 
               ncol = n_states, 
               dimnames = list(0:n_t, v_n))     

# initialize Markov trace
m_TR[1, ] <- c(1, 0, 0, 0)   

# PROCESS (Markov model in cycles 1 to n_t)

for (t in 1:n_t){ # throughout the number of cycles
  # estimate next cycle (t+1) of Markov trace
  m_TR[t + 1, ] <- m_TR[t, ] %*% m_P           
}

# OUTPUTS

# create vectors of utility and costs for each state
v_u_trt    <- c(df_psa_1$u_H, df_psa_1$u_Trt, df_psa_1$u_S2, df_psa_1$u_D)
v_u_no_trt <- c(df_psa_1$u_H, df_psa_1$u_S1, df_psa_1$u_S2, df_psa_1$u_D)
v_c_trt
v_c_no_trt <- c(df_psa_1$c_H, df_psa_1$c_S1, df_psa_1$c_S2, df_psa_1$c_D)
# estimate mean QALYs and costs
v_E_no_trt <- m_TR %*% v_u_no_trt

v_E_trt    <- m_TR %*% v_u_trt

v_C_no_trt <- m_TR %*% v_c_no_trt

v_C_trt    <- m_TR %*% v_c_trt

### discount costs and QALYs
#   1x31 %*% 31x1 -> 1x1

te_no_trt <- t(v_E_no_trt) %*% v_dwe  
te_trt    <- t(v_E_trt) %*% v_dwe

tc_no_trt
tc_trt

results <- c(
  "Cost_NoTrt" = tc_no_trt, 
  "Cost_Trt"   = tc_trt, 
  "QALY_NoTrt" = te_no_trt, 
  "QALY_Trt" = te_trt,
  "ICER" = (tc_trt - tc_no_trt)/
    (te_trt - te_no_trt)
)

