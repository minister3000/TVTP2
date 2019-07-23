

#******************************************************************
#*****                  Dissertation Effort                ********
#*****  Spread Dynamics proxied with Regime switch models  ********
#*****  with Time Varying Transition Probabilities         ********
#*****                Markov TVTP Estimation               ********
#******************************************************************


#******************************************************************
#                 Data Import, Variable Declaration
#******************************************************************


options(max.print = 100000)

# load package 'TVTP'. This R package has been written exclusively
# to accompany this dissertation effort by providing the neccessary
# data and algorithms to estimate unknown model parameters by maxi-
# mizing the incomplete data log-likelihood using the Expectation
# Maximization - Maximum Likelihood (EM-ML) algorithm

library(TVTP) # Make sure package 'TVTP' is installed

# load data frame with monthly time series of raw economic data
# and define as R time series object
# 'help(RawEconomicData)' provides additional information
D_raw <-  TVTP::RawEconomicData; help(RawEconomicData)
D     <-  stats::ts(D_raw, frequency = 12, start = c(1920, 1))
obs   <-  nrow(D)
print(colnames(D)); View(D)
help(RawEconomicData)


# calculate and transform raw time series
AAA_spreads         <-  D[,"rawAAAyld"] - D[,"raw20TsryYld"]
BBB_spreads         <-  D[,"rawBBByld"] - D[,"raw20TsryYld"]
treasury_bill       <-  D[,"raw3moBill"]
log_SPX             <-  log(D[,"rawSPX"])
log_production      <-  log(D[,"rawIndPro"])
yield_curve_slope   <-  D[,"raw20TsryYld"] - D[,"raw3moBill"]
ann_infl_rate       <-  D[,"rawAnnInflRate"]


#******************************************************************
#              Diagnostic Data Series Plots
#******************************************************************

# provide diagnostic plots
graphics.off()
windows()
plot(cbind(BBB_spreads,AAA_spreads),  plot.type = "single",
     lty=1:3,
     main = "BBB and AAA Credit Spreads",
     ylab = "% Credit Spread over risk-less alternative")
windows()
plot(cbind(BBB_spreads,AAA_spreads),
     lty=1:3,
     main = "BBB and AAA Credit Spreads",
     ylab = "% Credit Spread over risk-less alternative")


#******************************************************************
#          Data Preparations and Transformations
#******************************************************************

# extrapolate dependent and explanatory variable data sets
# > dependent variable: BBB_spreads or AAA_spreads
delta_spread        <-  diff(BBB_spreads)
# > explanatory variables governing regime transitions
trans_p_explanatory <-  cbind( rep(1,obs), # <- constant
                               lag(     ann_infl_rate,     -1))
# > explanatory variables governing the credit spread dynamics
process_explanatory <-  cbind( rep(1,obs), # <- constant
                               lag(diff(BBB_spreads),      -1),
                               lag(diff(treasury_bill),    -1),
                               lag(diff(log_SPX),          -1),
                               lag(diff(log_production),   -1),
                               lag(     yield_curve_slope, -1))


# Final data set for BBB spreads as dependent, and explanatory
# variables governing regime switch dynamics and the credit spread
# process. Note: The BBB and AAA dataset is also available
# precalculated accompanying this package.
# Type: 'help(BBB_dataset)' and 'help(AAA_dataset)'
BBB_dataset <-  stats::window(cbind(process_explanatory,
                                    trans_p_explanatory,
                                    delta_spread),
                              start = c(1920, 3),
                              end   = c(2018,12))
help(BBB_dataset)


#******************************************************************
#                Initializing Estimation Parameters
#******************************************************************

# Initializing starting values 'beta' which gowern the credit
# spread process:
#
#   Regime 0:                    Regime 1:
#    [1] constant                [7]  constant
#    [2] delta spread            [8]  delta spread
#    [3] delta t-bill            [9]  delta t-bill
#    [4] delta S&P               [10] delta S&P
#    [5] delta production        [11] delta production
#    [6] slope                   [12] slope
#
#               [13] volatility regime 0
#               [14] volatility regime 1
#
# Governing factors for the transition process between regimes.
# For each point in time, starting from regime:
#
#    Regime 0:                   Regime 1:
#    [15] constant               [17] constant
#    [16] inflation              [18] inflation
#
#                           t
#                state 0        state 1
#                --                   --
#     state 0   |   p_00        p_01    |
# t-1           |                       |
#     state 1   |   p_10        p_11    |
#                --                   --


# Initializing parameters:
paramInit  <-  c( #---phi_0------------------------------
                   0.077,    # 1   const_s0
                   0.150,    # 2   d_spread_s0
                   0.052,    # 3   d_bill_s0
                  -0.495,    # 4   d_log(spx)_s0
                   0.058,    # 5   d_log(prod)_s0
                  -0.004,    # 6   slope_s0
                  #---phi_1------------------------------
                   0.508,    # 7   const_s1
                   0.091,    # 8   d_spread_s1
                   0.017,    # 9   d_bill_s1
                  -0.879,    # 10  d_log(spx)_s1
                  -2.620,    # 11  d_log(prod)_s1
                   0.028,    # 12  slope_s1
                  #---simga------------------------------
                   0.075,    # 13  sigma_s0
                   0.349,    # 14  sigma_s1
                  #---beta_0-----------------------------
                   0.996,    # 15  theta_const_s0
                 -12.811,    # 16  theta_infl_s0
                  #---beta_1-----------------------------
                   0.976,    # 17  theta_const_s1
                   2.043)    # 18  theta_infl_s1


#******************************************************************
#   Initial log-Likelihood Given Data and Parmeter Estimates
#******************************************************************


# number of explanatory variables for the credit spread process
# including constant
n_explanatory  <- 6
# number of explanatory variables for the regime transition process
# including constant
n_transition   <- 2
# probability that initial state is in regime (0)
prob_initial_state_0  <- 0.99

# Calculate Log-Likelihood using fast C++ implementation:
res_C <- TVTP::LogLikeliHood_FULL(
                    param         =  paramInit,
                    datax         =  BBB_dataset,
                    n_explanatory =  n_explanatory,
                    n_transition  =  n_transition,
                    prob_initial_state_0  =  prob_initial_state_0,
                    printLogLike  =  1)

res_C$like_t_T_temp
res_C$LogLikelihoodVal
# LogLikelihoodVal (BBB): -313.10135
# LogLikelihoodVal (AAA):  148.4767



# Calculate Log-Likelihood using R implementation:
res_R <- TVTP::LogLikeliHood_FULL_R(
                    param         =  paramInit,
                    datax         =  BBB_dataset,
                    n_explanatory =  n_explanatory,
                    n_transition  =  n_transition,
                    prob_initial_state_0  =  prob_initial_state_0,
                    printLogLike  =  1)

res_R$filtered_joined_state_prob
res_R$smoothed_joint_state
res_R$smoothed_marginal_state_prob
res_R$cond_density
res_R$LogLikelihoodVal
# LogLikelihoodVal (BBB): -313.10135
# LogLikelihoodVal (AAA):  148.4767


#******************************************************************
#          Maximize Incomplete Data log-Likelihood
#******************************************************************

# load package 'optimx'. It provides the 'Nelder-Mead' (NM),
# 'Broyden-Fletcher-Goldfarb-Shanno' (BFGS), the
# 'Bounded Broyden-Fletcher-Goldfarb-Shanno' (L-BFGS-B) and the
# 'Conjugate Gradients' (CG) nonlinear optimization algorithm,
# among others
library(optimx)

# load package 'maxLik'. It provides the 'Brendt-Hall-Hall-Hausman'
# (BHHH) nonlinear optimization method among others
library(maxLik)



#---------------------------#
#  Derivative Free methods  #
#---------------------------#

# ==> NR [Newton-Rhapson]
# numerically maximize the log-likelyhood function using the EM-ML
# algorithm with the 'Newton-Rhapson' method
require(maxLik)
system.time(
res_maxLik_NR  <- maxLik::maxLik(
                      method  =  "NR",
                      logLik  =  LogLikeliHood, # C++ fun()
                      control =  list(iterlim = 1000),
                      start   =  paramInit,
                      datax         =  BBB_dataset,
                      n_explanatory =  n_explanatory,
                      n_transition  =  n_transition,
                      prob_initial_state_0 = prob_initial_state_0,
                      printLogLike  =  1
                    )
)
# max(Log-Like) BBB: -313.1014 user: 86.72

# ==> NM [Nelder-Mead]
# numerically maximize the log-likelyhood function using the EM-ML
# algorithm with the 'Nelder-Mead' method
require(stats)
system.time(
res_optim_NM   <- stats::optim(
                      method  =  "Nelder-Mead",
                      fn      =  LogLikeliHood, # C++ fun()
                      control =  list(  fnscale  =  -1,
                                        maxit    =  12000),
                      par     =  paramInit,
                      datax   =  BBB_dataset,
                      n_explanatory =  n_explanatory,
                      n_transition  =  n_transition,
                      prob_initial_state_0 = prob_initial_state_0,
                      printLogLike  =  1
                    )
)
# max(Log-Like) BBB: 762.385 user: 331.62


#---------------------------#
#  Gradient based methods   #
#---------------------------#

# ==> BFGS [Broyden-Fletcher-Goldfarb-Shanno]
# numerically maximize the log-likelyhood function using the EM-ML
# algorithm with the 'Broyden-Fletcher-Goldfarb-Shanno' method
require(stats) #Alternative: optimx::optimx which is double as slow
system.time(
res_optim_BFGS   <- stats::optim(
                      method  =  "BFGS",
                      fn      =  LogLikeliHood, # C++ fun()
                      control =  list(  fnscale  =  -1,
                                        maxit    =  300),
                      par     =  paramInit,
                      datax   =  BBB_dataset,
                      n_explanatory =  n_explanatory,
                      n_transition  =  n_transition,
                      prob_initial_state_0 = prob_initial_state_0,
                      printLogLike  =  1
                    )
  )
# max(Log-Like) BBB: 823.5807 user: 40.46


# ==> BHHH [Brendt-Hall-Hall-Hausman]
# numerically maximize the log-likelyhood function using the EM-ML
# algorithm with the 'Brendt-Hall-Hall-Hausman' method
require(maxLik)
system.time(
res_maxLik_BHHH  <- maxLik::maxLik(
                      method  =  "BHHH",
                      logLik  =  LogLikeliHood_byObs, # C++ fun()
                      control =  list(iterlim = 1000),
                      start   =  paramInit,
                      datax         =  BBB_dataset,
                      n_explanatory =  n_explanatory,
                      n_transition  =  n_transition,
                      prob_initial_state_0 = prob_initial_state_0,
                      printLogLike  =  1
                    )
)
# max(Log-Like) BBB: 823.5844 user: 391.67


# ==> CG [Conjugated Gradients]
# numerically maximize the log-likelyhood function using the EM-ML
# algorithm with the 'Conjugated Gradients' method
require(stats)  #Alternative optimx::optimx same speed as stats::optim
system.time(
res_optim_CG   <- stats::optim(
                      method  =  "CG",
                      fn      =  LogLikeliHood, # C++ fun()
                      control =  list(  fnscale  =  -1,
                                        maxit    =  12000),
                      par     =  paramInit,
                      datax   =  BBB_dataset,
                      n_explanatory =  n_explanatory,
                      n_transition  =  n_transition,
                      prob_initial_state_0 = prob_initial_state_0,
                      printLogLike  =  1
                    )
)
#max(Log-Like) BBB: 823.2737 user: 2755.58


# ==> NLM [Non-Linear Minimization]
# numerically maximize the log-likelyhood function using the EM-ML
# algorithm with the 'Non-Linear Minimization' method
require(stats)
system.time(
res_nlm          <- stats::nlm(
                      f  =  LogLikeliHood_min,
                      p  =  paramInit,
                      iterlim  = 1000,
                      datax         =  BBB_dataset,
                      n_explanatory =  n_explanatory,
                      n_transition  =  n_transition,
                      prob_initial_state_0 = prob_initial_state_0,
                      printLogLike  =  1
                    )
)
# max(Log-Like) BBB: 823.585 user: 109.29



