

#************************************************************************************************************
#********                                        DS-Equity R&D                                       ********
#********      Spread Dynamics proxied with Regime switch models with Time Varying Transition Prob.  ********
#********                                   Markov TVTP Estimation                                   ********
#************************************************************************************************************



#***************************************************************************************************
#                        Reference code to be used similar to the C++ implementation
#***************************************************************************************************


# options(max.print=110000) # max print of individual values, not rows!


# Sample call
#-------------------
# n_explanatory = 6
# n_transition  = 2
# prob_initial_state_0 = 0.99

# res <- LogLikeliHood_FULL_R(
#   param         =  paramInit,
#   datax         =  BBB_dataset,
#   n_explanatory =  n_explanatory,
#   n_transition  =  n_transition,
#   prob_initial_state_0  =  prob_initial_state_0,
#   printLogLike  =  1)


#' @title LogLikeliHood_FULL_R
#'
#' @description
#' \code{LogLikeliHood_FULL_R} returns a numeric scalar with the log-likelihood of the Expectation-Maximization
#' Maximum-Likelihood (EM-ML) algorithm given observations for the explanatory variables governing direct
#' effects on the the dependent variable and observations for explanatory variables that govern the regime
#' transition process. This implemenation is several thousand times slower compared to the fast C++ implementation.
#'     Author: Tobias Gummersbach
#' @param x,y two numeric values to be added up
#' @return The addition of the inputs
#' @examples
#' ## add up two numbers
#' GoTobiGo(1,3)
#'
#' ## add up two values
#' x <- 1
#' y <- 3
#' GoTobiGo(x,y)
#' @section Author:
#' Tobias Gummersbach
#' @export



LogLikeliHood_FULL_R <- function(
  param, datax, n_explanatory, n_transition,
  prob_initial_state_0, printLogLike=1){


# Extract variables
#-------------------
obs                 <- NROW(datax)
process_explanatory <- datax[,1:n_explanatory]
trans_p_explanatory <- datax[,(n_explanatory+1):(n_explanatory+n_transition)]
delta_spread        <- datax[,NCOL(datax)]

beta_               <- param[1:(2*n_explanatory+2)]
beta1               <- param[1:n_explanatory]
beta2               <- param[(n_explanatory+1):(n_explanatory*2)]
vol1                <- param[(2*n_explanatory+1)]
vol2                <- param[(2*n_explanatory+2)]
theta_              <- param[(length(param)-n_transition*2+1):length(param)]


# Step 1. Initialize
#-------------------
# Obtain an initial set of transition probabilities, depending on explanatory variables --> [T x 4]
help1       <-  exp(trans_p_explanatory %*% (cbind(theta_[1:2],theta_[3:4])))
help2       <-  help1 / (1 + help1) #contains p_00 and p_11
trans_prob  <-  cbind(help2[,1], 1-help2[,2], 1-help2[,1], help2[,2]) # contains p_00 ~ p_10 ~ p_01 ~ p_11. Read p_01: transition from state 0 to state 1

# Obtain an initial set of conditional densities, depending upon one of the two regimes: Density of observing y_t given regime is in state 0|1  --> [T x 2]
cond_density  <-  cbind( 1/(sqrt(2*pi)*vol1) * exp(-0.5*((delta_spread - process_explanatory %*% beta1)/vol1)^2),
                         1/(sqrt(2*pi)*vol2) * exp(-0.5*((delta_spread - process_explanatory %*% beta2)/vol2)^2)  )
colnames(cond_density) <- c("cond_desity_s0", "cond_desity_s1")

#cond_density[cond_density == 0]  <-  .Machine$double.eps  # where cond_density == 0, replace by a very small value to prevent ln(.) from becoming -INF which stops cml routine
# Type: ?.Machine  for info: smallest value x such that 1+x != 0: here: 2.220446e-16
# Note, this is not equivalent to the check 2.220446e-20 == 0 which will still say FALSE!
cond_density[cond_density < 2^(-1074)]  <-  2^(-1074)      # 2^(-1074) == 5e-324;   5e-324 < .Machine$double.eps;
# According to ?.Machine, On a typical R platform the smallest positive double is about 5e-324


# Step 2. Calculation of filtered joint state probabilities (Recursive algorithm, for each t until T)
#----------------------------------------------------------------------------------------------------

joint_density               <-  matrix(0,obs,4)  # mat_d
filtered_joined_state_prob  <-  matrix(0,obs,4)  # mat_e

#prob_initial_state_0       <-  0.5     # a-priori assumption for state 0 in t=1
prob_initial_state_1        <-  1-prob_initial_state_0     # a-priori assumption for state 1 in t=1

# 2.a
# Calculate [filtered_joint_state_prob] for t = 2
joint_density[2,1]  <-  cond_density[2,1] * trans_prob[2,1] * prob_initial_state_0;   # 0 -> 0
joint_density[2,2]  <-  cond_density[2,1] * trans_prob[2,2] * prob_initial_state_1;   # 1 -> 0
joint_density[2,3]  <-  cond_density[2,2] * trans_prob[2,3] * prob_initial_state_0;   # 0 -> 1
joint_density[2,4]  <-  cond_density[2,2] * trans_prob[2,4] * prob_initial_state_1;   # 1 -> 1

# Filtered state probabilities for t = 2 (see step 2.c)
filtered_joined_state_prob[2,] = joint_density[2,] / sum(joint_density[2,])

# Repeat steps 2a to 2c for t = 3 to T
get_right_a <- c(1,1,2,2)
get_right_b <- matrix(c(1, 2, 3, 4, 1, 2, 3, 4), 4, 2, byrow = TRUE)
for (i in 3:obs){
  #step 2.a
  for (j in 1:4){
    joint_density[i,j]  <-  cond_density[i,get_right_a[j]] * trans_prob[i,j] *
      (filtered_joined_state_prob[i-1, get_right_b[j,1]] +  filtered_joined_state_prob[i-1, get_right_b[j,2]])
  }
  #step 2.b and 2.c
  filtered_joined_state_prob[i,] <- joint_density[i,] / sum(joint_density[i,])
}



#  Step 3. Calculation of smoothed joint state probabilities
#-----------------------------------------------------------
#  recursive iteration: steps 3a-->3c, for each t = 2 until T
#  recursive iteration: steps 3a-->3b, for each set of states at t: [s_t=0,s_(t-1)=0], [s_t=0 ,s_(t-1)=1], [s_t=1 ,s_(t-1)=0], [s_t=1,s_(t-1)=1]
#  recursive iteration: for each tau = t+2(=4),....,T

smoothed_joint_state  <-  matrix(1,obs,4)  # mat_g
get_right_d           <-  c(1,1,2,2)
get_right_e           <-  matrix(c(1, 2, 3, 4, 1, 2, 3, 4), 4, 2, byrow = TRUE)
time_t=2
state_set = 1
#if(is.na(Sys.getenv("BINPREF", unset=NA))) {Sys.setenv("BINPREF"="C:/Rtools/mingw_$(WIN)/bin/")}


for (time_t in 2:obs){ #loops: 3a-->3c, start at t = 2 to T, fills matrix: smoothed_joint_state

  # 3.a Initialize smoothed_joint_state_probability at each t = 2,...,T for four different state_sets

  for (state_set in 1:4){ # loops: 3a-->3b, [s_t=0,s_(t-1)=0], [s_t=0 ,s_(t-1)=1], [s_t=1 ,s_(t-1)=0], [s_t=1,s_(t-1)=1]

    if (time_t != obs){
      ####################
      # tau = t+1:       #  We calculate tau = t+1 only as an input to being able to calculating tau = t+2.
      ####################
      get_right_c = c(0,0,1,1) # helps to choose the right index from trans_prob
      joint_prob_initial   <-  cond_density[time_t+1,] * trans_prob[time_t+1, c(1+get_right_c[state_set], 3+get_right_c[state_set])] * filtered_joined_state_prob[time_t,state_set]
      joint_prob_initial   <-  joint_prob_initial / sum(joint_density[time_t+1,])  # P(0->0->0) & P(0->0->1)
    }


    if (time_t == obs){
      smoothed_joint_state[time_t,state_set] <-  filtered_joined_state_prob[time_t,state_set]
      # smoothed_joint_state probabilities equal the filtered_joined_state probabilities for the last observation as
      # filtered_joint_state probabilities equal the probabilities for the states at t and t-1 conditioning
      # upon all observations up to T.
    } else if (time_t == obs-1){
      smoothed_joint_state[time_t,state_set] <-  sum(joint_prob_initial)
    } else {
      ####################
      # tau = t+2        #
      ####################
      joint_prob              <-  matrix(1,obs,4)
      joint_prob[time_t+2,1]  <-  cond_density[time_t+2,1] * trans_prob[time_t+2,1] * joint_prob_initial[1]
      joint_prob[time_t+2,2]  <-  cond_density[time_t+2,1] * trans_prob[time_t+2,2] * joint_prob_initial[2]
      joint_prob[time_t+2,3]  <-  cond_density[time_t+2,2] * trans_prob[time_t+2,3] * joint_prob_initial[1]
      joint_prob[time_t+2,4]  <-  cond_density[time_t+2,2] * trans_prob[time_t+2,4] * joint_prob_initial[2]
      joint_prob[time_t+2, ]  <-  joint_prob[time_t+2,] / sum(joint_density[time_t+2,])


      ####################
      # tau = t+3,...,T  #
      ####################
      if (time_t < (obs-2)){
        for (i in (time_t+3):obs){
          for (j in 1:4){
            joint_prob[i,j] <-  cond_density[i,get_right_d[j]] * trans_prob[i,j] * (joint_prob[i-1,get_right_e[j,1]] + joint_prob[i-1,get_right_e[j,2]])
          }
          joint_prob[i,] <-  joint_prob[i,] / sum(joint_density[i,])
        }
      }




      # Step: 3b
      smoothed_joint_state[time_t,state_set]  <-  sum(joint_prob[obs,])
    } # end if

  } # next state set

} # next time t



#  Step 4. Calculation of smoothed marginal state probabilities
#--------------------------------------------------------------

smoothed_marginal_state_prob      <-  cbind(smoothed_joint_state[,1] + smoothed_joint_state[,2], smoothed_joint_state[,3] + smoothed_joint_state[,4])
smoothed_marginal_state_prob[1,]  <-  cbind(prob_initial_state_0, prob_initial_state_1)


#  Plot Initial results
#--------------------------------------------------------------

windows()
plot(ts(data.frame(  delta_spread      =  cumsum(delta_spread),
                     Prob_BeingIn0     =  smoothed_marginal_state_prob[,1],
                     Prop_Trans_0_To_1 =  trans_prob[,3],
                     Prop_Trans_1_0    =  trans_prob[,2]), frequency = 12, start = c(1921, 1))   , plot.type="multiple", lty=1:3,
     main = "[AAA delta Spreads], [Prob of being in state 0] and [Probability to transition 0 to 1 and 1 to 0]")

cat(paste("Count of estimated observations in state 0:", sum(smoothed_marginal_state_prob[,1] > 0.5) ,"(out of", obs, "observations total).\n"))


#  Estimate complete data Log-Likelihood
#--------------------------------------------------------------

# in previous versions, also in GAUSS, cond_density[2,1] and cond_density[2,2] were used, incorrectly
like_t_1         <-  prob_initial_state_0 * (log(cond_density[1,1]) + log(prob_initial_state_0)) +
  prob_initial_state_1 * (log(cond_density[1,2]) + log(prob_initial_state_1))
like_t_T_temp    <-  smoothed_marginal_state_prob[,1] * log(cond_density[,1]) +
  smoothed_marginal_state_prob[,2] * log(cond_density[,2]) +
  rowSums(smoothed_joint_state * log(trans_prob))

like_t_T_temp[1] <-  like_t_1

return(list(  LogLikelihoodVal             =  sum(like_t_T_temp),
              filtered_joined_state_prob   =  filtered_joined_state_prob,
              smoothed_joint_state         =  smoothed_joint_state,
              smoothed_marginal_state_prob =  smoothed_marginal_state_prob,
              cond_density                 =  cond_density))}


