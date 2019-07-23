
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXX                                                                                                                   XXXX
//XXXX                                              Dissertation Effort                                                  XXXX
//XXXX                                                 Working Title:                                                    XXXX
//XXXX                                   The Application of Regime Switch Models with                                    XXXX
//XXXX                              Time Varying Transition Probabilities to Solving the                                 XXXX
//XXXX                                            Credit Spread Puzzle                                                   XXXX
//XXXX                                                                                                                   XXXX
//XXXX                         Required C++ code accompanying R-package for package compilation                          XXXX
//XXXX                                                                                                                   XXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



#include <Rcpp.h>    // This is important, always has to be there
#include <math.h>
//--// #include <omp.h>           // both lines need to be enabled for parallelization. Search in the code
//--// [[Rcpp::plugins(openmp)]]  // down below for 'pragma'


using namespace Rcpp;

// General Notes:
// NEVER HAVE A FUNCTION DECLARATION FUN() WHERE THE IMMEDIATE NEXT LINE IS COMMENTED OUT. THIS WILL THROW ERRORS



//---------------------------------------------------------------------------------------------------------------------------
// LogLikeliHood_FULL
// Returns the log-likelihood of observing the supplied 'datax' given 'param'. Additional values are provided for diagnostics
// such as each individual observation's log-likelihood contribution. If 'printLogLike'=1, the resulting log-likelihood is
// printed to the R console. Date of creation: 04/08/2017
// Date of change: 07/13/2019: Added new functionality: Allowed for up to 10 process_explanatory variables, formerly only 7.
// Date of change: 07/21/2019: Added roxygen comments below.
//---------------------------------------------------------------------------------------------------------------------------

// roxygen comments:

//' LogLikeliHood_FULL
//'
//' \code{LogLikeliHood_FULL} is a fast \code{C++} implementation that returns a numeric scalar \code{c(0)} with the log-likelihood of the
//' Expectation-Maximization Maximum-Likelihood (EM-ML) algorithm given observations for the explanatory variables governing direct
//' effects on the dependent variable, and observations for explanatory variables that govern the regime
//' transition process. It also returns various diagnostic statistics such as each individual observation's log-likelihood
//' contribution. If parameter 'printLogLike' = 1, the sum of the individual observation's log-likelihood contributions
//' is printed to the console. This fast \code{C++} implementation is several thousand times fast to the \R analog.
//' @export
// [[Rcpp::export]]
Rcpp::List LogLikeliHood_FULL(
  Rcpp::NumericVector param,  // param
  Rcpp::NumericMatrix datax,  // data
  int n_explanatory,          // number of explanatory variables governing 'process_explanatory' including constant [1..7]
  int n_transition,           // number of explanatory variables governing 'trans_p_explanatory' including constant [1,2]
  double prob_initial_state_0,// prob_initial_state_0: Probability that the first observation is in state '0' [0.00001, 0.99999] We take the log(prob_initial_state_0) in the final stage which must not become -Inf
  int printLogLike            // 0: Do not print output to console
){

  //======================================
  // Variable declaration, data extraction
  //======================================

  // Input variables 'param' and 'datax' always follow the same format:
  int obs = datax.nrow();
  Rcpp::NumericVector beta_(n_explanatory * 2);  // parameter set for process_explanatory (excl. volatilities)
  Rcpp::NumericVector sigma_(2);                 // parameter set for volatilities for process_explanatory
  Rcpp::NumericVector theta_(n_transition * 2);  // parameter set for trans_p_explanatory

  // int n_param = param.size();

  //Extract explanatory data from 'datax'
  // ... process_explanatory:
  Rcpp::NumericMatrix process_explanatory = datax( _, Range(0, n_explanatory-1));
  // ... trans_p_explanatory:
  Rcpp::NumericMatrix trans_p_explanatory = datax( _, Range(n_explanatory, n_explanatory + n_transition - 1));
  // ... delta_spread:
  Rcpp::NumericVector delta_spread = datax( _, n_explanatory + n_transition);


  // Extract parameters theta_ and beta_:
  switch(n_explanatory)
    {
    case 10:                          // regime 0              regime 1
    beta_  = param[IntegerVector::create(0,1,2,3,4,5,6,7,8,9,  10,11,12,13,14,15,16,17,18,19)];
    sigma_ = param[IntegerVector::create(20,21)];
    if(n_transition == 2){
      theta_ = param[IntegerVector::create(22,23,  24,25)];
    } else {
      theta_ = param[IntegerVector::create(22,23)];
    }
    break;

    case 9:
    beta_  = param[IntegerVector::create(0,1,2,3,4,5,6,7,8,  9,10,11,12,13,14,15,16,17)];
    sigma_ = param[IntegerVector::create(18,19)];
    if(n_transition == 2){
      theta_ = param[IntegerVector::create(20,21,  22,23)];
    } else {
      theta_ = param[IntegerVector::create(20,21)];
    }
    break;

    case 8:
      beta_  = param[IntegerVector::create(0,1,2,3,4,5,6,7  ,8,9,10,11,12,13,14,15)];
      sigma_ = param[IntegerVector::create(16,17)];
    if(n_transition == 2){
      theta_ = param[IntegerVector::create(18,19,  20,21)];
    } else {
      theta_ = param[IntegerVector::create(18,19)];
    }
    break;

    case 7:
      beta_  = param[IntegerVector::create(0,1,2,3,4,5,6,  7,8,9,10,11,12,13)];
      sigma_ = param[IntegerVector::create(14,15)];
      if(n_transition == 2){
      theta_ = param[IntegerVector::create(16,17,  18,19)];
      } else {
      theta_ = param[IntegerVector::create(16,17)];
      }
      break;

    case 6:
      beta_  = param[IntegerVector::create(0,1,2,3,4,5,  6,7,8,9,10,11)];
      sigma_ = param[IntegerVector::create(12,13)];
      if(n_transition == 2){
      theta_ = param[IntegerVector::create(14,15,  16,17)];
      } else {
      theta_ = param[IntegerVector::create(14,15)];
      }
      break;

    case 5:
      beta_  = param[IntegerVector::create(0,1,2,3,4,  5,6,7,8,9)];
      sigma_ = param[IntegerVector::create(10,11)];
      if(n_transition == 2){
      theta_ = param[IntegerVector::create(12,13,  14,15)];
      } else {
      theta_ = param[IntegerVector::create(12,13)];
      }
      break;

    case 4:
      beta_  = param[IntegerVector::create(0,1,2,3,  4,5,6,7)];
      sigma_ = param[IntegerVector::create(8,  9)];
      if(n_transition == 2){
      theta_ = param[IntegerVector::create(10,11,  12,13)];
      } else {
      theta_ = param[IntegerVector::create(10,11)];
      }
      break;

    case 3:
      beta_  = param[IntegerVector::create(0,1,2,  3,4,5)];
      sigma_ = param[IntegerVector::create(6,7)];
      if(n_transition == 2){
      theta_ = param[IntegerVector::create(8,9,  10,11)];
      } else {
      theta_ = param[IntegerVector::create(8,9)];
      }
      break;

    case 2:
      beta_  = param[IntegerVector::create(0,1,  2,3)];
      sigma_ = param[IntegerVector::create(4,5)];
      if(n_transition == 2){
      theta_ = param[IntegerVector::create(6,7,  8, 9)];
      } else {
      theta_ = param[IntegerVector::create(6,7)];
      }
      break;

    case 1:
      beta_  = param[IntegerVector::create(0,  1)];
      sigma_ = param[IntegerVector::create(2,3)];
      if(n_transition == 2){
      theta_ = param[IntegerVector::create(4,5,  6, 7)];
      } else {
      theta_ = param[IntegerVector::create(4,5)];
      }
      break;
    }

  // Alternative construction:
  //  IntegerVector idx = IntegerVector::create(0,1,2,3,4,5,6,7);
  //  beta_ = param[idx];


  //=========
  // Step: 1
  //=========

  // Obtain an initial set of transition probabilities, depending on explanatory variables --> [T x 4] matrix 'mat_b'
  // not used: Rcpp::NumericMatrix trans_prob = GET_trans_prob(trans_p_explanatory, theta_, obs);

  Rcpp::NumericMatrix help1(obs,2);
  Rcpp::NumericMatrix help2(obs,2);  // contains p_00 and p_11
  Rcpp::NumericMatrix mat_b(obs,4);  // trans_prob: Contains p_00 ~ p_10 ~ p_01 ~ p_11. Read p_01: transition from state 0 to state 1

  for(int i = 0; i < obs; i++){

    switch(n_transition)
      {
      case 2:
        help1(i,0) = exp(trans_p_explanatory(i,0) * theta_(0) + trans_p_explanatory(i,1) * theta_(1));
        help1(i,1) = exp(trans_p_explanatory(i,0) * theta_(2) + trans_p_explanatory(i,1) * theta_(3));
        break;

      case 1:
        help1(i,0) = exp(trans_p_explanatory(i,0) * theta_(0));
        help1(i,1) = exp(trans_p_explanatory(i,0) * theta_(1));
        break;
      }

    help2(i,0) = help1(i,0) / (1+help1(i,0));
    help2(i,1) = help1(i,1) / (1+help1(i,1));

    mat_b(i,0) = help2(i,0);     // p_00
    mat_b(i,1) = 1 - help2(i,1); // p_10
    mat_b(i,2) = 1 - help2(i,0); // p_01
    mat_b(i,3) = help2(i,1);     // p_11
  }

  // Obtain an initial set of conditional densities, depending upon one of the two regimes: Density of observing y_t given regime is in state 0|1  --> [T x 2] matrix 'mat_a'
  // not used: Rcpp::NumericMatrix cond_density = GET_cond_density(process_explanatory, delta_spread, beta_, obs);

  Rcpp::NumericMatrix mat_a(obs,2);  // cond_density: 'mat_a'
  double ExpMean1 = 0;
  double ExpMean2 = 0;

  for(int i = 0; i < obs; i++){
    switch(n_explanatory)
      {
      case 10:
        ExpMean1 = process_explanatory(i,0) * beta_(0)  +
                   process_explanatory(i,1) * beta_(1)  +
                   process_explanatory(i,2) * beta_(2)  +
                   process_explanatory(i,3) * beta_(3)  +
                   process_explanatory(i,4) * beta_(4)  +
                   process_explanatory(i,5) * beta_(5)  +
                   process_explanatory(i,6) * beta_(6)  +
                   process_explanatory(i,7) * beta_(7)  +
                   process_explanatory(i,8) * beta_(8)  +
                   process_explanatory(i,9) * beta_(9);
        ExpMean2 = process_explanatory(i,0) * beta_(10) +
                   process_explanatory(i,1) * beta_(11) +
                   process_explanatory(i,2) * beta_(12) +
                   process_explanatory(i,3) * beta_(13) +
                   process_explanatory(i,4) * beta_(14) +
                   process_explanatory(i,5) * beta_(15) +
                   process_explanatory(i,6) * beta_(16) +
                   process_explanatory(i,7) * beta_(17) +
                   process_explanatory(i,8) * beta_(18) +
                   process_explanatory(i,9) * beta_(19);
        break;

      case 9:
        ExpMean1 = process_explanatory(i,0) * beta_(0)  +
                   process_explanatory(i,1) * beta_(1)  +
                   process_explanatory(i,2) * beta_(2)  +
                   process_explanatory(i,3) * beta_(3)  +
                   process_explanatory(i,4) * beta_(4)  +
                   process_explanatory(i,5) * beta_(5)  +
                   process_explanatory(i,6) * beta_(6)  +
                   process_explanatory(i,7) * beta_(7)  +
                   process_explanatory(i,8) * beta_(8);
        ExpMean2 = process_explanatory(i,0) * beta_(9)  +
                   process_explanatory(i,1) * beta_(10) +
                   process_explanatory(i,2) * beta_(11) +
                   process_explanatory(i,3) * beta_(12) +
                   process_explanatory(i,4) * beta_(13) +
                   process_explanatory(i,5) * beta_(14) +
                   process_explanatory(i,6) * beta_(15) +
                   process_explanatory(i,7) * beta_(16) +
                   process_explanatory(i,8) * beta_(17);
        break;

      case 8:
        ExpMean1 = process_explanatory(i,0) * beta_(0)  +
                   process_explanatory(i,1) * beta_(1)  +
                   process_explanatory(i,2) * beta_(2)  +
                   process_explanatory(i,3) * beta_(3)  +
                   process_explanatory(i,4) * beta_(4)  +
                   process_explanatory(i,5) * beta_(5)  +
                   process_explanatory(i,6) * beta_(6)  +
                   process_explanatory(i,7) * beta_(7);
        ExpMean2 = process_explanatory(i,0) * beta_(8)  +
                   process_explanatory(i,1) * beta_(9)  +
                   process_explanatory(i,2) * beta_(10) +
                   process_explanatory(i,3) * beta_(11) +
                   process_explanatory(i,4) * beta_(12) +
                   process_explanatory(i,5) * beta_(13) +
                   process_explanatory(i,6) * beta_(14) +
                   process_explanatory(i,7) * beta_(15);
        break;

      case 7:
        ExpMean1 = process_explanatory(i,0) * beta_(0)  +
                   process_explanatory(i,1) * beta_(1)  +
                   process_explanatory(i,2) * beta_(2)  +
                   process_explanatory(i,3) * beta_(3)  +
                   process_explanatory(i,4) * beta_(4)  +
                   process_explanatory(i,5) * beta_(5)  +
                   process_explanatory(i,6) * beta_(6);
        ExpMean2 = process_explanatory(i,0) * beta_(7)  +
                   process_explanatory(i,1) * beta_(8)  +
                   process_explanatory(i,2) * beta_(9)  +
                   process_explanatory(i,3) * beta_(10) +
                   process_explanatory(i,4) * beta_(11) +
                   process_explanatory(i,5) * beta_(12) +
                   process_explanatory(i,6) * beta_(13);
        break;

      case 6:
        ExpMean1 = process_explanatory(i,0) * beta_(0)  +
                   process_explanatory(i,1) * beta_(1)  +
                   process_explanatory(i,2) * beta_(2)  +
                   process_explanatory(i,3) * beta_(3)  +
                   process_explanatory(i,4) * beta_(4)  +
                   process_explanatory(i,5) * beta_(5);
        ExpMean2 = process_explanatory(i,0) * beta_(6)  +
                   process_explanatory(i,1) * beta_(7)  +
                   process_explanatory(i,2) * beta_(8)  +
                   process_explanatory(i,3) * beta_(9)  +
                   process_explanatory(i,4) * beta_(10) +
                   process_explanatory(i,5) * beta_(11);
        break;

      case 5:
        ExpMean1 = process_explanatory(i,0) * beta_(0) +
                   process_explanatory(i,1) * beta_(1) +
                   process_explanatory(i,2) * beta_(2) +
                   process_explanatory(i,3) * beta_(3) +
                   process_explanatory(i,4) * beta_(4);
        ExpMean2 = process_explanatory(i,0) * beta_(5) +
                   process_explanatory(i,1) * beta_(6) +
                   process_explanatory(i,2) * beta_(7) +
                   process_explanatory(i,3) * beta_(8) +
                   process_explanatory(i,4) * beta_(9);
        break;

      case 4:
        ExpMean1 = process_explanatory(i,0) * beta_(0) +
                   process_explanatory(i,1) * beta_(1) +
                   process_explanatory(i,2) * beta_(2) +
                   process_explanatory(i,3) * beta_(3);
        ExpMean2 = process_explanatory(i,0) * beta_(4) +
                   process_explanatory(i,1) * beta_(5) +
                   process_explanatory(i,2) * beta_(6) +
                   process_explanatory(i,3) * beta_(7);
        break;

      case 3:
        ExpMean1 = process_explanatory(i,0) * beta_(0) +
                   process_explanatory(i,1) * beta_(1) +
                   process_explanatory(i,2) * beta_(2);
        ExpMean2 = process_explanatory(i,0) * beta_(3) +
                   process_explanatory(i,1) * beta_(4) +
                   process_explanatory(i,2) * beta_(5);
        break;

      case 2:
        ExpMean1 = process_explanatory(i,0) * beta_(0) +
                   process_explanatory(i,1) * beta_(1);
        ExpMean2 = process_explanatory(i,0) * beta_(2) +
                   process_explanatory(i,1) * beta_(3);
        break;

      case 1:
        ExpMean1 = process_explanatory(i,0) * beta_(0);
        ExpMean2 = process_explanatory(i,0) * beta_(1);
        break;
      }

    mat_a(i,0) = 1 / (sqrt(2*M_PI)*sigma_[0]) * exp(-0.5 * pow(((delta_spread(i) - ExpMean1)/sigma_(0)),2));  //conditional density observation i is in state '0'
    mat_a(i,1) = 1 / (sqrt(2*M_PI)*sigma_[1]) * exp(-0.5 * pow(((delta_spread(i) - ExpMean2)/sigma_(1)),2));  //conditional density observation i is in state '1'
  }


  // where cond_density == 0, replace by a very small negative value to prevent log(.) from becoming -INF which might stop
  // optimization algorithms
  for(int i = 0; i < obs; i++){
    for(int j = 0; j < 2; j++){
      if(mat_a(i,j) == 0){
        mat_a(i,j) = pow(2, -1074); // This is the smallest floating number that can still be interpreted in R. For consitency reasons with R, we keep it this way
      }
    }
  }

  // Step 2:
  // Calculation of filtered joint state probabilities (Recursive algorithm, for each t until T)
  // not used:  Rcpp::NumericMatrix R_cond_density = GET_filtered_joined_state_prob(R_cond_density, R_trans_prob, prob_initial_state_0, prob_initial_state_1, obs);

  //===============================================
  // Step: 2a, 2b and 2c for t=2 and then t=3 to T
  //===============================================

  double RowSum_i = 0;
  double prob_initial_state_1 = 1 - prob_initial_state_0; // prob_initial_state_1
  Rcpp::NumericMatrix mat_d(obs,4); // joint_density:              'mat_d'
  Rcpp::NumericMatrix mat_e(obs,4); // filtered_joined_state_prob: 'mat_e'

  // Initialize mat_e
  // mat_e = mat_e + 0.25; // <- This is probably not used at all. delete later

  //********
  //  t = 2
  //********

  // 2a: Calculate [joint_density] for t = 2
  mat_d(1,0)  =  mat_a(1,0) * mat_b(1,0) * prob_initial_state_0;   // switch 0 -> 0
  mat_d(1,1)  =  mat_a(1,0) * mat_b(1,1) * prob_initial_state_1;   // switch 1 -> 0
  mat_d(1,2)  =  mat_a(1,1) * mat_b(1,2) * prob_initial_state_0;   // switch 0 -> 1
  mat_d(1,3)  =  mat_a(1,1) * mat_b(1,3) * prob_initial_state_1;   // switch 1 -> 1

  // 2b and 2c: Filtered state probabilities for t = 2
  RowSum_i = mat_d(1,0) + mat_d(1,1) + mat_d(1,2) + mat_d(1,3);
  mat_e(1,0) = mat_d(1,0) / RowSum_i;
  mat_e(1,1) = mat_d(1,1) / RowSum_i;
  mat_e(1,2) = mat_d(1,2) / RowSum_i;
  mat_e(1,3) = mat_d(1,3) / RowSum_i;

  //*************
  //  t = 3 to T
  //*************

  // Repeat steps 2a, 2b and 2c for t = 3 to T
  for(int i = 2; i < obs; i++){
    mat_d(i,0)  =  mat_a(i,0) * mat_b(i,0) * (mat_e(i-1, 0) + mat_e(i-1, 1));
    mat_d(i,1)  =  mat_a(i,0) * mat_b(i,1) * (mat_e(i-1, 2) + mat_e(i-1, 3));
    mat_d(i,2)  =  mat_a(i,1) * mat_b(i,2) * (mat_e(i-1, 0) + mat_e(i-1, 1));
    mat_d(i,3)  =  mat_a(i,1) * mat_b(i,3) * (mat_e(i-1, 2) + mat_e(i-1, 3));

    RowSum_i = mat_d(i,0) + mat_d(i,1) + mat_d(i,2) + mat_d(i,3);

    mat_e(i,0) = mat_d(i,0) / RowSum_i;
    mat_e(i,1) = mat_d(i,1) / RowSum_i;
    mat_e(i,2) = mat_d(i,2) / RowSum_i;
    mat_e(i,3) = mat_d(i,3) / RowSum_i;
  }


  // Step 3:
  // Calculation of smoothed joint state probabilities
  // recursive iteration: steps 3a-->3c, for each t = 2 until T
  // recursive iteration: steps 3a-->3b, for each set of states at t: [s_t=0,s_(t-1)=0], [s_t=0 ,s_(t-1)=1], [s_t=1 ,s_(t-1)=0], [s_t=1,s_(t-1)=1]
  // recursive iteration: for each tau = t+2(=4),....,T

  //=============
  // Step: 3a & b
  //=============

  Rcpp::NumericMatrix mat_g(obs,4); // smoothed_joint_state: 'mat_g'
  Rcpp::NumericMatrix mat_f(1,2);   // joint_prob_initial:   'mat_f'
  Rcpp::NumericMatrix mat_c(obs,4); // joint_prob:           'mat_c'

  //int state_set = 1;  // Stay true to R code
  //int time_t = 1;     // substracted one compared to R code

  //--// http://bisqwit.iki.fi/story/howto/openmp/#ExampleInitializingATableInParallelMultipleThreads:
  //--// private: Each thread has it's own copy of that variable
  //--// shared: Each thread accesses the same variable
  //--// #pragma omp parallel for  shared(mat_a, mat_b, mat_e) private(i, mat_d,  RowSum_i)
  for(int time_t = 1; time_t < obs; time_t++){           //loops: 3a-->3c, start at t = 2 to T, fills matrix: smoothed_joint_state 'mat_g'
    for(int state_set = 1; state_set < 5; state_set++){  //loops: 3a-->3b, [s_t=0,s_(t-1)=0], [s_t=0 ,s_(t-1)=1], [s_t=1 ,s_(t-1)=0], [s_t=1,s_(t-1)=1]

      // Get a single smoothed joint state probability: The content of this loop can alternatively be accessed by this function, which would be much slower, though.
      // not used:  mat_g(time_t,state_set-1) = GET_SINGLE_smoothed_joint_state(time_t, state_set, mat_a, mat_b, mat_d, mat_e);

      //====================
      //  tau = t+1
      //====================

      if (time_t < (obs-1)){
        RowSum_i = mat_d(time_t+1,0) + mat_d(time_t+1,1) + mat_d(time_t+1,2) + mat_d(time_t+1,3);

        switch(state_set)
        {
        case 1:
          mat_f(0,0) = mat_a(time_t+1,0) * mat_b(time_t+1, 0) * mat_e(time_t,state_set-1) / RowSum_i;  // joint probability 0->0 -> 0
          mat_f(0,1) = mat_a(time_t+1,1) * mat_b(time_t+1, 2) * mat_e(time_t,state_set-1) / RowSum_i;  // joint probability 0->0 -> 1
          break;

        case 2: //identical to case 1:
          mat_f(0,0) = mat_a(time_t+1,0) * mat_b(time_t+1, 0) * mat_e(time_t,state_set-1) / RowSum_i;  // joint probability 1->0 -> 0
          mat_f(0,1) = mat_a(time_t+1,1) * mat_b(time_t+1, 2) * mat_e(time_t,state_set-1) / RowSum_i;  // joint probability 1->0 -> 1
          break;

        case 3:
          mat_f(0,0) = mat_a(time_t+1,0) * mat_b(time_t+1, 1) * mat_e(time_t,state_set-1) / RowSum_i;  // joint probability 0->1 -> 0
          mat_f(0,1) = mat_a(time_t+1,1) * mat_b(time_t+1, 3) * mat_e(time_t,state_set-1) / RowSum_i;  // joint probability 0->1 -> 1
          break;

        case 4: //identical to case 3:
          mat_f(0,0) = mat_a(time_t+1,0) * mat_b(time_t+1, 1) * mat_e(time_t,state_set-1) / RowSum_i;  // joint probability 1->1 -> 0
          mat_f(0,1) = mat_a(time_t+1,1) * mat_b(time_t+1, 3) * mat_e(time_t,state_set-1) / RowSum_i;  // joint probability 1->1 -> 1
          break;
        }
      }

      if (time_t == (obs-1)){

        mat_g(time_t,state_set-1) =  mat_e(time_t,state_set-1); //smoothed_joint_state

      } else if (time_t == (obs-2)){

        mat_g(time_t,state_set-1) =  mat_f(0,0) + mat_f(0,1);  //smoothed_joint_state

      } else {

        //*************
        //  tau = t+2
        //*************

        RowSum_i = mat_d(time_t+2,0) + mat_d(time_t+2,1) + mat_d(time_t+2,2) + mat_d(time_t+2,3);
                                                                                            // Assume we are in state_set = 1:
        mat_c(time_t+2,0) = mat_a(time_t+2,0) * mat_b(time_t+2,0) * mat_f(0,0) / RowSum_i;  // joint probability 0->0->0 -> 0
        mat_c(time_t+2,1) = mat_a(time_t+2,0) * mat_b(time_t+2,1) * mat_f(0,1) / RowSum_i;  // joint probability 0->0->1 -> 0
        mat_c(time_t+2,2) = mat_a(time_t+2,1) * mat_b(time_t+2,2) * mat_f(0,0) / RowSum_i;  // joint probability 0->0->0 -> 1
        mat_c(time_t+2,3) = mat_a(time_t+2,1) * mat_b(time_t+2,3) * mat_f(0,1) / RowSum_i;  // joint probability 0->0->1 -> 1


        //*******************
        //  tau = t+3,...,T
        //*******************

        for(int i = (time_t+3); i < obs; i++){
          RowSum_i = mat_d(i,0) + mat_d(i,1) + mat_d(i,2) + mat_d(i,3);
          mat_c(i,0) =  mat_a(i,0) * mat_b(i,0) * (mat_c(i-1,0) + mat_c(i-1,1)) / RowSum_i;
          mat_c(i,1) =  mat_a(i,0) * mat_b(i,1) * (mat_c(i-1,2) + mat_c(i-1,3)) / RowSum_i;
          mat_c(i,2) =  mat_a(i,1) * mat_b(i,2) * (mat_c(i-1,0) + mat_c(i-1,1)) / RowSum_i;
          mat_c(i,3) =  mat_a(i,1) * mat_b(i,3) * (mat_c(i-1,2) + mat_c(i-1,3)) / RowSum_i;
        }

        mat_g(time_t,state_set-1) = mat_c(obs-1,0) + mat_c(obs-1,1) + mat_c(obs-1,2) + mat_c(obs-1,3);  //smoothed_joint_state
      }

    }
  }

  //=============
  // Step: 4
  //=============

  Rcpp::NumericMatrix mat_h(obs,2); // smoothed_marginal_state_prob: 'mat_h'

  mat_h(0,0) =  prob_initial_state_0;
  mat_h(0,1) =  prob_initial_state_1;

  for(int i = 1; i < obs; i++){
    mat_h(i, 0)  =  mat_g(i,0) + mat_g(i,1);
    mat_h(i, 1)  =  mat_g(i,2) + mat_g(i,3);
  }


  //=========================
  // Calculate LogLikelihood
  //=========================

  Rcpp::NumericVector like_t_T_temp(obs); // like_t_T_temp


                          // in previous versions, mat_a(1,0) and mat_a(1,1) were used incorrectly
    like_t_T_temp(0)  =   prob_initial_state_0 * (log(mat_a(0,0)) + log(prob_initial_state_0)) +
                          prob_initial_state_1 * (log(mat_a(0,1)) + log(prob_initial_state_1));


  for(int i = 1; i < obs; i++){
    like_t_T_temp(i)  =   mat_h(i,0) * log(mat_a(i,0)) +
                          mat_h(i,1) * log(mat_a(i,1)) +
                          ( mat_g(i,0) * log(mat_b(i,0)) +
                            mat_g(i,1) * log(mat_b(i,1)) +
                            mat_g(i,2) * log(mat_b(i,2)) +
                            mat_g(i,3) * log(mat_b(i,3)));
  }

  if (printLogLike > 0){
    double LogLikelihoodVal = std::accumulate(like_t_T_temp.begin(), like_t_T_temp.end(), 0.0);
    Rcpp::Rcout.precision(5);
    // Rcpp::Rcout << "param: " << std::fixed << param << std::endl;
    // Rcpp::Rcout << \"LogLikelihoodVal: \" << std::fixed << LogLikelihoodVal << std::endl;
    Rcpp::Rcout << "LogLikelihoodVal: " << std::fixed << LogLikelihoodVal << std::endl;
  }

  // sum the individual Log Likelihood contributions for each observation
  double LogLikelihoodVal = std::accumulate(like_t_T_temp.begin(), like_t_T_temp.end(), 0.0);

  return Rcpp::List::create(Rcpp::Named("like_t_T_temp", like_t_T_temp),  // individual observation's log likelihood contribution
                            Rcpp::Named("LogLikelihoodVal",  LogLikelihoodVal),  // sum of all log likelihood contributions
                            Rcpp::Named("mat_a",  mat_a),  // cond_density
                            Rcpp::Named("mat_b",  mat_b),  // trans_prob
                            Rcpp::Named("mat_c",  mat_c),  // joint_prob
                            Rcpp::Named("mat_d",  mat_d),  // joint_density
                            Rcpp::Named("mat_e",  mat_e),  // filtered_joined_state_prob
                            Rcpp::Named("mat_h",  mat_h),  // smoothed_marginal_state_prob
                            Rcpp::Named("mat_g",  mat_g)); // smoothed_joint_state

}



//=========================
//=========================
// Convenience Wrappers
//=========================
//=========================


//' LogLikeliHood
//'
//' \code{LogLikeliHood} is a convenience wrapper for function 'LogLikeliHood_FULL' that returns a single numeric scalar
//' with the log-likelihood of the Expectation-Maximization Maximum-Likelihood (EM-ML) algorithm as this is required by some
//' optimizers. This is a fast C++ implementation.
//' @export
// [[Rcpp::export]]
double LogLikeliHood(
    Rcpp::NumericVector param,  // param
    Rcpp::NumericMatrix datax,  // data
    int n_explanatory,          // number of explanatory variables governing 'process_explanatory' including constant [1..7]
    int n_transition,           // number of explanatory variables governing 'trans_p_explanatory' including constant [1,2]
    double prob_initial_state_0,// prob_initial_state_0: Probability that the first observation is in state '0'
    int printLogLike            // 0: Do not print output to console
){
    // Convenience wrapper to return a single log-likelihood value as required by certain optimizers
    return LogLikeliHood_FULL(param, datax, n_explanatory, n_transition, prob_initial_state_0, printLogLike)["LogLikelihoodVal"];
}



//' LogLikeliHood_byObs
//'
//' \code{LogLikeliHood_byObs} is a convenience wrapper for function 'LogLikeliHood_FULL' that returns all observations'
//' log-likelihood constributions as required by certain optimizers.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector LogLikeliHood_byObs(
    Rcpp::NumericVector param,  // param
    Rcpp::NumericMatrix datax,  // data
    int n_explanatory,          // number of explanatory variables governing 'process_explanatory' including constant [1..7]
    int n_transition,           // number of explanatory variables governing 'trans_p_explanatory' including constant [1,2]
    double prob_initial_state_0,// prob_initial_state_0: Probability that the first observation is in state '0'
    int printLogLike            // 0: Do not print output to console
){
    // Convenience wrapper to return single observation's log-likelihood constribution as required by certain optimizers
    // Rcpp::NumericVector LikelihoodContribution = LogLikeliHood_FULL(param3, datax3, printLogLike)["like_t_T_temp"];
    // return LikelihoodContribution;
    return LogLikeliHood_FULL(param, datax, n_explanatory, n_transition, prob_initial_state_0, printLogLike)["like_t_T_temp"];
}



//' LogLikeliHood_min
//'
//' \code{LogLikeliHood_min} is a convenience wrapper for function 'LogLikeliHood_FULL' that returns -1 * resulting log-likelihood
//' value as required by certain optimizers
//' @export
// [[Rcpp::export]]
double LogLikeliHood_min(
    Rcpp::NumericVector param,  // param
    Rcpp::NumericMatrix datax,  // data
    int n_explanatory,          // number of explanatory variables governing 'process_explanatory' including constant [1..7]
    int n_transition,           // number of explanatory variables governing 'trans_p_explanatory' including constant [1,2]
    double prob_initial_state_0,// prob_initial_state_0: Probability that the first observation is in state '0'
    int printLogLike            // 0: Do not print output to console
){
  // Convenience wrapper to return -1 * resulting log-likelihood value as required by certain optimizers
  double LogLike = LogLikeliHood_FULL(param, datax, n_explanatory, n_transition, prob_initial_state_0, printLogLike)["LogLikelihoodVal"];
  double res = -1 * LogLike;
  Rcpp::Rcout << "The negative of the printed LogLikelihoodVal is returned!";
  return res;
}




//' LogLikeliHood_RunTimeCheck
//'
//' \code{LogLikeliHood_RunTimeCheck} is a convenience wrapper for function 'LogLikeliHood_FULL' that accepts the number of repetitive log-likelihood
//' calculations as an additional input. This is used to test the systems's computational speed.
//' value as required by certain optimizers
//' @export
// [[Rcpp::export]]
double LogLikeliHood_RunTimeCheck(
    Rcpp::NumericVector param,  // param
    Rcpp::NumericMatrix datax,  // data
    int n_explanatory,          // number of explanatory variables governing 'process_explanatory' including constant [1..7]
    int n_transition,           // number of explanatory variables governing 'trans_p_explanatory' including constant [1,2]
    double prob_initial_state_0,// prob_initial_state_0: Probability that the first observation is in state '0'
    int printLogLike,           // 0: Do not print output to console
    int n_it                    // number of iterations
){
  for(int k = 0; k < n_it; k++){
     LogLikeliHood_FULL(param, datax, n_explanatory, n_transition, prob_initial_state_0, printLogLike)["like_t_T_temp"];
  }
  return n_it;
}


//---------------------------------------------------------------------------------------------------------------------------
// ReturnTiny
// Returns the smallest floating number that R can interpret on a regular machine.
// Last changed: 04/08/2017
//---------------------------------------------------------------------------------------------------------------------------

//' ReturnTiny
//'
//' \code{ReturnTiny} returns the smallest floating number that R can interpret on a regular machine.
//' @export
// [[Rcpp::export]]
double ReturnTiny(
){
  double tiny_1 = pow(2, -1074);
  return(log(tiny_1));
}
