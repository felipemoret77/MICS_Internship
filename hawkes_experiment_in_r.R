library(hawkesbow)
library(latex2exp)

##SUPPPORT FUNCTIONS

##Function to truncate time variation by a fixed delta
return_nearest_delta_multiple <- function(timestamp, delta){
  if(delta == 0){
    return(timestamp)
  }
  integer_quotient = as.integer(timestamp/delta)
  return(delta*integer_quotient)
}
  
## EXPERIMENT - TRUNCATING TIME. WHAT IF WE PERMIT ONLY VARIATIONS OF TIME OF A VALUE DELTA ?

##OBSERVATION1: Possible inputs to duplicates argument: "supression" and "randomize" 
##OBSERVATION2: Possible inputs to opt_method: "mle" and "whittle" (for a discretized Hawkes proccess)
##OBSERVATION3: Possible inputs to init_estimation: "original_parameters" or "random"

experiment_truncate_time <- function(test_number, time_interval, number_of_simulations, baseline, alpha, beta, number_of_deltas, min_delta, max_delta, duplicates, opt_method, init_estimation) {
  
  ##Vector of deltas for truncating time
  delta_vector = seq(from = min_delta, to = max_delta , length.out = number_of_deltas)
  
  for(i in 1:number_of_simulations){
    original_process = hawkes(end = time_interval, fun = baseline, repr = alpha, family = "exp", rate = beta)
    original_process_timestamps <- original_process$p
    number_of_timestamps <- length(original_process_timestamps)
    
    ##CODING THE DELTA TRUNCATING
    for(truncate_delta_experiment in 1:number_of_deltas){
      if(opt_method == "mle"){
        truncated_timestamps <- rep(c(0), number_of_timestamps) 
        print(delta_vector[truncate_delta_experiment])
        for(k in 1:number_of_timestamps)
          truncated_timestamps[k] <- return_nearest_delta_multiple(original_process_timestamps[k], delta_vector[truncate_delta_experiment])
          ##ELIMINATES DUPLICATES IF DUPLICATES == "suppression"
          if(duplicates == "supression")
            truncated_timestamps <- truncated_timestamps[!duplicated(truncated_timestamps)]
          ##RANDOMIZE TIMESTAMPS IF DUPLICATES == "randomize"
          else{
            random_uniform_numbers = runif(number_of_timestamps, min = 0, max = 1)*delta_vector[truncate_delta_experiment] 
            truncated_timestamps <- truncated_timestamps + random_uniform_numbers
          }
        if(init_estimation == "original_parameters")
          parameters_estimation = try(mle(truncated_timestamps , "Exponential", original_process$end, opts = list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel" = 1.0e-8), init = c(baseline, alpha, beta)))
        else
          parameters_estimation = try(mle(original_process_timestamps , "Exponential", original_process$end, opts = list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel" = 1.0e-8), init = NULL))
      }
      else{
        y = discrete(original_process, binsize = truncate_delta_experiment)
        if(init_estimation == "original_parameters")
          parameters_estimation  = try(whittle(y, "Exponential", init = c(baseline, alpha, beta)))
        else
          parameters_estimation  = try(whittle(y, "Exponential", init = NULL))
      }
      
      ## Getting estimated parameters
      estimated_parameters <- parameters_estimation$par
      
      ## Writing data to a csv file
      list_to_csv <- list(test_number, i, alpha, beta, baseline, alpha/beta, delta_vector[truncate_delta_experiment], estimated_parameters[2], estimated_parameters[3], estimated_parameters[1], estimated_parameters[2]/estimated_parameters[3], duplicates, opt_method, init_estimation)
      
      ## Writing results to a xlsx filw
      write.table(list_to_csv ,  
                  file = "results.csv", 
                  append = T, 
                  sep = ',', 
                  row.names = F, 
                  col.names = F)
    }
  }
}

##Defining parameters for massive tests
number_of_tests <- 1
number_of_simulations_per_delta = 1000
end_time <- 100
number_of_deltas <- 1
estimation_method <- "mle"  ## Could be "mle" for maximum-likelihood maximization estimation or "whittle" for bin aggregation estimation
duplicates <- "randomize" ## Could be "randomize" for adding an appropriate uniform distribution random number to the timestamps or "supression" for removing duplicated timestamps
init_estimation <- "original_parameters" ## Could be "original_parameters"for using original baseline/alpha/beta or "random" for an arbitrary starting point when maximizing MLE

min_baseline <- 1
max_baseline <- 1
min_alpha <- 0.5
max_alpha <- 0.5
min_beta <- 1.5
max_beta <- 1.5
min_delta <- 0
max_delta <- 0

##Defining parameters vectors for different tests
"BaselineVector" = baseline_vector <- seq(from = min_baseline, to = max_baseline , length.out = number_of_tests)
"AlphaVector" = alpha_vector <-  seq(from = min_alpha, to = max_alpha , length.out = number_of_tests)
"BetaVector" = beta_vector <-  seq(from = min_beta, to = max_beta , length.out = number_of_tests)

##Writing column names to the CSV file
list_column_names_to_csv <- list("Test Index", "Simulation Index", "Real Alpha", "Real Beta", "Real Baseline", "Real Endogeneity", "Delta", "Estimated Alpha", "Estimated Beta", "Estimated Baseline", "Estimated Endogeneity", "Duplicates", "Estimation Method", "Initial Estimation Point")

write.table(list_column_names_to_csv,  
            file = "results.csv", 
            append = F, 
            sep =',', 
            row.names = F, 
            col.names = F)

##Executing tests
for(test in 1:number_of_tests){
  experiment_truncate_time(test, end_time, number_of_simulations_per_delta, baseline_vector[test], alpha_vector[test], beta_vector[test], number_of_deltas, min_delta, max_delta, duplicates, estimation_method, init_estimation)
}

