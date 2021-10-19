library(hawkesbow)
library(latex2exp)
library(xlsx)

##SUPPPORT FUNCTIONS

return_vector_mean_from_vector_of_vectors <- function(vector_of_vectors){
  number_of_vectors <- length(vector_of_vectors)
  number_of_elements_in_each_vector <- length(vector_of_vectors[1])
  mean_vector <- c(1:number_of_elements_in_each_vector)
  
  for(i in 1:number_of_elements_in_each_vector){
    sum_variable <- 0.0
    for(j in 1:number_of_vectors){
      sum_variable <- sum_variable + vector_of_vectors[j][i]
    }
    mean_property_i <- sum_variable/number_of_vectors 
    mean_vector[i] <- mean_property_i
  }
  return(mean_vector)
}

##Function to truncate time variation by a fixed delta
return_nearest_delta_multiple <- function(timestamp, delta){
  nearest_truncate_timestamp <- 0.0
  while(nearest_truncate_timestamp <= timestamp){
    nearest_truncate_timestamp <- nearest_truncate_timestamp + delta
  }
  return(nearest_truncate_timestamp)
}
  

## EXPERIMENT - TRUNCATING TIME. WHAT IF WE PERMIT ONLY VARIATIONS OF TIME OF A VALUE DELTA ?

##OBSERVATION1: Possible inputs to duplicates argument: "supression" and "randomize" 
##OBSERVATION2: Possible inputs to opt_method: "mle" and "whittle" (for a discretized hawkes proccess)

experiment_truncate_time <- function(time_interval, number_of_simulations, baseline, alpha, beta, number_of_deltas, min_delta, max_delta, duplicates, opt_method) {
  
  ##Vector of deltas for truncating time
  delta_vector = seq(from = min_delta, to = max_delta , length.out = number_of_deltas)
  
  ##Vectors for the parameters (and mean of the parameters) for the simulations when delta = 0
  zero_delta_MLE_baseline_vector <- rep(c(0), number_of_simulations)
  zero_delta_MLE_alpha_vector <- rep(c(0), number_of_simulations)
  zero_delta_MLE_beta_vector <- rep(c(0), number_of_simulations)
  zero_delta_MLE_parameters_mean <- rep(c(0), 3)
  
  ##Matrices for storing results for different delta truncated timestamps simulations
  non_zero_delta_MLE_baseline_matrix <- matrix(0, number_of_deltas, number_of_simulations)
  non_zero_delta_MLE_alpha_matrix <- matrix(0, number_of_deltas, number_of_simulations)
  non_zero_delta_MLE_beta_matrix <- matrix(0, number_of_deltas, number_of_simulations) 
  
  ##Vectors for storing the mean result for each delta truncated timestamp experiment
  non_zero_delta_MLE_parameters_mean <- matrix(0, number_of_deltas, 3) 
  
  for(i in 1:number_of_simulations){
    original_process = hawkes(end = time_interval, fun = baseline, repr = alpha, family = "exp", rate = beta)
    original_process_timestamps <- original_process$p
    if(opt_method == "mle")
      zero_delta_MLE_estimation = mle(original_process_timestamps , "Exponential", original_process$end, opts = list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8))
    else{
      y = discrete(original_process, binsize = 1)
      zero_delta_MLE_estimation = whittle(y, "Exponential")
    }
    zero_delta_MLE_parameters <- zero_delta_MLE_estimation$par
    zero_delta_MLE_baseline_vector[i] <- zero_delta_MLE_parameters[1]
    zero_delta_MLE_alpha_vector[i] <- zero_delta_MLE_parameters[2]
    zero_delta_MLE_beta_vector[i] <- zero_delta_MLE_parameters[3]
    number_of_timestamps <- length(original_process_timestamps)
    
    ##CODING THE DELTA TRUNCATING
    for(truncate_delta_experiment in 1:number_of_deltas){
      truncated_timestamps <- rep(c(0), number_of_timestamps) 
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
      if(opt_method == "mle")
        current_delta_MLE_estimation = mle(truncated_timestamps , "Exponential", original_process$end, opts = list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8))
      else{
        original_process$p <- truncated_timestamps  
        y = discrete(original_process, binsize = 1)
        current_delta_MLE_estimation  = whittle(y, "Exponential")
      }
      current_delta_MLE_parameters <- current_delta_MLE_estimation$par
      non_zero_delta_MLE_baseline_matrix[truncate_delta_experiment, i] = current_delta_MLE_parameters[1]
      non_zero_delta_MLE_alpha_matrix[truncate_delta_experiment, i] = current_delta_MLE_parameters[2]
      non_zero_delta_MLE_beta_matrix[truncate_delta_experiment, i] = current_delta_MLE_parameters[3]
    }
  }
      
  ##EVALUATING MEAN RESULTS
  zero_delta_MLE_parameters_mean[1] <- mean(zero_delta_MLE_baseline_vector)
  zero_delta_MLE_parameters_mean[2] <- mean(zero_delta_MLE_alpha_vector)
  zero_delta_MLE_parameters_mean[3] <- mean(zero_delta_MLE_beta_vector)
  
  for(truncate_delta_experiment in 1:number_of_deltas){
    non_zero_delta_MLE_parameters_mean[truncate_delta_experiment, 1] = mean(non_zero_delta_MLE_baseline_matrix[truncate_delta_experiment,]) 
    non_zero_delta_MLE_parameters_mean[truncate_delta_experiment, 2] = mean(non_zero_delta_MLE_alpha_matrix[truncate_delta_experiment,])
    non_zero_delta_MLE_parameters_mean[truncate_delta_experiment, 3] = mean(non_zero_delta_MLE_beta_matrix[truncate_delta_experiment,])
  }

  ##print("The estimation of the parameters for Delta = 0 is:\n")
  ##print(zero_delta_MLE_parameters_mean)
  
  ##print("The estimation of the parameters for the different delta values are: ")
  
  ##for(truncate_delta_experiment in 1:number_of_deltas){
  ##  print(paste("Delta: ", delta_vector[truncate_delta_experiment]))
  ##  print(paste("Baseline: ", non_zero_delta_MLE_parameters_mean[truncate_delta_experiment, 1]))
  ##  print(paste("Alpha: ", non_zero_delta_MLE_parameters_mean[truncate_delta_experiment, 2]))
  ##  print(paste("Beta: ", non_zero_delta_MLE_parameters_mean[truncate_delta_experiment, 3]))
  ##}
  
  results <- list(zero_delta_MLE_parameters_mean, non_zero_delta_MLE_parameters_mean, delta_vector)

  return(results)
}

## EXPERIMENT - TRANSLATING ALL THE TIME SERIES BY A FIXED AMOUNT DELTA

##OBSERVATION1: Possible inputs to opt_method: "mle" and "whittle" (for a discretized hawkes proccess)

experiment_translating_timeseries <- function(time_interval, number_of_simulations, baseline, alpha, beta, number_of_deltas, min_delta, max_delta, opt_method) {
  
  ##Vector of deltas for truncating time
  delta_vector = seq(from = min_delta, to = max_delta , length.out = number_of_deltas)
  
  ##Vectors for the parameters (and mean of the parameters) for the simulations when delta = 0
  zero_delta_MLE_baseline_vector <- rep(c(0), number_of_simulations)
  zero_delta_MLE_alpha_vector <- rep(c(0), number_of_simulations)
  zero_delta_MLE_beta_vector <- rep(c(0), number_of_simulations)
  zero_delta_MLE_parameters_mean <- rep(c(0), 3)
  
  ##Matrices for storing results for different delta trasnlated timestamps simulations
  non_zero_delta_MLE_baseline_matrix <- matrix(0, number_of_deltas, number_of_simulations)
  non_zero_delta_MLE_alpha_matrix <- matrix(0, number_of_deltas, number_of_simulations)
  non_zero_delta_MLE_beta_matrix <- matrix(0, number_of_deltas, number_of_simulations) 
  
  ##Vectors for storing the mean result for each delta truncated timestamp experiment
  non_zero_delta_MLE_parameters_mean <- matrix(0, number_of_deltas, 3) 
  
  for(i in 1:number_of_simulations){
    original_process = hawkes(end = time_interval, fun = baseline, repr = alpha, family = "exp", rate = beta)
    original_process_timestamps <- original_process$p
    if(opt_method == "mle")
      zero_delta_MLE_estimation = mle(original_process_timestamps , "Exponential", original_process$end)
    else{
      y = discrete(original_process, binsize = 1)
      zero_delta_MLE_estimation = whittle(y, "Exponential")
    }
    zero_delta_MLE_parameters <- zero_delta_MLE_estimation$par
    zero_delta_MLE_baseline_vector[i] <- zero_delta_MLE_parameters[1]
    zero_delta_MLE_alpha_vector[i] <- zero_delta_MLE_parameters[2]
    zero_delta_MLE_beta_vector[i] <- zero_delta_MLE_parameters[3]
    number_of_timestamps <- length(original_process_timestamps)
    
    ##CODING THE DELTA TRANSLATION
    for(translate_delta_experiment in 1:number_of_deltas){
      translated_timestamps <- rep(c(0), number_of_timestamps) 
      translated_timestamps <- original_process_timestamps + rep(c(translate_delta_experiment), number_of_timestamps)
      if(opt_method == "mle")
        current_delta_MLE_estimation = mle(translated_timestamps , "Exponential", original_process$end)
      else{
        original_process$p <- translated_timestamps 
        y = discrete(original_process, binsize = 1)
        urrent_delta_MLE_parameters  = whittle(y, "Exponential")
      }
      current_delta_MLE_parameters <- current_delta_MLE_estimation$par
      non_zero_delta_MLE_baseline_matrix[translate_delta_experiment, i] = current_delta_MLE_parameters[1]
      non_zero_delta_MLE_alpha_matrix[translate_delta_experiment, i] = current_delta_MLE_parameters[2]
      non_zero_delta_MLE_beta_matrix[translate_delta_experiment, i] = current_delta_MLE_parameters[3]
    }
  }
  
  ##EVALUATING MEAN RESULTS
  zero_delta_MLE_parameters_mean[1] <- mean(zero_delta_MLE_baseline_vector)
  zero_delta_MLE_parameters_mean[2] <- mean(zero_delta_MLE_alpha_vector)
  zero_delta_MLE_parameters_mean[3] <- mean(zero_delta_MLE_beta_vector)
  
  for(translate_delta_experiment in 1:number_of_deltas){
    non_zero_delta_MLE_parameters_mean[translate_delta_experiment, 1] = mean(non_zero_delta_MLE_baseline_matrix[translate_delta_experiment,]) 
    non_zero_delta_MLE_parameters_mean[translate_delta_experiment, 2] = mean(non_zero_delta_MLE_alpha_matrix[translate_delta_experiment,])
    non_zero_delta_MLE_parameters_mean[translate_delta_experiment, 3] = mean(non_zero_delta_MLE_beta_matrix[translate_delta_experiment,])
  }
  
  ##print("The estimation of the parameters for Delta = 0 is:\n")
  ##print(zero_delta_MLE_parameters_mean)
  
  ##print("The estimation of the parameters for the different delta values are: ")
  
  ##for(translate_delta_experiment in 1:number_of_deltas){
  ##  print(paste("Delta: ", delta_vector[translate_delta_experiment]))
  ##  print(paste("Baseline: ", non_zero_delta_MLE_parameters_mean[translate_delta_experiment, 1]))
  ##  print(paste("Alpha: ", non_zero_delta_MLE_parameters_mean[translate_delta_experiment, 2]))
  ##  print(paste("Beta: ", non_zero_delta_MLE_parameters_mean[translate_delta_experiment, 3]))
  ##}
  
  results <- list(zero_delta_MLE_parameters_mean, non_zero_delta_MLE_parameters_mean, delta_vector)
  
  return(results)
}
 

##Defining parameters for massive tests
number_of_tests <- 1
number_of_simulations_per_delta = 10
end_time <- 1000
number_of_deltas <- 10
estimation_method <- "whittle"  ##Could be "mle" for maximum-likelihood maximization estimation
duplicates <- "suppression" ##Could be "randomize" for adding an appropriate uniform distribution random number to the timestamps

min_baseline <- 1
max_baseline <- 10
min_alpha <- 0.1
max_alpha <- 1
min_beta <- 1
max_beta <- 2
min_delta <- 0.01
max_delta <- 0.1

##Defining parameters vectors for different tests
"BaselineVector" = baseline_vector <- seq(from = min_baseline, to = max_baseline , length.out = number_of_tests)
"AlphaVector" = alpha_vector <-  seq(from = min_alpha, to = max_alpha , length.out = number_of_tests)
"BetaVector" = beta_vector <-  seq(from = min_beta, to = max_beta , length.out = number_of_tests)
"DeltaVector" =  delta_vector <- seq(from = min_delta, to = max_delta , length.out = number_of_deltas) 

##Defining vectors for storing results
##results_truncate_time <- rep(c(0), number_of_tests)
##zero_delta_MLE_parameters_mean_truncate <- rep(c(0), number_of_tests)
##non_zero_delta_MLE_parameters_mean_truncate <- rep(c(0), number_of_tests) 

##Executing tests and storing/plotting results
for(test in 1:number_of_tests){
  results_truncate_time <- experiment_truncate_time(end_time, number_of_simulations_per_delta, baseline_vector[test], alpha_vector[test], beta_vector[test], number_of_deltas, min_delta, max_delta, duplicates, estimation_method)
  zero_delta_MLE_parameters_mean_truncate <- results_truncate_time[1]
  non_zero_delta_MLE_parameters_mean_truncate <- results_truncate_time[2]
  plot(delta_vector, non_zero_delta_MLE_parameters_mean_truncate[[1]][,1],  xlab = TeX("$\\delta$"), ylab= TeX("$\\lambda_{0}$"), pch = 15, col = "red", type = "o", main = paste("Real lambda0 = ", baseline_vector[test]))
  plot(delta_vector, non_zero_delta_MLE_parameters_mean_truncate[[1]][,2],  xlab = TeX("$\\delta$"), ylab= TeX("$\\alpha$"), pch = 15, col = "red", type = "o", main = paste("Real alpha = ", alpha_vector[test]))
  plot(delta_vector, non_zero_delta_MLE_parameters_mean_truncate[[1]][,3],  xlab = TeX("$\\delta$"), ylab= TeX("$\\beta$"), pch = 15, col = "red", type = "o", main = paste("Real beta = ", beta_vector[test]))
  plot(delta_vector, non_zero_delta_MLE_parameters_mean_truncate[[1]][,2]/non_zero_delta_MLE_parameters_mean_truncate[[1]][,3],  xlab = TeX("$\\delta$"), ylab= TeX("$\\alpha \\ \\beta$"), pch = 15, col = "red", type = "b", main = paste("Real endogeneity = ", alpha_vector[test]/beta_vector[test]))
  ##List to wirte to xlsx file
  list_to_xlsx <- list("Baseline " = non_zero_delta_MLE_parameters_mean_truncate[[1]][,1], "Alpha" = non_zero_delta_MLE_parameters_mean_truncate[[1]][,2], "Beta" = non_zero_delta_MLE_parameters_mean_truncate[[1]][,3], "Endogeneity" = non_zero_delta_MLE_parameters_mean_truncate[[1]][,2]/non_zero_delta_MLE_parameters_mean_truncate[[1]][,3])
  ##Writing results to a xlsx filw
  if(test == 1){
    write.xlsx2(list_to_xlsx, file = "results.xlsx", sheetName = paste("TEST", test), append = FALSE, col.names = TRUE )
  }
  else{
    write.xlsx2(list_to_xlsx, file = "results.xlsx", sheetName = paste("TEST", test), append = TRUE, col.names = TRUE )
  }
}

write.xlsx2(delta_vector, file = "results.xlsx", sheetName = "DeltaValues", append = TRUE, col.names = TRUE)
write.xlsx2(alpha_vector, file = "results.xlsx", sheetName = "RealAlphaValues", append = TRUE, col.names = TRUE)
write.xlsx2(beta_vector, file = "results.xlsx", sheetName = "RealBetaValues", append = TRUE, col.names = TRUE)
write.xlsx2(baseline_vector, file = "results.xlsx", sheetName = "RealBaselineIntensities", append = TRUE, col.names = TRUE)



##results_truncate_time <- experiment_truncate_time(1000, 100, 1, 0.5, 2, 30, 0.01 ,0.05, "randomize", "mle")
##zero_delta_MLE_parameters_mean_truncate <- results_truncate_time[1]
##non_zero_delta_MLE_parameters_mean_truncate <- results_truncate_time[2]
##delta_vector_truncate <- results_truncate_time[3]
##plot(delta_vector_truncate[[1]], non_zero_delta_MLE_parameters_mean_truncate[[1]][,1],  xlab = TeX("$\\delta$"), ylab= TeX("$\\lambda_{0}$"), pch = 15, col = "red", type = "o")
##plot(delta_vector_truncate[[1]], non_zero_delta_MLE_parameters_mean_truncate[[1]][,2],  xlab = TeX("$\\delta$"), ylab= TeX("$\\alpha$"), pch = 15, col = "red", type = "o")
##plot(delta_vector_truncate[[1]], non_zero_delta_MLE_parameters_mean_truncate[[1]][,3],  xlab = TeX("$\\delta$"), ylab= TeX("$\\beta$"), pch = 15, col = "red", type = "o")
##plot(delta_vector_truncate[[1]], non_zero_delta_MLE_parameters_mean_truncate[[1]][,2]/non_zero_delta_MLE_parameters_mean_truncate[[1]][,3],  xlab = TeX("$\\delta$"), ylab= TeX("$\\alpha \\ \\beta$"), pch = 15, col = "red", type = "b")


##LOOPING FOR RESULTS OF THE TRANSLATING EXPERIMENT (EXPERIMENT 2)
##results_translate_time <- experiment_translating_timeseries(100, 100, 1, 0.5, 2, 20, 0.00001 ,000.1)
##zero_delta_MLE_parameters_mean_translate <- results_translate_time[1]
##non_zero_delta_MLE_parameters_mean_translate <- results_translate_time[2]
##delta_vector_translate <- results_translate_time[3]
##plot(delta_vector_translate[[1]], non_zero_delta_MLE_parameters_mean_translate[[1]][,1],  xlab = TeX("$\\delta$"), ylab= TeX("$\\lambda_{0}$"), pch = 15, col = "red", type = "b")
##plot(delta_vector_translate[[1]], non_zero_delta_MLE_parameters_mean_translate[[1]][,2],  xlab = TeX("$\\delta$"), ylab= TeX("$\\alpha$"), pch = 15, col = "red", type = "b")
##plot(delta_vector_translate[[1]], non_zero_delta_MLE_parameters_mean_translate[[1]][,3],  xlab = TeX("$\\delta$"), ylab= TeX("$\\beta$"), pch = 15, col = "red", type = "b")
##plot(delta_vector_translate[[1]], non_zero_delta_MLE_parameters_mean_translate[[1]][,2]/non_zero_delta_MLE_parameters_mean_translate[[1]][,3],  xlab = TeX("$\\delta$"), ylab= TeX("$\\alpha \\ \\beta$"), pch = 15, col = "red", type = "b")


