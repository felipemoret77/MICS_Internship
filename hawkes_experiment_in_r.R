options(warn=-1)
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

##Function to update the bins_vector for each simulation
update_bins_vector <- function(bins_vector, timestamp, delta){
  bin = ceiling(timestamp/delta) 
  bins_vector[bin] <- bins_vector[bin] + 1
  return(bins_vector)
}


## EXPERIMENT - TRUNCATING TIME. WHAT IF WE PERMIT ONLY VARIATIONS OF TIME OF A VALUE DELTA ?

##OBSERVATION1: Possible inputs to duplicates argument: "suppression" and "randomize" 
##OBSERVATION2: Possible inputs to opt_method: "mle" and "whittle" (for a discretized Hawkes proccess)
##OBSERVATION3: Possible inputs to init_estimation: "original_parameters" or "random"

experiment_truncate_time <- function(test_number, time_interval, number_of_simulations, baseline, alpha, beta, number_of_deltas, min_delta, max_delta, duplicates, opt_method, init_estimation){
  
  ##Vector of deltas for truncating time
  delta_vector = as.numeric(sprintf("%.2f", seq(from = min_delta, to = max_delta , length.out = number_of_deltas)))
  
  intensity_vector <- rep(c(0), number_of_simulations) 

  for(i in 1:number_of_simulations){
    
    original_process = hawkes(end = time_interval, fun = baseline, repr = alpha, family = "exp", rate = beta)
    original_process_timestamps <- original_process$p
    number_of_timestamps <- length(original_process_timestamps)
    intensity_vector[i] <- number_of_timestamps/time_interval
    
    ##CODING THE DELTA TRUNCATING
    for(truncate_delta_experiment in 1:number_of_deltas){
      
      if(delta_vector[truncate_delta_experiment] != 0){
        number_of_bins = as.integer(time_interval/delta_vector[truncate_delta_experiment])
        bins_vector <- rep(c(0), number_of_bins)
        
        for(k in 1:number_of_timestamps)
          bins_vector <- update_bins_vector(bins_vector, original_process_timestamps[k], delta_vector[truncate_delta_experiment])
        
        mean_timestamp_per_bin <- mean(bins_vector)
        
        list_to_csv_mean_bin_count <- list(test_number, i, delta_vector[truncate_delta_experiment], mean_timestamp_per_bin)
        
        ## Writing bin count mean results to a csv file
        write.table(list_to_csv_mean_bin_count,  
                    file = "bin_mean_count.csv", 
                    append = T, 
                    sep = ',', 
                    row.names = F, 
                    col.names = F)
      
      list_to_csv_bin_histograms <- list(test_number, i, delta_vector[truncate_delta_experiment], bins_vector)
      
      ## Writing bin histograms to a csv file
      write.table(list_to_csv_bin_histograms,  
                  file = "bin_count_histograms.csv", 
                  append = T, 
                  sep = ',', 
                  row.names = F, 
                  col.names = F)
    }
      
      if(opt_method == "mle" || opt_method == "both"){
        truncated_timestamps <- rep(c(0), number_of_timestamps) 
        for(k in 1:number_of_timestamps){
          truncated_timestamps[k] <- return_nearest_delta_multiple(original_process_timestamps[k], delta_vector[truncate_delta_experiment])
        }
        ##ELIMINATES DUPLICATES IF DUPLICATES == "suppression"
        if(duplicates == "supression")
          truncated_timestamps <- truncated_timestamps[!duplicated(truncated_timestamps)]
        ##RANDOMIZE TIMESTAMPS IF DUPLICATES == "randomize"
        else{
          random_uniform_numbers = runif(number_of_timestamps, min = 0, max = 1)*delta_vector[truncate_delta_experiment] 
          truncated_timestamps <- truncated_timestamps + random_uniform_numbers
        }
        if(init_estimation == "original_parameters")
          parameters_estimation_mle = try(mle(truncated_timestamps , "Exponential", original_process$end, opts = list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel" = 1.0e-8), init = c(baseline, alpha, beta)))
        else
          parameters_estimation_mle = try(mle(original_process_timestamps , "Exponential", original_process$end, opts = list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel" = 1.0e-8), init = NULL))
      }
      
      if(opt_method == "whittle" || opt_method == "both"){
        if(truncate_delta_experiment > 0){
          y = discrete(original_process, binsize = truncate_delta_experiment)
          if(init_estimation == "original_parameters")
            parameters_estimation_whittle  = try(whittle(y, "Exponential", init = c(baseline, alpha, beta)))
          else
            parameters_estimation_whittle  = try(whittle(y, "Exponential", init = NULL))
        }
      }
        
      if(opt_method == "both"){
        estimated_parameters_mle <- parameters_estimation_mle$par
        estimated_parameters_whittle <- parameters_estimation_whittle$par
        
        ## Writing MLE data to the csv file
        list_to_csv_mle <- list(test_number, i, alpha, beta, baseline, alpha/beta, delta_vector[truncate_delta_experiment], estimated_parameters_mle[2], estimated_parameters_mle[3], estimated_parameters_mle[1], estimated_parameters_mle[2]/estimated_parameters_mle[3], duplicates, "mle", init_estimation)
        
        ## Writing Whittle data to the csv file
        list_to_csv_whittle <- list(test_number, i, alpha, beta, baseline, alpha/beta, delta_vector[truncate_delta_experiment], estimated_parameters_whittle[2], estimated_parameters_whittle[3], estimated_parameters_whittle[1], estimated_parameters_whittle[2]/estimated_parameters_whittle[3], duplicates, "whittle", init_estimation)
        
        ## Writing MLE results to a csv file
        write.table(list_to_csv_mle ,  
                    file = "results.csv", 
                    append = T, 
                    sep = ',', 
                    row.names = F, 
                    col.names = F)
        
        ## Writing Whittle results to a csv file
        write.table(list_to_csv_whittle ,  
                    file = "results.csv", 
                    append = T, 
                    sep = ',', 
                    row.names = F, 
                    col.names = F)
      }
      
      if(opt_method == "mle" || opt_method == "whittle"){
        ## Getting estimated parameters
        if(opt_method == "mle")
          parameters_estimation <- parameters_estimation_mle
        else
          parameters_estimation <- parameters_estimation_whittle
          
        estimated_parameters <- parameters_estimation$par
      
        ## Writing data to the csv file
        list_to_csv <- list(test_number, i, alpha, beta, baseline, alpha/beta, delta_vector[truncate_delta_experiment], estimated_parameters[2], estimated_parameters[3], estimated_parameters[1], estimated_parameters[2]/estimated_parameters[3], duplicates, opt_method, init_estimation)
      
        ## Writing results to a csv file
        write.table(list_to_csv ,  
                    file = "results.csv", 
                    append = T, 
                    sep = ',', 
                    row.names = F, 
                    col.names = F)
      }
    }
  }
  print(mean(intensity_vector)) 
}


##Defining parameters for massive tests
number_of_tests <- 1
number_of_simulations_per_delta = 50
end_time <- 1000
number_of_deltas <- 10
estimation_method <- "mle"  ## Could be "mle" for maximum-likelihood maximization estimation or "whittle" for bin aggregation estimation or "both" for running the two methods in the same samples  
duplicates <- "randomize" ## Could be "randomize" for adding an appropriate uniform distribution random number to the timestamps or "supression" for removing duplicated timestamps
init_estimation <- "original_parameters" ## Could be "original_parameters"for using original baseline/alpha/beta or "random" for an arbitrary starting point when maximizing MLE

min_baseline <- 1
max_baseline <- 1
min_alpha <- 0.5
max_alpha <- 0.5
min_beta <- 1.5
max_beta <- 1.5
min_delta <- 0
max_delta <- 0.2

##Defining parameters vectors for different tests
"BaselineVector" = baseline_vector <- as.numeric(sprintf("%.2f", seq(from = min_baseline, to = max_baseline, length.out = number_of_tests)))
"AlphaVector" = alpha_vector <-  as.numeric(sprintf("%.2f", seq(from = min_alpha, to = max_alpha, length.out = number_of_tests)))
"BetaVector" = beta_vector <-  as.numeric(sprintf("%.2f", seq(from = min_beta, to = max_beta, length.out = number_of_tests)))

##Writing column names to the CSV file containing estimates of the parameters for each simulation
list_column_names_to_csv_parameters <- list("Test Index", "Simulation Index", "Real Alpha", "Real Beta", "Real Baseline", "Real Endogeneity", "Delta", "Estimated Alpha", "Estimated Beta", "Estimated Baseline", "Estimated Endogeneity", "Duplicates", "Estimation Method", "Initial Estimation Point")

write.table(list_column_names_to_csv_parameters,  
            file = "results.csv", 
            append = F, 
            sep =',', 
            row.names = F, 
            col.names = F)

##Writing column names to the CSV file containing the mean count of timestamps per bin per simulation per delta
list_column_names_to_csv_bin_mean_count <- list("Test Index", "Simulation Index", "Delta", "Mean Bin Count")

write.table(list_column_names_to_csv_bin_mean_count,  
            file = "bin_mean_count.csv", 
            append = F, 
            sep =',', 
            row.names = F, 
            col.names = F)

##Writing the column names to the CSV file containing the histograms of timestamps per bin per simulation per delta
list_column_names_to_csv_bin_count_histograms <- list("Test Index", "Simulation Index", "Delta", "Bin Count Histogram")

write.table(list_column_names_to_csv_bin_count_histograms,  
            file = "bin_count_histograms.csv", 
            append = F, 
            sep =',', 
            row.names = F, 
            col.names = F)


##Executing tests
for(test in 1:number_of_tests){
  experiment_truncate_time(test, end_time, number_of_simulations_per_delta, baseline_vector[test], alpha_vector[test], beta_vector[test], number_of_deltas, min_delta, max_delta, duplicates, estimation_method, init_estimation)
}

