# This R Script contains functions for doing cross validation of multiple races in the 2025 season

########################################################
# stint_CV gets RMSPE Cross validation results for Gaussian SSMs

stint_CV <- function(model_string, data_list) {
  
  # Store the full data/time series
  y <- data_list[["y"]]
  # Store the full compound vector
  compound <- data_list[["Compound"]]
  # Store the full pit vector
  pit <- data_list[["Pit"]]
  # Store the Fuel vector
  fuel <- data_list[["fuel_mass"]]
  
  
  # Get the index of laps that end a stint
  end_stint <- which(data_list[["Pit"]] == 1)
  # Add the last lap because that also marks the end of a stint
  end_stint <- c(end_stint,length(y))
  
  # Number of Stints
  num_stints <- length(end_stint)
  
  # Will get RMSPE for each stint so need an empty vector
  RMSPE_vec <- rep(NA, num_stints)
  
  # Will get plots for each stint
  plots <- vector("list", length = num_stints)
  
  posteriors <- vector("list", length = num_stints)
  
  for(j in 1:num_stints) {
    
    # Do CV on last quarter of each stint        
    # Number of folds
    if(j == 1) {
      K <- round(end_stint[j]/4)
    } else {
      K <- round((end_stint[j] - end_stint[j-1])/4)
    }    
    
    stan.models <- vector("list", length = K)
    posteriors[[j]] <- vector("list", length = K)
    
    # Vector to Store one-step ahead predictions and CI bounds
    ypred_hat <- rep(NA, K)
    ypred_l <- rep(NA, K)
    ypred_u <- rep(NA, K)
    
    for(i in 1:K) {
      
      # Update data vector to allow for next step prediction
      data_list[["TT"]] <- length(y[1:(i+end_stint[j]-K-1)])
      data_list[["y"]] <- y[1:(i+end_stint[j]-K-1)]
      data_list[["Compound"]] <- compound[1:(i+end_stint[j]-K-1)]
      c_map <- sort(unique(compound[1:(i+end_stint[j]-K-1)]))
      data_list[["compound_map"]] <- array(c_map, dim = length(c_map))
      data_list[["C_used"]] <- length(c_map)
      data_list[["Pit"]] <- pit[1:(i+end_stint[j]-K-1)]
      data_list[["fuel_mass"]] <- fuel[1:(i+end_stint[j]-K-1)]
      
      # Fit model on subset of data then get one step ahead prediction
      stan.models[[i]] <- stan(file = model_string,
                               data = data_list,
                               chains = 3, iter = 30000,
                               control = list(adapt_delta = .99, max_treedepth = 12))
      posteriors[[j]][[i]] <- extract(stan.models[[i]])
      ypred_hat[i]  <- mean(posteriors[[j]][[i]]$y_pred)
      ypred_l[i] <- quantile(posteriors[[j]][[i]]$y_pred, probs = .05)
      ypred_u[i] <- quantile(posteriors[[j]][[i]]$y_pred, probs = .95)
      
    }    
    
    # Plots
    df <- tibble(laps = 1:end_stint[j],
                 y = y[1:end_stint[j]])
    
    pred_df <- tibble(time = (end_stint[j] - K + 1):end_stint[j],
                      pred = ypred_hat,
                      pred_l = ypred_l,
                      pred_u = ypred_u)
    
    plots[[j]] <- ggplot() +
      geom_point(data = df,
                 aes(x = laps, y = y, color = "Observations")) +
      geom_point(data = pred_df,
                 aes(x = time, y = pred, color = "Predictions")) +
      geom_ribbon(data = pred_df, aes(x = time, ymin = pred_l, ymax = pred_u), fill = "black",color = "gray", alpha = .1) +
      scale_color_manual(values = c("Observations" = "red",
                                    "Predictions" = "blue")) +
      labs(title = paste("Scatter Plot of Lap vs Lap Times with Predictions for Stint", j),
           x = "Lap",
           y = "Laptime in Seconds",
           color = "Legend")
    
    # Calculate Summary Statistic
    MSPE <- mean((y[(end_stint[j]-K+1):end_stint[j]] - ypred_hat)^2)
    RMSPE_vec[j] <- sqrt(MSPE)
    
  }
  return(list(RMSPE_vec,plots,posteriors))
}

##################################################
# stint_CV is a function that takes in the model to use and a data list to pass to the stan code.
# It returns RMSPE, some plots, and the posterior samples for parameters obtained in stan

skewt_CV <- function(model_string, data_list) {
  
  # Store the full data/time series
  y <- data_list[["y"]]
  # Store the full compound vector
  compound <- data_list[["Compound"]]
  # Store the full pit vector
  pit <- data_list[["Pit"]]
  # Store the Fuel vector
  fuel <- data_list[["fuel_mass"]]
  
  
  # Get the index of laps that end a stint
  end_stint <- which(data_list[["Pit"]] == 1)
  # Add the last lap because that also marks the end of a stint
  end_stint <- c(end_stint,length(y))
  
  # Number of Stints
  num_stints <- length(end_stint)

  
  posteriors <- vector("list", length = num_stints)
  
  for(j in 1:num_stints) {
    
    # Do CV on last quarter of each stint        
    # Number of folds
    if(j == 1) {
      K <- round(end_stint[j]/4)
    } else {
      K <- round((end_stint[j] - end_stint[j-1])/4)
    }    
    
    stan.models <- vector("list", length = K)
    posteriors[[j]] <- vector("list", length = K)
    
    for(i in 1:K) {
      
      # Update data vector to allow for next step prediction
      data_list[["TT"]] <- length(y[1:(i+end_stint[j]-K-1)])
      data_list[["y"]] <- y[1:(i+end_stint[j]-K-1)]
#      data_list[["Compound"]] <- compound[1:(i+end_stint[j]-K-1)]
#     c_map <- sort(unique(compound[1:(i+end_stint[j]-K-1)]))
#      data_list[["compound_map"]] <- array(c_map, dim = length(c_map))
#      data_list[["C_used"]] <- length(c_map)
      data_list[["Pit"]] <- pit[1:(i+end_stint[j]-K-1)]
      data_list[["fuel_mass"]] <- fuel[1:(i+end_stint[j]-K-1)]
      
      # Fit model on subset of data then get one step ahead prediction
      stan.models[[i]] <- stan(file = model_string,
                               data = data_list,
                               chains = 2, iter = 10000,
                               seed = 1234,
                               control = list(adapt_delta = .99, max_treedepth = 12))
      posteriors[[j]][[i]] <- extract(stan.models[[i]])
      posteriors[[j]][[i]]$y_pred <- rsgt(n = length(posteriors[[j]][[i]]$z_pred), posteriors[[j]][[i]]$z_pred + posteriors[[j]][[i]]$gamma * fuel[ncol(posteriors[[j]][[i]]$z)+1], posteriors[[j]][[i]]$sdo, posteriors[[j]][[i]]$skew,2,2)
      
    }    
    
  }
  
  CRPS <- skewt_CRPS(posteriors, y)
  RMSPE <- skewt_RMSPE(posteriors, y)
  
  return(list(RMSPE,CRPS,posteriors))
}


##################################################
# arima_CV Does the same as above but it fits an Arima model instead

arima_cv <- function(y, pit) {
  
  # Get the index of laps that end a stint
  end_stint <- which(pit == 1)
  
  # Add First lap for laps that start a stint
  start_stint <- c(1, end_stint + 1)
  
  # Add the last lap because that also marks the end of a stint
  end_stint <- c(end_stint,length(y))
  
  # Number of Stints
  num_stints <- length(end_stint)
  
  # Will get RMSPE for each stint so need an empty vector
  RMSPE_vec <- rep(NA, num_stints)
  
  # Will get CRPS for each stint too
  CRPS_vec <- rep(NA, num_stints)
  
  # Will get plots for each stint
  plots <- vector("list", length = num_stints)
  
  for(i in 1:num_stints) {
    
    # Number of folds
    if(i == 1) {
      K <- round(end_stint[i]/4)
    } else {
      K <- round((end_stint[i] - end_stint[i-1])/4)
    }    
    
    arma.models <- vector("list", length = K)
    
    # Vector to Store one-step ahead predictions
    ypred_hat <- rep(NA, K)
    ypred_se <- rep(NA, K)
    
    for(j in 1:K) {
      
      arma.models[[j]] <- arima(y[start_stint[i]:(j+end_stint[i]-K-1)], order = c(1,0,0))
      pr <- predict(arma.models[[j]],n.ahead = 1)
      ypred_hat[j] <- as.numeric(pr$pred)
      ypred_se[j] <- as.numeric(pr$se)
      
    }
    
    # Calculate Summary Statistic
    MSPE <- mean((y[(end_stint[i]-K+1):end_stint[i]] - ypred_hat)^2)
    RMSPE_vec[i] <- sqrt(MSPE)
    CRPS_vec[i] <- mean(scoringRules::crps(y[(end_stint[i]-K+1):end_stint[i]], family = "normal", mean = ypred_hat, sd = ypred_se))
    
    
  }
  
  return(list(RMSPE_vec,CRPS_vec))
  
}


##########################################################
# A special function is needed for calculating the RMSPE of the skew t model
# posteriors is a list of posteriors obtained in skewt_CV

skewt_RMSPE <- function(posteriors,hamilton_laps){
  
  num_stints <- length(posteriors)
  
  RMSPE <- rep(NA,num_stints)
  
  plots <- vector("list", num_stints)
  
  for(i in 1:num_stints){
    
    temp <- rep(NA, length(posteriors[[i]]))
    ypred_hat <- rep(NA, length(posteriors[[i]]))
    
    pred_laps <- rep(NA, length(posteriors[[i]]))
    
    for(j in 1:length(posteriors[[i]])){
      
      #posteriors[[i]][[j]]$y_pred <- rsgt(n = length(posteriors[[i]][[j]]$z_pred), posteriors[[i]][[j]]$z_pred + posteriors[[i]][[j]]$gamma * fuel.kg[ncol(posteriors[[i]][[j]]$z)+1], posteriors[[i]][[j]]$sdo, posteriors[[i]][[j]]$skew,2,2)
      
      ypred_hat[j] <- mean(posteriors[[i]][[j]]$y_pred)
      
      temp[j] <- (hamilton_laps[ncol(posteriors[[i]][[j]]$z)+1] - mean(posteriors[[i]][[j]]$y_pred))^2
      
      pred_laps[j] <- ncol(posteriors[[i]][[j]]$z)+1
    }
    
    print(ypred_hat)
    
    RMSPE[i] <- sqrt(mean(temp))
  }
  
  return(RMSPE)
}

######################################################################
# skewt_CRPS calculates the CRPS for the posteriors obtained in skewt_CV
# Contrary to its name, it can be used for all the state-space models since it takes in a generic posterior


skewt_CRPS <- function(posteriors, hamilton_laps){
  
  num_stints <- length(posteriors)
  
  crps_skewt <- rep(NA,num_stints)
  
  for(i in 1:num_stints) {
    crps_temp <- rep(NA,length(posteriors[[i]]))
    for(j in 1:length(posteriors[[i]])) {
      crps_temp[j] <- scoringRules::crps_sample(y = hamilton_laps[ncol(posteriors[[i]][[j]]$z)+1], dat = as.vector(posteriors[[i]][[j]]$y_pred))
    }
    crps_skewt[i] <- mean(crps_temp)
  }
  
  return(crps_skewt)
}

###################################################################
# load_race_data creates the data_list that is passed in to our model
# for a specified race
#
# allLaps is passed implicitly from the global environment


load_race_data <- function(race) {
  
  hamilton_df <- allLaps %>%
    filter(Driver == "HAM", 
           is.na(PitOutTime) & is.na(PitInTime),
           Race == race) %>%
    mutate(Compound_Code = as.integer(factor(Compound,levels = c("HARD","MEDIUM","SOFT"))))
  
  hamilton_laps <- hamilton_df$LapTime
  
  hamilton_compound <- hamilton_df$Compound_Code
  
  hamilton_pit <- rep(0,length(hamilton_laps))
  
  for(i in 2:length(hamilton_df$Stint)) {
    if(hamilton_df$Stint[i] != hamilton_df$Stint[i-1]) {
      hamilton_pit[i-1] <- 1
    }    
  }                  
  
  fuel.kg <- seq(110,1,length.out = length(hamilton_laps))
  
  hamilton_map <- sort(unique(hamilton_compound))
  hamilton_used <- length(hamilton_map)
  
  # Data for stan code
  hamilton_data_base <- list(TT = length(hamilton_laps), 
                             y = hamilton_laps,
                             Pit = hamilton_pit,
                             v_reset0 = 0, 
                             sdo0 = .3, 
                             fuel_mass = fuel.kg,
                             C = 3,
                             Compound = hamilton_compound,
                             C_used = hamilton_used,
                             compound_map = hamilton_map)
  
  return(hamilton_data_base)
}



###################################################################
# main_CV will use the above functions to get both the RMSPE and CRPS for a single race
# for both AR(1) and skewt

main_CV <- function(model_string,data_list,racename){
  
  skewt_results <- skewt_CV(model_string = model_string, data_list = data_list)
  
  arima_results <- arima_cv(data_list$y,data_list$Pit)
  
  results <- tibble(skewt_RMSPE = skewt_results[[1]],
                    arima_RMSPE = arima_results[[1]],
                    skewt_CRPS = skewt_results[[2]],
                    arima_CRPS = arima_results[[2]])
  
  # Create filename
  filename <- paste0(racename, "_results2.csv")
  
  # Optional: make filename safe
  filename <- gsub(" ", "_", filename)
  
  # Write file
  write_csv(results, filename)
  
  summary_tbl <- tibble(
    Race = racename,
    
    skewt_RMSPE_total = sum(results$skewt_RMSPE),
    arima_RMSPE_total = sum(results$arima_RMSPE),
    
    skewt_RMSPE_mean = mean(results$skewt_RMSPE),
    arima_RMSPE_mean = mean(results$arima_RMSPE),
    
    skewt_CRPS  = mean(results$skewt_CRPS),
    arima_CRPS  = mean(results$arima_CRPS),
    
    n_stints = length(skewt_results[[1]])
  )
  
  return(summary_tbl)
  
}

