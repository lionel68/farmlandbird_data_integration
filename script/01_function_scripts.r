#################
# Helper function to fit the models presented in:
# Hertzog et al Model-based data integration
# problems can be reported by mailing:
# lionel.hertzog@thuenen.de
#############

# the functions take as inputs:
# @species, a character, the species name
# @path, a character, the path to save the files
# ..., additional parameters to pass to the fitting function

# the functions returns as files saved on disk:
# the fitted model object
# the R2 and RMSE
# the Rhat of the model parameters
# the predicted annual abundance index for the CBBS routes
# two posterior predictive checks graphs

# to run the functions assume that the following data frames
# are loaded in the environment:
# mhb, ornitho_list, ornitho_unstr, list_all and unstr_all
# see the example script for the required structure of the 
# data frames

## function to fit model 1

fit_model1 <- function(species, path = ""){
  
  # filter mhb data
  dat_mhb <-  filter(mhb, Artname == species)
  # turn the year effect into a factor
  dat_mhb$fJahr <- factor(dat_mhb$Jahr)
  
  m <- trim(count ~ ROUTENCODE + Jahr, dat_mhb, overdisp = TRUE, verbose = TRUE, model = 3)
  # re-fit the model without overdisp if needed
  if(overdispersion(m) < 1.0){
    m <- trim(count ~ ROUTENCODE + Jahr, dat_mhb, overdisp = FALSE, verbose = TRUE, model = 3)
  }
  # the total breeding bird abundance over all Mhb plots
  # with 95% confidence band
  tt <- totals(m, level = 0.95) # save this
  
  write.table(tt, paste0(path, species, "/totals_modeltrim.csv"), sep = ",", row.names = FALSE)
  
  # rmse and r square for TRIM
  fitted <- results(m)
  # add the actual obs
  ff <- merge(fitted, dat_mhb, by.x = c("site", "time"), by.y = c("ROUTENCODE", "Jahr"))
  # the residuals
  ff$res <- ff$count - ff$fitted
  # the normalized root mean squared error
  rmse <- sqrt(mean(ff$res ** 2, na.rm = TRUE)) / mean(dat_mhb$count) # save this
  # R2
  r2 <- 1 - (sum(ff$res ** 2, na.rm = TRUE) / sum((ff$count - mean(ff$count)) ** 2)) # save this
  
  r2_rmse <- data.frame(metric = c("r2", "rmse"), value = c(r2, rmse))
  
  write.table(r2_rmse, paste0(path, species, "/r2_modeltrim.csv"),
              sep = ",", row.names = FALSE)
  
}

## function to fit model 2

fit_model2 <- function(species, path = "", ...){
  
  # filter Mhb data
  dat_mhb <-  filter(mhb, Artname == species)
  
  # set the model matrix
  dat_mhb$fJahr <- factor(dat_mhb$Jahr)
  modmat <- model.matrix(~fJahr, dat_mhb, contrasts.arg = list(fJahr = "contr.sum"))
  route_id <- as.numeric(factor(dat_mhb$ROUTENCODE))
  
  # priors
  intercept <- student(3, 0, 10)
  yr_eff <- normal(0, 1, dim = ncol(modmat) - 1)
  sd_route <- exp(student(3, 0, 10))
  alpha <- normal(0, 1, dim= max(route_id))
  route_avg <- intercept + alpha * sd_route
  sd_res <- gamma(0.01, 0.01)
  
  # linear predictor
  mu <- exp(modmat[,-1] %*% yr_eff + route_avg[route_id])
  # the response
  prob <- sd_res / (sd_res + mu)
  
  count <- as_data(dat_mhb$count)
  distribution(count) <- negative_binomial(size = sd_res, prob = prob)
  
  # set the model
  m_1 <- model(yr_eff, sd_route, sd_res, intercept)
  
  # run the model
  t1_1 <- Sys.time()
  d_1 <- mcmc(m_1, ...)
  t1_2 <- Sys.time()
  
  
  # save the model (save other stuff for later use?)
  saveRDS(d_1, paste0(path, "", species, "/fit_model1.rds"))
  cat(paste0("Fitting time for ", species,
             ": ", round(as.numeric(difftime(t1_2, t1_1, units = "min")), 3),
             " min\n"))
  
  # compute rhat
  tmp <- as.data.frame(coda::gelman.diag(d_1)$psrf)
  tmp$species <- species
  names(tmp)[1:2] <- c("point", "upper")
  tmp$param <- rownames(tmp)
  rownames(tmp) <- NULL
  tmp$ESS <- coda::effectiveSize(d_1)
  write.table(tmp, paste0(path,  species, "/rhat_model1.csv"),
              sep = ",", row.names = FALSE)
  
  cat(paste0("Maximum rhat value is: ", tmp[which.max(tmp$point),"point"], "\n"))
  

  ## r2 and rmse computation
  r2 <- greta.checks::bayes_R2(count, mu, d_1)
  error <- greta.checks::rmse(count, mu, d_1, norm = TRUE)
  # put the two together
  r2_error <- as.data.frame(rbind(r2, error))
  r2_error$metric <- c("r2", "rmse")
  
  write.table(r2_error, paste0(path,  species, "/rmse_model1.csv"),
              sep = ",", row.names = FALSE)
  
  
  # now get se from the predictions between 2012 and 2017
  # in a manner similar to TRIM, we compute the totals
  
  # new data frame to compute
  newdat <- expand.grid(fJahr = factor(2005:2018, levels = 2005:2018),
                        route_id = unique(route_id))
  # predicted model matrix
  modmat_pred <- model.matrix(~ fJahr, 
                              newdat,
                              contrasts.arg = list(fJahr = "contr.sum"))
  mu_pred <- exp(modmat_pred[,-1] %*% yr_eff + route_avg[newdat$route_id])
  # calculate expectations
  d_pred <- t(as.matrix(calculate(mu_pred, values = d_1)))
  # some data wraggling
  d_df <- as.data.frame(d_pred)
  d_df$Jahr <- newdat$fJahr
  
  tt <- gather(d_df, key = "iteration", value = "val", -Jahr)
  
  tt %>%
    group_by(Jahr, iteration) %>%
    summarise(S = sum(val)) %>%
    ungroup() %>%
    group_by(Jahr) %>%
    summarise(M = mean(S), SD = sd(S),
              lo = quantile(S, probs = 0.025),
              hi = quantile(S, probs = 0.975)) -> tt_d
  
  write.table(tt_d, paste0(path,  species, "/totals_model1.csv"),
              sep = ",", row.names = FALSE)
  # get posterior predictive checks
  pp1 <- greta.checks::pp_check(d_1, count, nsim = 50, type = "dens_overlay")
  pp2 <- greta.checks::pp_check(d_1, count, nsim = 50, type = "error_scatter_avg")
  
  ggsave(paste0(path,  species, "/pp1_model1.png"), pp1)
  ggsave(paste0(path,  species, "/pp2_model1.png"), pp2)
  
}


## function to fit model 3

fit_model3 <- function(species,path="", ...){
  
  # filter Mhb data
  dat_mhb <-  filter(mhb, Artname == species)
  
  ## ornitho list
  dat_ornitho_list <- subset(ornitho_list, SPECIES_NAME == species)
 
  ## set the model matrices
  dat_mhb$fJahr2 <- factor(dat_mhb$Jahr, levels = 2005:2019)
  dat_ornitho_list$fyear <- factor(dat_ornitho_list$year, levels = 2005:2019)
  
  yr_eff_mhb <- model.matrix(~fJahr2, dat_mhb,
                             contrasts.arg = list(fJahr2 = "contr.sum"))
  yr_eff_list <- model.matrix(~fyear, dat_ornitho_list,
                              contrasts.arg = list(fyear = "contr.sum"))
  
  # prior for common year effect
  yr_eff <- normal(0, 1, dim = ncol(yr_eff_mhb) - 1)
  
  ## model for mhb
  # index for routes
  route_id <- as.numeric(factor(dat_mhb$ROUTENCODE))
  
  # priors 
  intercept_mhb <- student(3, 0, 10)
  sd_route <- exp(student(3, 0, 10))
  alpha_route <- normal(0, 1, dim = max(route_id))
  route_avg <- intercept_mhb + alpha_route * sd_route
  sd_res_mhb <- gamma(0.01, 0.01)
  
  # linear predictor
  mu_mhb <- exp(yr_eff_mhb[,-1] %*% yr_eff + route_avg[route_id])
  
  # re-param for negative binomial
  prob_mhb <- sd_res_mhb / (mu_mhb + sd_res_mhb)
  
  ## model for ornitho lists
  # index for locality
  locality_id <- as.numeric(factor(dat_ornitho_list$ID_PLACE))
  # prior fixed effect
  intercept_list <- student(3, 0, 10)

  # prior for random effect
  sd_locality <- exp(student(3, 0, 10))
  alpha_locality <- normal(0, 1, dim = max(locality_id))
  locality_avg <- intercept_list + alpha_locality * sd_locality
  
  # prior for zero-inflation
  theta <- beta(1, 1) # uniform between 0 and 1
  
  
  mu_list <- dat_ornitho_list$effort_h * exp(yr_eff_list[,-1] %*% yr_eff + locality_avg[locality_id])
  
  # the response
  y_mhb <- as_data(dat_mhb$count)
  y_ornitho <- as_data(dat_ornitho_list$count)
  
  distribution(y_mhb) <- negative_binomial(size = sd_res_mhb,
                                           prob = prob_mhb)
  distribution(y_ornitho) <- greta.multivariate::zero_inflated_poisson(theta, 
                                                                       mu_list,
                                                                       dim = length(y_ornitho))
  
  m_2 <- model(intercept_mhb, intercept_list,
               yr_eff, sd_route, sd_locality,
               sd_res_mhb, theta)
  
  cat("Starting fitting model 2")
  
  t2_1 <- Sys.time()
  d_2 <- mcmc(m_2, ...)
  t2_2 <- Sys.time()
  
  #cat(paste0("Model 2: ", t2 - t1))
  
  # save the model
  saveRDS(d_2, paste0(path,species,"/fit_model2.rds"))
  cat(paste0("Fitting time for ", species,
             ": ", round(as.numeric(difftime(t2_2, t2_1, units = "min")), 3),
             " min"))
  # compute rhat
  tmp <- as.data.frame(coda::gelman.diag(d_2)$psrf)
  tmp$species <- species
  names(tmp)[1:2] <- c("point", "upper")
  tmp$param <- rownames(tmp)
  tmp$ESS <- coda::effectiveSize(d_2)
  rownames(tmp) <- NULL
  write.table(tmp, paste0(path, species, "/rhat_model2.csv"),
              sep = ",", row.names = FALSE)
  
  cat(paste0("Maximum rhat value is: ", tmp[which.max(tmp$point),"point"], "\n"))
  
  # compute the r2 and rmse
  r2 <- greta.checks::bayes_R2(y_mhb, mu_mhb, d_2)
  error <- greta.checks::rmse(y_mhb, mu_mhb, d_2, norm = TRUE)
  
  # put the two together
  r2_error <- as.data.frame(rbind(r2, error))
  r2_error$metric <- c("r2", "rmse")
  
  write.table(r2_error, paste0(path, species, "/rmse_model2.csv"),
              sep = ",", row.names = FALSE)
  
  # get the totals over Mhb plots
  newdat <- expand.grid(fJahr = factor(2005:2019, levels = 2005:2019),
                        route_id = unique(route_id))
  # predicted model matrix
  modmat_pred <- model.matrix(~ fJahr, 
                              newdat,
                              contrasts.arg = list(fJahr = "contr.sum"))
  mu_pred <- exp(modmat_pred[,-1] %*% yr_eff + route_avg[newdat$route_id])
  # calculate expectations
  d_pred <- t(as.matrix(calculate(mu_pred, values = d_2)))
  # some data wraggling
  d_df <- as.data.frame(d_pred)
  d_df$Jahr <- newdat$fJahr
  
  tt <- gather(d_df, key = "iteration", value = "val", -Jahr)
  
  tt %>%
    group_by(Jahr, iteration) %>%
    summarise(S = sum(val)) %>%
    ungroup() %>%
    group_by(Jahr) %>%
    summarise(M = mean(S), SD = sd(S),
              lo = quantile(S, probs = 0.025),
              hi = quantile(S, probs = 0.975)) -> tt_d
  
  write.table(tt_d, paste0(path, species, "/totals_model2.csv"),
              sep = ",", row.names = FALSE)
  
}


## function to fit model 4
fit_model4 <- function(species, path = "", ...){
  
  # filter Mhb data
  dat_mhb <-  filter(mhb, Artname == species)
  
  
  ## ornitho list
  dat_ornitho_list <- subset(ornitho_list, SPECIES_NAME == species)

  # ornitho unstr
  dat_ornitho_unstr <- subset(ornitho_unstr, SPECIES_NAME == species)
  # std the effort
  dat_ornitho_unstr$effort2 <- scale(dat_ornitho_unstr$effort2)
  
  
  
  dat_mhb$fJahr2 <- factor(dat_mhb$Jahr, levels = 2005:2019)
  dat_ornitho_list$fyear <- factor(dat_ornitho_list$year,
                                   levels = 2005:2019)
  
  dat_ornitho_unstr$fyear <- factor(dat_ornitho_unstr$year,
                                    levels = 2005:2019)
  
  # the common year effect
  yr_eff_mhb <- model.matrix(~fJahr2, dat_mhb,
                             contrasts.arg = list(fJahr2 = "contr.sum"))
  yr_eff_list <- model.matrix(~fyear, dat_ornitho_list,
                              contrasts.arg = list(fyear = "contr.sum"))
  yr_eff_unstr <- model.matrix(~fyear, dat_ornitho_unstr,
                               contrasts.arg = list(fyear = "contr.sum"))
  
  # prior for common year effect
  yr_eff <- normal(0, 1, dim = ncol(yr_eff_mhb) - 1)
  
  ## model for mhb
  # index for routes
  route_id <- as.numeric(factor(dat_mhb$ROUTENCODE))
  
  # priors
  intercept_mhb <- student(3, 0, 10)
  sd_route <- exp(student(3, 0, 10))
  alpha_route <- normal(0, 1, dim = max(route_id))
  route_avg <- intercept_mhb + sd_route * alpha_route
  sd_res_mhb <- gamma(0.01, 0.01)
  
  # linear predictor
  mu_mhb <- exp(yr_eff_mhb[,-1] %*% yr_eff + route_avg[route_id])
  
  # re-param for negative binomial
  prob_mhb <- sd_res_mhb / (mu_mhb + sd_res_mhb)
  
  
  ## model for ornitho lists
  # index for locality
  locality_id <- as.numeric(factor(dat_ornitho_list$ID_PLACE))
  # prior fixed effect
  intercept_list <- student(3, 0, 10)
  #dens_eff_list <- normal(0, 1)
  
  # prior for random effect
  sd_locality <- exp(student(3, 0, 10))
  alpha_locality <- normal(0, 1, dim = max(locality_id))
  locality_avg <- intercept_list + alpha_locality * sd_locality
  #zero inflation
  theta <- beta(1, 1)
  
  
  mu_list <- dat_ornitho_list$effort_h * exp(yr_eff_list[,-1] %*% yr_eff +
                                               locality_avg[locality_id])
  

  ## model for ornitho unstr
  # index for raster cell
  raster_id <- as.numeric(factor(dat_ornitho_unstr$ID_PLACE))
  # prior fixed effect
  intercept_unstr <- greta::student(3, 0, 10)
  slp <- normal(0, 1)
  #dens_eff_unstr <- normal(0, 1)
  
  # prior for random effect
  sd_raster <- exp(greta::student(3, 0, 10))
  alpha_raster <- normal(0, 1, dim = max(raster_id))
  raster_avg <- intercept_unstr + alpha_raster * sd_raster
  #sd_res_unstr <- gamma(0.01, 0.01)
  
  mu_unstr <- exp(slp * dat_ornitho_unstr$effort2 +
                    yr_eff_unstr[,-1] %*% yr_eff + raster_avg[raster_id])

  # the response
  y_mhb <- as_data(dat_mhb$count)
  y_list <- as_data(dat_ornitho_list$count)
  y_unstr <- as_data(dat_ornitho_unstr$count)
  
  distribution(y_mhb) <- negative_binomial(prob = prob_mhb, size = sd_res_mhb)
  distribution(y_list) <- zero_inflated_poisson(theta = theta, lambda = mu_list, dim = length(y_list))
  distribution(y_unstr) <- poisson(mu_unstr)
  
  
  m_3 <- model(intercept_mhb, intercept_list, intercept_unstr, slp,
               yr_eff, sd_res_mhb, theta, 
               sd_route, sd_locality, sd_raster)
  
  cat("Starting fitting model 3")
  t3_1 <- Sys.time()
  d_3 <- mcmc(m_3, ...)
  t3_2 <- Sys.time()
  
  #cat(paste0("Model 3: ", t2 - t1))
  
  # save the model
  saveRDS(d_3, paste0(path,species,"/fit_model3.rds"))
  cat(paste0("Fitting time for ", species,
             ": ", round(as.numeric(difftime(t3_2, t3_1, units = "min")), 3),
             " min"))
  
  # compute rhat
  tmp <- as.data.frame(coda::gelman.diag(d_3)$psrf)
  tmp$species <- species
  names(tmp)[1:2] <- c("point", "upper")
  tmp$param <- rownames(tmp)
  rownames(tmp) <- NULL
  tmp$ESS <- coda::effectiveSize(d_3)
  write.table(tmp, paste0(path,  species, "/rhat_model3.csv"),
              sep = ",", row.names = FALSE)
  
  cat(paste0("Maximum rhat value is: ", tmp[which.max(tmp$point),"point"], "\n"))
  
  # do some standard pp_check and bayes_R2 computation
  r2 <- greta.checks::bayes_R2(y_mhb, mu_mhb, d_3)
  error <- greta.checks::rmse(y_mhb, mu_mhb, d_3, norm = TRUE)
  # put the two together
  r2_error <- as.data.frame(rbind(r2, error))
  r2_error$metric <- c("r2", "rmse")
  write.table(r2_error, paste0(path,  species, "/rmse_model3.csv"),
              sep = ",", row.names = FALSE)
  
  ## compute breeding totals
  newdat <- expand.grid(fJahr = factor(2005:2019, levels = 2005:2019),
                        route_id = unique(route_id))
  # predicted model matrix
  modmat_pred <- model.matrix(~ fJahr, 
                              newdat,
                              contrasts.arg = list(fJahr = "contr.sum"))
  mu_pred <- exp(modmat_pred[,-1] %*% yr_eff + route_avg[newdat$route_id])
  # calculate expectations
  d_pred <- t(as.matrix(calculate(mu_pred, values = d_3)))
  # some data wraggling
  d_df <- as.data.frame(d_pred)
  d_df$Jahr <- newdat$fJahr
  
  tt <- gather(d_df, key = "iteration", value = "val", -Jahr)
  
  tt %>%
    group_by(Jahr, iteration) %>%
    summarise(S = sum(val)) %>%
    ungroup() %>%
    group_by(Jahr) %>%
    summarise(M = mean(S), SD = sd(S),
              lo = quantile(S, probs = 0.025),
              hi = quantile(S, probs = 0.975)) -> tt_d
  
  write.table(tt_d, paste0(path,  species, "/totals_model3.csv"),
              sep = ",", row.names = FALSE)

}

## function to fit model 5

fit_model5 <- function(species, path = "", ...){
  # filter Mhb data
  dat_mhb <-  filter(mhb, Artname == species)
  
  
  ## complete list (eBird and ornitho)
  dat_list_all <- subset(list_all, SPECIES_NAME == species)

  dat_mhb$fJahr2 <- factor(dat_mhb$Jahr, levels = 2005:2019)
  dat_list_all$fyear <- factor(dat_list_all$year,
                               levels = 2005:2019)
  
  # the common year effect
  yr_eff_mhb <- model.matrix(~fJahr2, dat_mhb,
                             contrasts.arg = list(fJahr2 = "contr.sum"))
  yr_eff_list <- model.matrix(~fyear, dat_list_all,
                              contrasts.arg = list(fyear = "contr.sum"))
  
  # prior for common year effect
  yr_eff <- normal(0, 1, dim = ncol(yr_eff_mhb) - 1)
  
  ## model for mhb
  # index for routes
  route_id <- as.numeric(factor(dat_mhb$ROUTENCODE))
  
  # priors
  intercept_mhb <- student(3, 0, 10)
  sd_route <- exp(student(3, 0, 10))
  alpha_route <- normal(0, 1, dim = max(route_id))
  route_avg <- intercept_mhb + sd_route * alpha_route
  sd_res_mhb <- gamma(0.01, 0.01)
  
  # linear predictor
  mu_mhb <- exp(yr_eff_mhb[,-1] %*% yr_eff + route_avg[route_id])
  
  # re-param for negative binomial
  prob_mhb <- sd_res_mhb / (mu_mhb + sd_res_mhb)
  
  
  ## model for ornitho lists
  # index for locality
  locality_id <- as.numeric(factor(dat_list_all$ID_PLACE))
  # prior fixed effect
  intercept_list <- student(3, 0, 10)
  #dens_eff_list <- normal(0, 1)
  
  # prior for random effect
  sd_locality <- exp(student(3, 0, 10))
  alpha_locality <- normal(0, 1, dim = max(locality_id))
  locality_avg <- intercept_list + sd_locality * alpha_locality
  
  # prior for zero-inflation
  theta <- beta(1, 1)
  
  
  mu_list <- dat_list_all$effort_h * exp(yr_eff_list[,-1] %*% yr_eff +
                                           locality_avg[locality_id])
  

  # the response
  y_mhb <- as_data(dat_mhb$count)
  y_list <- as_data(dat_list_all$count)
  
  
  distribution(y_mhb) <- negative_binomial(prob = prob_mhb, size = sd_res_mhb)
  distribution(y_list) <- greta.multivariate::zero_inflated_poisson(theta,
                                                                    mu_list,
                                                                    dim = length(y_list))
  
  
  m_4 <- model(intercept_mhb, intercept_list, 
               yr_eff, sd_res_mhb, theta, 
               sd_route, sd_locality)
  
  cat("Starting fitting model 4")
  t4_1 <- Sys.time()
  d_4 <- mcmc(m_4, ...)
  t4_2 <- Sys.time()
  
  #cat(paste0("Model 3: ", t2 - t1))
  
  # save the model
  saveRDS(d_4, paste0(path,species,"/fit_model4.rds"))
  cat(paste0("Fitting time for ", species,
             ": ", round(as.numeric(difftime(t4_2, t4_1, units = "min")), 3),
             " min\n"))
  # compute rhat
  tmp <- as.data.frame(coda::gelman.diag(d_4)$psrf)
  tmp$species <- species
  names(tmp)[1:2] <- c("point", "upper")
  tmp$param <- rownames(tmp)
  rownames(tmp) <- NULL
  tmp$ESS <- coda::effectiveSize(d_4)
  write.table(tmp, paste0(path, species, "/rhat_model4.csv"),
              sep = ",", row.names = FALSE)
  
  cat(paste0("Maximum rhat value is: ", tmp[which.max(tmp$point),"point"], "\n"))
  
  # do some standard pp_check and bayes_R2 computation
  r2 <- greta.checks::bayes_R2(y_mhb, mu_mhb, d_4)
  error <- greta.checks::rmse(y_mhb, mu_mhb, d_4, norm = TRUE)
  
  # put the two together
  r2_error <- as.data.frame(rbind(r2, error))
  r2_error$metric <- c("r2", "rmse")
  
  write.table(r2_error, paste0(path, species, "/rmse_model4.csv"),
              sep = ",", row.names = FALSE)
  

  ## compute breeding totals
  newdat <- expand.grid(fJahr = factor(2005:2019, levels = 2005:2019),
                        route_id = unique(route_id))
  # predicted model matrix
  modmat_pred <- model.matrix(~ fJahr, 
                              newdat,
                              contrasts.arg = list(fJahr = "contr.sum"))
  mu_pred <- exp(modmat_pred[,-1] %*% yr_eff + route_avg[newdat$route_id])
  # calculate expectations
  d_pred <- t(as.matrix(calculate(mu_pred, values = d_4)))
  # some data wraggling
  d_df <- as.data.frame(d_pred)
  d_df$Jahr <- newdat$fJahr
  
  tt <- gather(d_df, key = "iteration", value = "val", -Jahr)
  
  tt %>%
    group_by(Jahr, iteration) %>%
    summarise(S = sum(val)) %>%
    ungroup() %>%
    group_by(Jahr) %>%
    summarise(M = mean(S), SD = sd(S),
              lo = quantile(S, probs = 0.025),
              hi = quantile(S, probs = 0.975)) -> tt_d
  
  write.table(tt_d, paste0(path, species, "/totals_model4.csv"),
              sep = ",", row.names = FALSE)

}

## function to fit model 6
fit_model6 <- function(species, path = "", ...){
  
  # filter Mhb data
  dat_mhb <-  filter(mhb, Artname == species)
  
  ## ornitho list
  dat_list_all <- subset(list_all, SPECIES_NAME == species)

  # ornitho unstr
  dat_unstr_all <- subset(unstr_all, SPECIES_NAME == species)
  dat_unstr_all$effort2 <- scale(dat_unstr_all$effort2) # std effort
  
  
  
  dat_mhb$fJahr2 <- factor(dat_mhb$Jahr, levels = 2005:2019)
  dat_list_all$fyear <- factor(dat_list_all$year,
                               levels = 2005:2019)
  
  dat_unstr_all$fyear <- factor(dat_unstr_all$year,
                                levels = 2005:2019)
  
  # the common year effect
  yr_eff_mhb <- model.matrix(~fJahr2, dat_mhb,
                             contrasts.arg = list(fJahr2 = "contr.sum"))
  yr_eff_list <- model.matrix(~fyear, dat_list_all,
                              contrasts.arg = list(fyear = "contr.sum"))
  yr_eff_unstr <- model.matrix(~fyear, dat_unstr_all,
                               contrasts.arg = list(fyear = "contr.sum"))
  
  # prior for common year effect
  yr_eff <- normal(0, 1, dim = ncol(yr_eff_mhb) - 1)
  
  ## model for mhb
  # index for routes
  route_id <- as.numeric(factor(dat_mhb$ROUTENCODE))
  
  # priors
  intercept_mhb <- student(3, 1.6, 2.5)
  sd_route <- exp(student(3, 0, 2.5))
  alpha_route <- normal(0, 1, dim = max(route_id))
  route_avg <- intercept_mhb + sd_route * alpha_route
  sd_res_mhb <- gamma(0.01, 0.01)
  
  # linear predictor
  mu_mhb <- exp(yr_eff_mhb[,-1] %*% yr_eff + route_avg[route_id])
  
  # re-param for negative binomial
  prob_mhb <- sd_res_mhb / (mu_mhb + sd_res_mhb)
  
  
  ## model for ornitho lists
  # index for locality
  locality_id <- as.numeric(factor(dat_list_all$ID_PLACE))
  # prior fixed effect
  intercept_list <- student(3, 0, 3.4)
  #dens_eff_list <- normal(0, 1)
  
  # prior for random effect
  sd_locality <- exp(student(3, 0, 3.4))
  alpha_locality <- normal(0, 1, dim = max(locality_id))
  locality_avg <- intercept_list + alpha_locality * sd_locality
  # zero inflation
  theta <- beta(1, 1)
  
  
  mu_list <- dat_list_all$effort_h * exp(yr_eff_list[,-1] %*% yr_eff +
                                           locality_avg[locality_id])
  
  ## model for ornitho unstr
  # index for raster cell
  
  raster_id <- as.numeric(factor(dat_unstr_all$ID_PLACE))
  # prior fixed effect
  intercept_unstr <- student(3, 1.6, 2.5)
  slp <- normal(0, 1)

  # prior for random effect
  sd_raster <- exp(student(3, 0, 2.5))
  alpha_raster <- normal(0, 1, dim = max(raster_id))
  raster_avg <- intercept_unstr + alpha_raster * sd_raster

  mu_unstr <- exp(slp * dat_unstr_all$effort2 +
                    yr_eff_unstr[,-1] %*% yr_eff + raster_avg[raster_id])

  
  # the response
  y_mhb <- as_data(dat_mhb$count)
  y_list <- as_data(dat_list_all$count)
  y_unstr <- as_data(dat_unstr_all$count)
  
  distribution(y_mhb) <- negative_binomial(prob = prob_mhb, size = sd_res_mhb)
  distribution(y_list) <- zero_inflated_poisson(theta = theta, lambda = mu_list, dim = length(y_list))
  distribution(y_unstr) <- poisson(mu_unstr)
  
  m_5 <- model(intercept_mhb, intercept_list, intercept_unstr,
               yr_eff, sd_res_mhb, theta, slp,
               sd_route, sd_locality, sd_raster)
  
  
  
  
  cat("Starting fitting model 5")
  t5_1 <- Sys.time()
  d_5 <- mcmc(m_5, ...)
  t5_2 <- Sys.time()
  
  #cat(paste0("Model 3: ", t2 - t1))
  
  # save the model
  saveRDS(d_5, paste0(path,  species,"/fit_model5.rds"))
  cat(paste0("Fitting time for ", species,
             ": ", round(as.numeric(difftime(t5_2, t5_1, units = "min")), 3),
             " min"))
  
  # compute rhat
  tmp <- as.data.frame(coda::gelman.diag(d_5)$psrf)
  tmp$species <- species
  names(tmp)[1:2] <- c("point", "upper")
  tmp$param <- rownames(tmp)
  rownames(tmp) <- NULL
  tmp$ESS <- coda::effectiveSize(d_5)
  write.table(tmp, paste0(path,  species, "/rhat_model5.csv"),
              sep = ",", row.names = FALSE)
  
  cat(paste0("Maximum rhat value is: ", tmp[which.max(tmp$point),"point"], "\n"))
  
  # compute R2 and normalized RMSE
  r2 <- greta.checks::bayes_R2(y_mhb, mu_mhb, d_5)
  error <- greta.checks::rmse(y_mhb, mu_mhb, d_5, norm = TRUE)
  
  # put the two together
  r2_error <- as.data.frame(rbind(r2, error))
  r2_error$metric <- c("r2", "rmse")
  
  write.table(r2_error, paste0(path,  species, "/rmse_model5.csv"),
              sep = ",", row.names = FALSE)
  
  ## compute breeding totals
  newdat <- expand.grid(fJahr = factor(2005:2019, levels = 2005:2019),
                        route_id = unique(route_id))
  # predicted model matrix
  modmat_pred <- model.matrix(~ fJahr, 
                              newdat,
                              contrasts.arg = list(fJahr = "contr.sum"))
  mu_pred <- exp(modmat_pred[,-1] %*% yr_eff + route_avg[newdat$route_id])
  # calculate expectations
  d_pred <- t(as.matrix(calculate(mu_pred, values = d_5)))
  # some data wraggling
  d_df <- as.data.frame(d_pred)
  d_df$Jahr <- newdat$fJahr
  
  tt <- gather(d_df, key = "iteration", value = "val", -Jahr)
  
  tt %>%
    group_by(Jahr, iteration) %>%
    summarise(S = sum(val)) %>%
    ungroup() %>%
    group_by(Jahr) %>%
    summarise(M = mean(S), SD = sd(S),
              lo = quantile(S, probs = 0.025),
              hi = quantile(S, probs = 0.975)) -> tt_d
  
  write.table(tt_d, paste0(path,  species, "/totals_model5.csv"),
              sep = ",", row.names = FALSE)

}





