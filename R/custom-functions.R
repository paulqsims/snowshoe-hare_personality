################################################################################
# Function: Custom functions
# Author: Paul Q. Sims 
# Contact: paul.q.sims@gmail.com
# Date: 2020
# Purpose: Custom functions for analyses and formatting
################################################################################

# Reads in data and converts characters to factors
# Also mutates variables
#
# Only need to input file name without extension in quotes
# Assumes file is in the data directory and data are .csv files
# 
# read_mod_data <- function(data) {
#   read_csv(paste0("data/", data, ".csv", sep = "")) %>%
#     mutate(across(where(is.character), ~ as_factor(.x)),
#            trial = as.factor(trial),
#            across(.cols = c(goal_z_lat, body_length,  # log transformations
#                             learn_prop, tot_z),  
#                   ~ log(.x),
#                   .names = "{col}_LN"),
#            across(.cols = c(body_length, tot_z),  # mean center and scale variables by 1 SD
#                   ~ na_rm_scale(.x, na.rm = T),
#                   .names = "{col}_sc"))
#}

# Create dataset without NA values for variables of interest
#  - ONLY put variables that will be used for model

remove_dat_na <- function(data, variables) {
  # variables is a vector of strings of variable names you want to keep in the dataset
  data %>%
    select(all_of(variables)) %>%  # select which variables go into the data set
    drop_na()  # remove NA values
} 

# Mean centers and standardizes (1 SD) vector with NA value removal option
#
# Note that the default is False so removal must be specified
# scale() from base R does not allow this

na_rm_scale <- function(x, na.rm = FALSE) {
  (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
}

#### Rounding functions ####

# Rounds a vector of p-values
#
# Distinguishes between 0 and <0.001
# Output is character vector!

round_pval <- function(p_value) {
  p_value <- ifelse(p_value < 0.001, "<0.001",
                    ifelse(p_value > 0.001 & p_value <= 0.1,
                           round(p_value, 3),
                           round(p_value,2)))
  p_value <- as.character(p_value)  # change to character to make compatible w/ df that have P <0.001
}

# Rounds a continuous vector statistic
# 
# Distinguishes between 0 and < 0.01
# Output is character vector!

round_est <- function(estimate) {
  estimate <- ifelse(estimate < 0.01 & estimate > 0, "<0.01",
                     ifelse(estimate > -0.01 & estimate < 0,
                            "<-0.01",
                            round(estimate,2)))
  estimate <- as.character(estimate) # change to character to make compatible w/ df that have estimates < 0.01
}

# Round tidy model output  
#  - Rounds estimates to 2, p-values to 3
rd_tidy_out <- function(tidyOutput) {
  tidyOutput %>%
    mutate(across(.cols = c(estimate:statistic), ~ round_est(.x)),  # round est
           p.value = round_pval(p.value)) %>%  # round p-value
    return(.)
}

# Round stepwise selection output  
#  - Rounds estimates to 2, p-values to 3

rd_stepwise_out <- function(stepwiseOutput, lm = FALSE, 
                            glm = FALSE, glmmTMB = FALSE) {
  stepwiseOutput %>%
    rownames_to_column(., var = "Variable") %>%  # col of variable names
    as_tibble(.) %>%
    {if (lm == TRUE) {
      mutate(., across(.cols = npar:LRT, ~ round_est(.x)),  # round est
             p.value = round_pval(`Pr(Chi)`)) %>%  # round p-value
        select(-`Pr(Chi)`) %>%
        return(.)
    } else if (glmmTMB == TRUE) {
      mutate(., across(.cols = Df:LRT, ~ round_est(.x)),  # round est
             p.value = round_pval(`Pr(>Chi)`)) %>%  # round p-value
        select(-`Pr(>Chi)`) %>%
        return(.)
    } else if (glm == TRUE) {
    mutate(., across(.cols = npar:LRT, ~ round_est(.x)),  # round est
           p.value = round_pval(`Pr(Chi)`)) %>%  # round p-value
    select(-`Pr(Chi)`) %>%
    return(.)
    }}
}



#### Kable and Tidy functions ####

# Kable title and alignment
#  - Creates title and aligns table location on page
kable_title <- function(dataframe, title = NULL) {
  knitr::kable(dataframe, align = "l",  # align to the left
               caption = if (!is.null(title)) { 
                 paste(title)
               } else { paste("") }) %>%  # paste blank if no title specified
    return(.)
}

# Tidy predictor tables
#  - Tidies model output depending on mixed or non-mixed
#  - Creates title and output can be df or kable for pretty kable
pretty_PredictTab <- function(modelOutput, title = NULL,
                              mixedModel = FALSE, kable = TRUE) {
  # Arguments
  #  kable = whether or not a kable should be printed
  if (mixedModel == TRUE) {
    broom.mixed::tidy(modelOutput,
                      effects = "fixed") %>%
      rd_tidy_out(.) %>%  # round tidy output
      {if (kable == TRUE) {  # return kable
        kable_title(., title)  # modify kable title and alignment
      } else { return(.) }}  # return dataframe, not kable
  } else {  # if not a mixed model
    broom::tidy(modelOutput) %>%
      rd_tidy_out(.) %>%  # round tidy output
      {if (kable == TRUE) {
        kable_title(., title)  # modify kable title and alignment
      } else { return(.) }}  # return dataframe, not kable
  }
}


#### Repeatability functions ####

# Repeatability and variance estimates for behaviors measured once per day
#  - These are customized to the structure of each model as random effects and 
#    model family varied
# - delta method DETAILS HERE
#  ONLY WORK BECAUSE THEY HAVE THE SAME STRUCTURE

rpt_day <- function(model) {
  ind.int <- VarCorr(model)$ID[1]  # individual variance
  data_temp <- model.frame(model)  # extract data frame from model for refitting null mod
  m_null <- glmer(as.formula(paste(formula(model)[[2]], "~ 1 + (1|ID)")),  # intercept only null model
                  family = "binomial", 
                  data = data_temp,
                  glmerControl(optimizer="bobyqa",
                               optCtrl = list(maxfun = 1e6)))
  VarTotNull <- sum(VarCorr(m_null)$ID[1]) # total pheno var. from null model
  pmean <- plogis(as.numeric(fixef(m_null)) - 0.5 * VarTotNull *  # mean prob
                    tanh(as.numeric((fixef(m_null)) *
                                      (1 + 2 * exp(-0.5 * VarTotNull)))/6))
  resid.var <- 1/(pmean * (1 - pmean))  # residual variance
  long.term.rpt <- ind.int/(sum(as.numeric(VarCorr(model))) + resid.var)  # rpt
  
  value <- c(ind.int, resid.var, long.term.rpt)  # ind. var, resid. var, and rtp
  value  # put results in a dataframe
}


# Repeatability and variance estimates for Open Field behavior

varEst_OF <- function(x) {
  ind_int <- VarCorr(x)$ID[1]
  series_int <- VarCorr(x)$trial_series[1]
  series_slope <- VarCorr(x)$trial_series.1[1]
  # No intercept-slope corr b/c converg. issues in bootstrap
  #  series.int.slp.cov <- attr(VarCorr(x)$trial.series,"correlation")[2,1]
  dat_temp <- model.frame(x) # extract data frame from model for refitting null mod
  # Build null model
  #  - Intercept only fixed effect
  #  - Custom RE formula based on random effect model selection
  m_null <- glmer(of_loc ~ 1 + (1|ID) + 
                    (0 + t_point_sc|trial_series) + (1|trial_series), 
                  family = "binomial", 
                  data = dat_temp,
                  glmerControl(optimizer="bobyqa",
                               optCtrl = list(maxfun = 100000))) 
  VarF <- var(as.vector(model.matrix(x) %*% fixef(x))) # Calculation of the variance in fitted values from full model
  # Total pheno. var. from null model
  Vt <- sum(VarCorr(m_null)$ID[1] +              # Ind intercept
              VarCorr(m_null)$trial_series[1] +  # series intercept
              VarCorr(m_null)$trial_series.1[1])   # series slope
  pmean <- plogis(as.numeric(fixef(m_null)) - 0.5 * Vt * tanh(as.numeric((fixef(m_null)) * 
                                                                           (1 + 2 * exp(-0.5 * Vt)))/6))
  resid_var <- 1/(pmean * (1 - pmean))
  series_int_rpt <- ind_int/(ind_int + series_int)
  long_term_rpt <- ind_int/(ind_int + series_int + series_slope + resid_var)
  # individual intercept, series intercept, series slope, residual variance, 
  # reaction norm rpt, long term repeatability
  c(ind_int, series_int, series_slope, resid_var,
    series_int_rpt, long_term_rpt)
}

# Repeatability and variance estimates for shelter behavior

varEst_shel <- function(x) {
  #ind_int <- VarCorr(x)$ID[1]
  series_int <- VarCorr(x)$trial_series[1]
  series_slope <- VarCorr(x)$trial_series.1[1]
  #series_IntSlp_cov <- attr(VarCorr(x)$trial_series,"correlation")[2,1]
  dat_temp <- model.frame(x) # extract data frame from model for refitting null mod
  # Build null model
  #  - Intercept only fixed effect
  #  - Custom RE formula based on random effect model selection
  m_null <- glmer(shelter ~ 1 + (1|trial_series) + (0 + t_point_sc|trial_series),
                  family = "binomial",
                  data = dat_temp,
                  glmerControl(optimizer="bobyqa",
                               optCtrl = list(maxfun = 1e6)))
  VarF <- var(as.vector(model.matrix(x) %*% fixef(x))) # Calculation of the variance in fitted values from full model
  # Total pheno. var. from null model
  Vt <- sum(VarCorr(m_null)$trial_series[1] +  # series intercept
              VarCorr(m_null)$trial_series.1[1])  # series slope
  #attr(VarCorr(m_null)$trial_series,"correlation")[2,1] + # series intercept-slope correlation
  #VarCorr(m_null)$obsID[1])  # obs ID
  pmean <- plogis(as.numeric(fixef(m_null)) - 0.5 * Vt * tanh(as.numeric((fixef(m_null)) *
                                                                           (1 + 2 * exp(-0.5 * Vt)))/6))
  resid_var <- 1/(pmean * (1 - pmean))
  # series_rpt <- ind_int/(ind_int + series_int)
  # long_term_rpt <- ind_int/(ind_int + series_int + series_slope + resid_var)
  # individual intercept, series intercept, series slope, residual variance, 
  # reaction norm rpt, long term repeatability
  c(series_int, series_slope, resid_var)
}

# Repeatability and variance estimates for background preferences

varEst_bkgrd <- function(x) {
  ind_int <- VarCorr(x)$ID[1]
  ind_slope <- VarCorr(x)$ID.1[1]
  series_int <- VarCorr(x)$trial_series[1]
  series_slope <- VarCorr(x)$trial_series.1[1]
  # series_IntSlp_cov <- attr(VarCorr(x)$trial_series,"correlation")[2,1]
  dat_temp <- model.frame(x) # extract data frame from model for refitting null mod
  # Build null model
  #  - Intercept only fixed effect
  #  - Custom RE formula based on random effect model selection
  m_null <- glmer(bkgrd ~ 1 + (1|ID) + (0 + t_point_sc|ID) +
                    (1|trial_series) + (0 + t_point_sc|trial_series),
                  family = "binomial",
                  data = dat_temp,
                  glmerControl(optimizer="bobyqa",
                               optCtrl = list(maxfun = 100000)))
  VarF <- var(as.vector(model.matrix(x) %*% fixef(x))) # Calculation of the variance in fitted values from full model
  # Total pheno. var. from null model
  Vt <- sum(VarCorr(m_null)$ID[1] +              # Ind intercept
            VarCorr(m_null)$ID.1[1] +  # Ind slope
              VarCorr(m_null)$trial_series[1] +  # series intercept
              VarCorr(m_null)$trial_series.1[1])   # series slope
  pmean <- plogis(as.numeric(fixef(m_null)) - 0.5 * Vt * tanh(as.numeric((fixef(m_null)) *
                                                                           (1 + 2 * exp(-0.5 * Vt)))/6))
  resid_var <- 1/(pmean * (1 - pmean))
  series_int_rpt <- ind_int/(ind_int + series_int)
  series_slope_rpt <- ind_slope/(ind_slope + series_slope)
  long_term_rpt <- ind_int/(ind_int + ind_slope + series_int + series_slope + resid_var)
  # individual intercept, series intercept, series slope, residual variance,
  # reaction norm rpt, long term repeatability
  c(ind_int, ind_slope, series_int, series_slope, resid_var,
    series_int_rpt, series_slope_rpt, long_term_rpt)
}

# Repeatability and variance estimates for activity

varEst_act <- function(x) {
  ind_int <- VarCorr(x)$cond$ID[1]
  series_int <- VarCorr(x)$cond$trial_series.1[1]
  series_slope <- VarCorr(x)$cond$trial_series[1]
  #series_IntSlp_cov <- attr(VarCorr(x)$cond$trial_series,"correlation")[2,1]
  dat_temp <- model.frame(x) # extract data frame from model for refitting null mod
  # Build null model
  #  - Intercept only fixed effect
  #  - Custom RE formula based on random effect model selection
  m_null <- 
    glmmTMB(line_cross ~ 1 +
              (1|ID) + (0 + t_point_sc|trial_series) + (1|trial_series),   # CHECK TO SEE IF NEED TO CHANGE RE FOR CONVERG FROM BOOT
            data = dat_temp,
            ziformula = ~ 1,
            family ="nbinom2")
  VarF <- var(as.vector(model.matrix(x) %*% fixef(x)$cond)) # Calculation of the variance in fitted values from full model
  # Total pheno. var. from null model
  lambda <- as.numeric(exp(fixef(m_null)$cond + 0.5 *  # total var, see eqtn 5.8 in Nakagawa et al. 2017
                             (as.numeric(VarCorr(m_null)$cond$ID[1]) +
                                as.numeric(VarCorr(m_null)$cond$trial_series[1]) +
                                as.numeric(VarCorr(m_null)$cond$trial_series.1[1]))))
  thetaF <- summary(x)$sigma  # dispersion parameter theta from full model
  resid_var <- trigamma((1/lambda + 1/thetaF)^(-1))
  series_int_rpt <- ind_int/(ind_int + series_int)
  long_term_rpt <- ind_int/(ind_int + series_int + series_slope +
                            resid_var)
  # individual intercept, series intercept, series slope, residual variance,
  # reaction norm rpt, long term repeatability
  c(ind_int, series_int, series_slope, resid_var, 
    series_int_rpt, long_term_rpt)
}




# PRETTY REPEATABILITY OUTPUT
## ALSO REQUIRES TIDYVERSE LOADED
pretty_DayRpt <- function(RptDay_output) {
  var.name <- c("Individual Intercept", "Residual", "Repeatability")
  data.frame(var.name, RptDay_output) %>%
    knitr::kable(col.names = c("Variable", "Estimate"),
                 digits = 2, align = "l",
                 caption = "Variance and repeatability estimates")
}

# Confidence interval calculations
#  - Only works for rpt with same model structure, see variable below
RptConfInt_day <- function(bootMerOutput) {
  # Convert bootstrap output into dataframe
  dat_boot <- as_tibble(bootMerOutput$t) %>%
    rename("Individual Intercept" = "V1",
           "Residual" = "V2",
           "Repeatability" = "V3")
  
  # Calculate C.I. 
  dat_boot %>%
    map_df(~ stats::quantile(.x,  # map 95% C.I. for each column
                             c((1 - 0.95)/2, 1 - (1 - 0.95)/2),
                             na.rm = TRUE)) %>%
    mutate(Estimate = bootMerOutput$t0,
           Variable = c("Individual Intercept", "Residual", "Repeatability")) %>%
    select(Variable, Estimate, Lower95CI = `2.5%`, Upper95CI = `97.5%`) %>%
    knitr::kable(., digits = 2, align = "l",
                 caption = "Variance and repeatability estimates ± 95% CI") 
}

# Likelihood ratio test for significance of random effect
# REQUIRES TIDYVERSE LOADED B?C OF PIPELINE

rpt_SigTest <- function(fullmod, redmod) {
  anova(fullmod, redmod, test = "LRT") %>%
    mutate(p_value = ifelse(!is.na(`Pr(>Chisq)`),  # divide p-value by 2 b/c testing on the boundary 
                            `Pr(>Chisq)`/2, `Pr(>Chisq)`),
           p_value = round_pval(p_value)) %>%  
    select(-`Pr(>Chisq)`) %>%
    knitr::kable(., digits = 2, align = "l",
                 caption = "Likelihood ratio test for significance of
                 random effect")
}


#### Random effects ####

# Refit model with new random effect structure but same fixed effects
refitModelRanEf <- function(oldModel, random = NULL) {
  if (is.null(random)) { return(warning("Missing random effects")) }
  if (!is.character(random)) { 
    return(warning("Random effects must be character string"))
  }
  # extract previous fixed effect structure
  old_formula <- nobars(formula(oldModel))  
  # combine old fixed effect with new random effect structure
  new_formula <- formula(paste(as.character(c(old_formula)), " + ",  
                               random))
  # Refit the model with new random effect structure
  new_model <- update(oldModel, new_formula)
  return(new_model)
} 