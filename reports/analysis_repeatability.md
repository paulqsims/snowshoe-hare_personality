
# Repeatability analyses

  - Author: Paul Q. Sims
  - Contact: <paul.q.sims@gmail.com>
  - Date: 2020
  - Purpose: Repeatability analyses and variance estimates for Lafferty
    et al. 2020

<!-- end list -->

``` r
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file(),
                     eval = TRUE, echo = TRUE, message = FALSE,
                     warning = FALSE, cache = TRUE)
knitr::opts_chunk$set(root.dir = rprojroot::find_rstudio_root_file(),
                     eval = TRUE, echo = TRUE, message = FALSE,
                     warning = FALSE, cache = TRUE)
```

## Setup

``` r
# Load libraries
library(tidyverse)  # for cleaning and modifying data
library(rptR)  # for repeatability analyses
library(lme4)
library(lmerTest)
library(DHARMa)
library(broom.mixed)  # for tidying model output
library(knitr)
library(glmmTMB)  # for activity analyses

# Load personalized functions
source("R/custom-functions.R")

# Read in data 
# Data collected once each day
data_day <-
  read_csv("data/data_day_LaffertyEtAl_2020.csv",
           col_types = list(exit_fast = col_factor(),
                            trial = col_factor(levels = c("1", "2")),
                            sex = col_factor(levels = c("Female", "Male")))) %>%
  mutate(obsID = 1:nrow(.))  # observation level random effect


# Data collected at time points throughout the day
data_time <-
  read_csv("data/data_time_LaffertyEtAl_2020.csv",
             col_types = list(trial = col_factor(levels = c("1", "2")),
                              sex = col_factor(levels = c("Female", "Male")))) %>%
  mutate(obsID = 1:nrow(.),  # observation level random effect
         across(where(is.character), ~ as_factor(.x)))  

# Data collected at time points throughout the day for lines crossed
data_act <-
  read_csv("data/data_activity_LaffertyEtAl_2020.csv",
             col_types = list(trial = col_factor(levels = c("1", "2")),
                              sex = col_factor(levels = c("Female", "Male")))) %>%
  mutate(obsID = 1:nrow(.),  # observation level random effect
         across(where(is.character), ~ as_factor(.x)))  
```

## Area explored

  - Proportion of unique tiles entered

**Create dataset without NAs**

``` r
data_expl_rpt <- 
  remove_dat_na(data_day, c("uniq_tiles", "tile_fails", "trial",
                            "sex", "ID", "day", "obsID")) %>%
  mutate(expl.y = cbind(uniq_tiles, tile_fails))
```

### Random effect selection

**Test significance of individual intercept (1|ID)**

``` r
# Fit model for adjusted repeatability
m_AdjRpt_expl <- 
  glmer(expl.y ~ trial + sex + (1|ID),  
        family = "binomial", 
        data = data_expl_rpt,
        contrasts = list(trial = c(-1,1), sex = c(-1,1)),
        glmerControl(optimizer="bobyqa",
        optCtrl = list(maxfun = 1e6)))

# Fit null model for testing significance of individual ID
m_AdjRpt_expl_null <- 
  glm(expl.y ~ trial + sex,  
        family = "binomial", 
        data = data_expl_rpt,
        contrasts = list(trial = c(-1,1), sex = c(-1,1)))

# Likelihood ratio test for significance of ind. random effect and rpt
#  Custom functions
rpt_SigTest(m_AdjRpt_expl, m_AdjRpt_expl_null)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 3    | 224.39 | 229.46 | \-109.20 | 218.39   | NA    | NA | NA       |
| 4    | 211.25 | 218.01 | \-101.63 | 203.25   | 15.14 | 1  | \<0.001  |

Likelihood ratio test for significance of random effect

### Variance and repeatability estimates

``` r
# Adjusted repeatability and variance estimates
# Bootstrap 95% C.I. 
res_boot <- bootMer(m_AdjRpt_expl,  # final rpt model
                    use.u = FALSE,  
                    FUN = rpt_day,  # custom rpt function
                    seed = 1,  # allows exact replication of results
                    nsim = 10,
                    type = "parametric",
                    parallel = "multicore",
                    .progress = "txt")
```

    ## ====================================================================================================================

``` r
# Calculate 95% CI
RptConfInt_day(res_boot)
```

| Variable             | Estimate | Lower95CI | Upper95CI |
| :------------------- | :------- | :-------- | :-------- |
| Individual Intercept | 0.17     | 0.02      | 0.16      |
| Residual             | 4.07     | 4.01      | 4.21      |
| Repeatability        | 0.04     | 0.01      | 0.04      |

Variance and repeatability estimates ± 95% CI

## Risk-taking

### Shelter exit latency

  - Exit immediately: yes or no

**Create data set with NA values removed**

``` r
data_shelExit_rpt <- 
  remove_dat_na(data_day, c("exit_fast", "trial",
                            "sex", "ID", "day",
                            "obsID"))
```

#### Random effect selection

**Test significance of individual intercept (1|ID)**

``` r
# Fit model for adjusted repeatability
m_AdjRpt_ExitLat <- 
  glmer(exit_fast ~ trial + sex + (1|ID),  
        family = "binomial", 
        data = data_shelExit_rpt,
        contrasts = list(trial = c(-1,1), sex = c(-1,1)),
        glmerControl(optimizer="bobyqa",
        optCtrl = list(maxfun = 1e6)))

# Fit null model for testing significance of individual ID
m_AdjRpt_ExitLat_null <- 
  glm(exit_fast ~ trial + sex,  
        family = "binomial", 
        data = data_shelExit_rpt,
        contrasts = list(trial = c(-1,1), sex = c(-1,1)))

# Likelihood ratio test for significance of ind. random effect and rpt
rpt_SigTest(m_AdjRpt_ExitLat, m_AdjRpt_ExitLat_null)
```

| npar | AIC   | BIC   | logLik  | deviance | Chisq | Df | p\_value |
| :--- | :---- | :---- | :------ | :------- | :---- | :- | :------- |
| 3    | 54.76 | 59.75 | \-24.38 | 48.76    | NA    | NA | NA       |
| 4    | 54.96 | 61.62 | \-23.48 | 46.96    | 1.79  | 1  | 0.09     |

Likelihood ratio test for significance of random effect

#### Variance and repeatability estimates

``` r
# Variance estimates, including Adjusted repeatability 
rpt_day(m_AdjRpt_ExitLat) %>%
  pretty_DayRpt(.)
```

| Variable             | Estimate |
| :------------------- | :------- |
| Individual Intercept | 2.59     |
| Residual             | 4.06     |
| Repeatability        | 0.39     |

Variance and repeatability estimates

  - 95% C.I. not computed since individual ID is non-significant

### Shelter preference

  - Inside or outside the shelter

**Create data set with NA values removed**

``` r
data_shel_rpt <- 
  remove_dat_na(data_time, c("shelter", "t_point_sc", "trial",
                             "sex", "ID", "trial_series", "day",
                             "obsID"))
```

#### Random-effect selection

**Test significance of individual intercept (1|ID)**

``` r
# Fit base model
m_shel_base <- 
  glmer(shelter ~ t_point_sc * trial + sex +
          (1|ID) + (t_point_sc|trial_series) + (1|obsID),
        data = data_shel_rpt,
        family = "binomial", 
        contrasts = list(trial = c(-1,1), sex = c(-1,1)),
        glmerControl(optimizer="bobyqa",
                    optCtrl = list(maxfun = 1e6)))

# Fit model without individual ID
m_shel_reduc <- refitModelRanEf(m_shel_base,
                                 "(t_point_sc|trial_series) + (1|obsID)")

# LRT for significance of individual ID
rpt_SigTest(m_shel_base, m_shel_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 9    | 401.87 | 445.04 | \-191.93 | 383.87   | NA    | NA | NA       |
| 10   | 403.86 | 451.83 | \-191.93 | 383.86   | 0.01  | 1  | 0.47     |

Likelihood ratio test for significance of random effect

**Test significance of individual intercept and time slope correlation
(t\_point\_sc|ID)**

  - Run before testing random slope beccause result of slope intercept
    will determine which slope model to use

<!-- end list -->

``` r
# Fit model with intercept-slope correlation
m_shel_full <- refitModelRanEf(m_shel_base,
                                 "(t_point_sc|ID) + (t_point_sc|trial_series) +
                                  (1|obsID)")  

# Fit model without intercept-slope m_shel_base
m_shel_reduc <- refitModelRanEf(m_shel_base,
                                 "(0 + t_point_sc|ID) + (1|ID) +
                                  (t_point_sc|trial_series) + (1|obsID)")

# LRT for significance of intercept-slope correlation
rpt_SigTest(m_shel_full, m_shel_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 11   | 405.73 | 458.49 | \-191.86 | 383.73   | NA    | NA | NA       |
| 12   | 407.73 | 465.29 | \-191.86 | 383.73   | 0     | 1  | 0.5      |

Likelihood ratio test for significance of random effect

``` r
# print non-significant intercept-slope correlation
kable(as_tibble_col(attr(VarCorr(m_shel_full)$ID,"correlation")[2,1],
              column_name = "Random intercept-slope variance estimate"),
      align = "l")
```

| Random intercept-slope variance estimate |
| :--------------------------------------- |
| NaN                                      |

**Test significance of time slope (t\_point\_sc)**

``` r
# Fit model with random slope
m_shel_full <- refitModelRanEf(m_shel_base,
                                 "(0 + t_point_sc|ID) + (1|ID) +
                                  (t_point_sc|trial_series) + (1|obsID)")  

# Fit model without random slope
m_shel_reduc <- refitModelRanEf(m_shel_base,
                                 "(1|ID) + (t_point_sc|trial_series) + (1|obsID)")

# LRT for significance of random slope 
rpt_SigTest(m_shel_full, m_shel_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 10   | 403.86 | 451.83 | \-191.93 | 383.86   | NA    | NA | NA       |
| 11   | 405.73 | 458.49 | \-191.86 | 383.73   | 0.13  | 1  | 0.36     |

Likelihood ratio test for significance of random effect

**Test significance of (reaction norm) trial series intercept**

``` r
# Fit model with random intercept
m_shel_full <- refitModelRanEf(m_shel_base,
                                 "(1|ID) + (1|trial_series)")

# Fit model without random intercept
m_shel_reduc <- refitModelRanEf(m_shel_base,
                                 "(1|ID)")

# LRT for significance of intercept
rpt_SigTest(m_shel_full, m_shel_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 6    | 573.24 | 602.02 | \-280.62 | 561.24   | NA    | NA | NA       |
| 7    | 494.46 | 528.03 | \-240.23 | 480.46   | 80.79 | 1  | \<0.001  |

Likelihood ratio test for significance of random effect

**Test significance of (reaction norm) trial series intercept and time
slope correlation**

  - Run before testing random slope because result of slope intercept
    will determine which slope model to use

<!-- end list -->

``` r
# Fit model with random intercept-slope correlation
m_shel_full <- refitModelRanEf(m_shel_base,
                                 "(1|ID) + (t_point_sc|trial_series) + (1|obsID)")  

# Fit model without random intercept-slope correlation
m_shel_reduc <- refitModelRanEf(m_shel_base,
                                 "(1|ID) + (0 + t_point_sc|trial_series) +
                                  (1|trial_series) + (1|obsID)")

# LRT for significance of intercept-slope correlation
rpt_SigTest(m_shel_full, m_shel_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 9    | 414.80 | 457.97 | \-198.40 | 396.80   | NA    | NA | NA       |
| 10   | 403.86 | 451.83 | \-191.93 | 383.86   | 12.93 | 1  | \<0.001  |

Likelihood ratio test for significance of random effect

**Test significance of (reaction norm) time slope**

``` r
# Fit model with random slope
m_shel_full <- refitModelRanEf(m_shel_base,
                                 "(1|ID) + (t_point_sc|trial_series) + (1|obsID)")  

# Fit model without random slope
m_shel_reduc <- refitModelRanEf(m_shel_base,
                                 "(1|ID) + (1|trial_series) + (1|obsID)")

# LRT for significance of random slope
rpt_SigTest(m_shel_full, m_shel_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 8    | 496.46 | 534.83 | \-240.23 | 480.46   | NA    | NA | NA       |
| 10   | 403.86 | 451.83 | \-191.93 | 383.86   | 96.6  | 2  | \<0.001  |

Likelihood ratio test for significance of random effect

**Final random effect model structure**

  - Removed trial series intercept and time slope correlation due to
    model convergence errors during bootstrapping

<!-- end list -->

``` r
m_shel_rpt_final <- 
  glmer(shelter ~ t_point_sc * trial + sex +
          (0 + t_point_sc|trial_series) + (1|trial_series),  
        data = data_shel_rpt,
        family = "binomial", 
        contrasts = list(trial = c(-1,1), sex = c(-1,1)),
        glmerControl(optimizer="bobyqa",
                    optCtrl = list(maxfun = 1e6)))

kable(as_tibble_col(as.character(c(formula(m_shel_rpt_final))),
              column_name = "Final formula"))
```

| Final formula                                                                                     |
| :------------------------------------------------------------------------------------------------ |
| shelter \~ t\_point\_sc \* trial + sex + (0 + t\_point\_sc | trial\_series) + (1 | trial\_series) |

#### Variance estimates

  - Repeatability not calculated since individual intercept and time
    slope estimates were non-significant
  - Removed random series intercept-slope covariance due to convergence
    errors

<!-- end list -->

``` r
# Bootstrap confidence intervals
res_bootOut <-
  bootMer(m_shel_rpt_final,  # final rpt model
          use.u = FALSE,
          FUN = varEst_shel,  # custom variance estimates function
          seed = 1,
          nsim = 10,
          type = "parametric",
          parallel = "multicore",
          .progress = "txt")
```

    ## ====================================================================================================================

``` r
# Compile variance estimates 
res_varEst <-
  as_tibble(res_bootOut$t) %>%
  rename("Series Intercept" = "V1",
         "Series Slope" = "V2",
         "Residual" = "V3") %>%
  map_df(~ quantile(.x, c((1 - 0.95)/2, 1 - (1 - 0.95)/2),
                  na.rm = TRUE), .id = "scale") %>%
  rename(Variable = scale, Lower95CI = `2.5%`, Upper95CI = `97.5%`) %>%
  mutate(Estimate = pluck(res_bootOut$t0)) %>%  # add R estimates
  relocate(Variable, Estimate) %>%  # reorder columns
  mutate_at(c(2:4), round, 2)  # # round estimates and errors to 2 digits

kable(res_varEst)
```

| Variable         |  Estimate | Lower95CI |   Upper95CI |
| :--------------- | --------: | --------: | ----------: |
| Series Intercept |      6.40 |      3.89 | 1.69700e+01 |
| Series Slope     |     17.24 |      9.39 | 3.05300e+01 |
| Residual         | 109445.31 |     14.50 | 1.05822e+11 |

### Open-Field location preference

  - Inner vs. outer area of the Open-Field

**Create data set with NA values removed**

``` r
data_of_rpt <- 
  remove_dat_na(data_time, c("of_loc", "t_point_sc", "trial",
                             "sex", "ID", "trial_series", "day",
                             "obsID"))
```

#### Random-effect selection

**Test significance of individual intercept (1|ID)**

``` r
# Fit base model
m_OFloc_base <- 
  glmer(of_loc ~ t_point_sc * trial + sex +
          (1|ID) + (t_point_sc|trial_series),  
        data = data_of_rpt,
        family = "binomial", 
        contrasts = list(trial = c(-1,1), sex = c(-1,1)),
        glmerControl(optimizer="bobyqa",
                    optCtrl = list(maxfun = 1e6)))

# Fit model without individual ID
m_OFloc_reduc <- refitModelRanEf(m_OFloc_base,
                                 "(t_point_sc|trial_series)")

# LRT for significance of individual ID
rpt_SigTest(m_OFloc_base, m_OFloc_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 8    | 905.14 | 943.51 | \-444.57 | 889.14   | NA    | NA | NA       |
| 9    | 902.03 | 945.20 | \-442.01 | 884.03   | 5.11  | 1  | 0.012    |

Likelihood ratio test for significance of random effect

**Test significance of individual intercept and time slope correlation
(t\_point\_sc|ID)**

  - Run before testing random slope because result of slope intercept
    will determine which slope model to use

<!-- end list -->

``` r
# Fit model with intercept-slope correlation
m_OFloc_full <- refitModelRanEf(m_OFloc_base,
                                 "(t_point_sc|ID) + (t_point_sc|trial_series)")  

# Fit model without intercept-slope m_OFloc_base
m_OFloc_reduc <- refitModelRanEf(m_OFloc_base,
                                 "(0 + t_point_sc|ID) + (1|ID) +
                                  (t_point_sc|trial_series)")

# LRT for significance of intercept-slope correlation
rpt_SigTest(m_OFloc_full, m_OFloc_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 10   | 904.03 | 952.00 | \-442.01 | 884.03   | NA    | NA | NA       |
| 11   | 904.53 | 957.29 | \-441.26 | 882.53   | 1.5   | 1  | 0.11     |

Likelihood ratio test for significance of random effect

**Test significance of time slope (t\_point\_sc)**

``` r
# Fit model with random slope
m_OFloc_full <- refitModelRanEf(m_OFloc_base,
                                 "(0 + t_point_sc|ID) + (1|ID) +
                                  (t_point_sc|trial_series)") 

# Fit model without random slope
m_OFloc_reduc <- refitModelRanEf(m_OFloc_base,
                                 "(1|ID) + (t_point_sc|trial_series)")

# LRT for significance of random slope 
rpt_SigTest(m_OFloc_full, m_OFloc_reduc)
```

| npar | AIC    | BIC   | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :---- | :------- | :------- | :---- | :- | :------- |
| 9    | 902.03 | 945.2 | \-442.01 | 884.03   | NA    | NA | NA       |
| 10   | 904.03 | 952.0 | \-442.01 | 884.03   | 0     | 1  | 0.5      |

Likelihood ratio test for significance of random effect

``` r
# Non-significant random slope variance estimate
kable(as_tibble_col(VarCorr(m_OFloc_full)$ID.1[1],
              column_name = "Random slope variance estimate"),
      align = "l")
```

| Random slope variance estimate |
| :----------------------------- |
| 0                              |

**Test significance of (reaction norm) trial series intercept**

``` r
# Fit model with random intercept
m_OFloc_full <- refitModelRanEf(m_OFloc_base,
                                 "(1|ID) + (1|trial_series)")  

# Fit model without random intercept
m_OFloc_reduc <- refitModelRanEf(m_OFloc_base,
                                 "(1|ID)")

# LRT for significance of intercept
rpt_SigTest(m_OFloc_full, m_OFloc_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 6    | 948.44 | 977.22 | \-468.22 | 936.44   | NA    | NA | NA       |
| 7    | 924.01 | 957.59 | \-455.01 | 910.01   | 26.43 | 1  | \<0.001  |

Likelihood ratio test for significance of random effect

**Test significance of (reaction norm) trial series intercept and time
slope correlation**

  - Run before testing random slope b/c result of slope intercept will
    determine which slope model to use

<!-- end list -->

``` r
# Fit model with random intercept-slope correlation
m_OFloc_full <- refitModelRanEf(m_OFloc_base,
                                 "(1|ID) + (t_point_sc|trial_series)")  

# Fit model without random intercept-slope correlation
m_OFloc_reduc <- refitModelRanEf(m_OFloc_base,
                                 "(1|ID) + (0 + t_point_sc|trial_series) +
                                  (1|trial_series)")

# LRT for significance of intercept-slope correlation
rpt_SigTest(m_OFloc_full, m_OFloc_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 8    | 907.46 | 945.84 | \-445.73 | 891.46   | NA    | NA | NA       |
| 9    | 902.03 | 945.20 | \-442.01 | 884.03   | 7.43  | 1  | 0.003    |

Likelihood ratio test for significance of random effect

**Test significance of (reaction norm) time slope**

``` r
# Fit model with random slope
m_OFloc_full <- refitModelRanEf(m_OFloc_base,
                                 "(1|ID) + (t_point_sc|trial_series)")  

# Fit model without random slope
m_OFloc_reduc <- refitModelRanEf(m_OFloc_base,
                                 "(1|ID) + (1|trial_series)")

# LRT for significance of random slope
rpt_SigTest(m_OFloc_full, m_OFloc_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 7    | 924.01 | 957.59 | \-455.01 | 910.01   | NA    | NA | NA       |
| 9    | 902.03 | 945.20 | \-442.01 | 884.03   | 25.98 | 2  | \<0.001  |

Likelihood ratio test for significance of random effect

**Final random effect model structure**

  - Removed trial series intercept and time slope correlation due to
    model convergence errors during bootstrapping

<!-- end list -->

``` r
m_OFloc_rpt_final <- 
  glmer(of_loc ~ t_point_sc * trial + sex +
          (1|ID) + (0 + t_point_sc|trial_series) + (1|trial_series),  
        data = data_of_rpt,
        family = "binomial", 
        contrasts = list(trial = c(-1,1), sex = c(-1,1)),
        glmerControl(optimizer="bobyqa",
                    optCtrl = list(maxfun = 1e6)))

kable(as_tibble_col(as.character(c(formula(m_OFloc_rpt_final))),
              column_name = "Final formula"))
```

| Final formula                                                                                                |
| :----------------------------------------------------------------------------------------------------------- |
| of\_loc \~ t\_point\_sc \* trial + sex + (1 | ID) + (0 + t\_point\_sc | trial\_series) + (1 | trial\_series) |

#### Variance and repeatability estimates

``` r
# Bootstrap confidence intervals
res_bootOut <- 
  bootMer(m_OFloc_rpt_final,  # final rpt model 
          use.u = FALSE, 
          FUN = varEst_OF,  # custom variance estimates function
          seed = 1,
          nsim = 10,  
          type = "parametric",
          parallel = "multicore",
          .progress = "txt")
```

    ## ====================================================================================================================

``` r
# Compile variance and repeatability estimates 
res_varEst <- 
  as_tibble(res_bootOut$t) %>%
  rename("Individual Intercept" = "V1",
         "Series Intercept" = "V2",
         "Series Slope" = "V3",
         "Residual" = "V4",
         "Reaction Norm Intercept Repeatability" = "V5",
         "Long-term Repeatability" = "V6") %>%
  map_df(~ quantile(.x, c((1 - 0.95)/2, 1 - (1 - 0.95)/2),  # calculate 95CI from bootstrap
                  na.rm = TRUE), .id = "scale") %>%
  rename(Variable = scale, Lower95CI = `2.5%`, Upper95CI = `97.5%`) %>%
  mutate(Estimate = pluck(res_bootOut$t0)) %>%  # add R estimates
  relocate(Variable, Estimate) %>%  # reorder columns
  mutate_at(c(2:4), round, 2)  # # round estimates and errors to 2 digits

kable(res_varEst)
```

| Variable                              | Estimate | Lower95CI | Upper95CI |
| :------------------------------------ | -------: | --------: | --------: |
| Individual Intercept                  |     1.24 |      0.07 |      2.26 |
| Series Intercept                      |     0.86 |      0.34 |      1.35 |
| Series Slope                          |     0.43 |      0.22 |      0.67 |
| Residual                              |     4.75 |      4.54 |      5.23 |
| Reaction Norm Intercept Repeatability |     0.59 |      0.07 |      0.86 |
| Long-term Repeatability               |     0.17 |      0.01 |      0.30 |

## Background preference

  - Dark background: yes or no

**Create data set with NA values removed**

``` r
data_bkgrd_rpt <- 
  remove_dat_na(data_time, c("bkgrd", "t_point_sc", "trial",
                             "sex", "ID", "trial_series", "day",
                             "obsID"))
```

### Random-effect selection

**Test significance of individual intercept (1|ID)**

``` r
# Fit base model
m_bkgrd_base <- 
  glmer(bkgrd ~ t_point_sc * trial + sex +
          (1|ID) + (t_point_sc|trial_series), 
        data = data_bkgrd_rpt,
        family = "binomial", 
        contrasts = list(trial = c(-1,1), sex = c(-1,1)),
        glmerControl(optimizer="bobyqa",
                    optCtrl = list(maxfun = 1e6)))

# Fit model without individual ID
m_bkgrd_reduc <- refitModelRanEf(m_bkgrd_base,
                                 "(t_point_sc|trial_series)")

# LRT for significance bkgrd individual ID
rpt_SigTest(m_bkgrd_base, m_bkgrd_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 8    | 713.64 | 752.02 | \-348.82 | 697.64   | NA    | NA | NA       |
| 9    | 705.56 | 748.74 | \-343.78 | 687.56   | 10.08 | 1  | \<0.001  |

Likelihood ratio test for significance of random effect

**Test significance of individual intercept and time slope correlation
(t\_point\_sc|ID)**

  - Run before testing random slope because result background slope
    intercept will determine which slope model to use

<!-- end list -->

``` r
# Fit model with intercept-slope correlation
m_bkgrd_full <- refitModelRanEf(m_bkgrd_base,
                                 "(t_point_sc|ID) + (t_point_sc|trial_series)")  

# Fit model without intercept-slope m_bkgrd_base
m_bkgrd_reduc <- refitModelRanEf(m_bkgrd_base,
                                 "(0 + t_point_sc|ID) + (1|ID) +
                                  (t_point_sc|trial_series)")

# LRT for significance bkgrd intercept-slope correlation
rpt_SigTest(m_bkgrd_full, m_bkgrd_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 10   | 707.56 | 755.53 | \-343.78 | 687.56   | NA    | NA | NA       |
| 11   | 704.22 | 756.98 | \-341.11 | 682.22   | 5.35  | 1  | 0.01     |

Likelihood ratio test for significance of random effect

``` r
# print non-significant intercept-slope correlation
kable(as_tibble_col(attr(VarCorr(m_bkgrd_full)$ID,"correlation")[2,1],
              column_name = "Random intercept-slope variance estimate"),
      align = "l")
```

| Random intercept-slope variance estimate |
| :--------------------------------------- |
| 0.8791246                                |

**Test significance of time slope (t\_point\_sc)**

``` r
# Fit model with random slope
m_bkgrd_full <- refitModelRanEf(m_bkgrd_base,
                                 "(t_point_sc|ID) +
                                  (t_point_sc|trial_series)") 

# Fit model without random slope
m_bkgrd_reduc <- refitModelRanEf(m_bkgrd_base,
                                 "(1|ID) + (t_point_sc|trial_series)")

# LRT for significance bkgrd random slope 
rpt_SigTest(m_bkgrd_full, m_bkgrd_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 9    | 705.56 | 748.74 | \-343.78 | 687.56   | NA    | NA | NA       |
| 11   | 704.22 | 756.98 | \-341.11 | 682.22   | 5.35  | 2  | 0.035    |

Likelihood ratio test for significance of random effect

**Test significance of (reaction norm) trial series intercept**

``` r
# Fit model with random intercept
m_bkgrd_full <- refitModelRanEf(m_bkgrd_base,
                                 "(t_point_sc|ID) + (1|trial_series)")

# Fit model without random intercept
m_bkgrd_reduc <- refitModelRanEf(m_bkgrd_base,
                                 "(t_point_sc|ID)")

# LRT for significance bkgrd intercept
rpt_SigTest(m_bkgrd_full, m_bkgrd_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 8    | 754.25 | 792.62 | \-369.12 | 738.25   | NA    | NA | NA       |
| 9    | 728.70 | 771.88 | \-355.35 | 710.70   | 27.54 | 1  | \<0.001  |

Likelihood ratio test for significance of random effect

**Test significance of (reaction norm) trial series intercept and time
slope correlation**

  - Run before testing random slope b/c result bkgrd slope intercept
    will determine which slope model to use

<!-- end list -->

``` r
# Fit model with random intercept-slope correlation
m_bkgrd_full <- refitModelRanEf(m_bkgrd_base,
                                 "(t_point_sc|ID) + (t_point_sc|trial_series)")  

# Fit model without random intercept-slope correlation
m_bkgrd_reduc <- refitModelRanEf(m_bkgrd_base,
                                 "(t_point_sc|ID) + (0 + t_point_sc|trial_series) +
                                  (1|trial_series)")

# LRT for significance bkgrd intercept-slope correlation
rpt_SigTest(m_bkgrd_full, m_bkgrd_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 10   | 718.32 | 766.29 | \-349.16 | 698.32   | NA    | NA | NA       |
| 11   | 704.22 | 756.98 | \-341.11 | 682.22   | 16.1  | 1  | \<0.001  |

Likelihood ratio test for significance of random effect

``` r
# print non-significant intercept-slope correlation
kable(as_tibble_col(attr(VarCorr(m_bkgrd_full)$trial_series,"correlation")[2,1],
              column_name = "Random intercept-slope variance estimate"),
      align = "l")
```

| Random intercept-slope variance estimate |
| :--------------------------------------- |
| 1                                        |

**Test significance of (reaction norm) time slope**

``` r
# Fit model with random slope
m_bkgrd_full <- refitModelRanEf(m_bkgrd_base,
                                 "(t_point_sc|ID) + (t_point_sc|trial_series)")  

# Fit model without random slope
m_bkgrd_reduc <- refitModelRanEf(m_bkgrd_base,
                                 "(t_point_sc|ID) + (1|trial_series)")

# LRT for significance bkgrd random slope
rpt_SigTest(m_bkgrd_full, m_bkgrd_reduc)
```

| npar | AIC    | BIC    | logLik   | deviance | Chisq | Df | p\_value |
| :--- | :----- | :----- | :------- | :------- | :---- | :- | :------- |
| 9    | 728.70 | 771.88 | \-355.35 | 710.70   | NA    | NA | NA       |
| 11   | 704.22 | 756.98 | \-341.11 | 682.22   | 28.49 | 2  | \<0.001  |

Likelihood ratio test for significance of random effect

**Final random effect model structure**

  - Removed trial series intercept and time slope correlation due to
    model convergence errors during bootstrapping

<!-- end list -->

``` r
# Refit without random effect covariances due to singular fits
m_bkgrd_rpt_final <- 
  glmer(bkgrd ~ t_point_sc * trial + sex +
          (0 + t_point_sc|ID) + (1|ID) + (0 + t_point_sc|trial_series) + (1|trial_series),  
        data = data_bkgrd_rpt,
        family = "binomial", 
        contrasts = list(trial = c(-1,1), sex = c(-1,1)),
        glmerControl(optimizer="bobyqa",
                    optCtrl = list(maxfun = 1e6)))

kable(as_tibble_col(as.character(c(formula(m_bkgrd_rpt_final))),
              column_name = "Final formula"))
```

| Final formula                                                                                                                        |
| :----------------------------------------------------------------------------------------------------------------------------------- |
| bkgrd \~ t\_point\_sc \* trial + sex + (0 + t\_point\_sc | ID) + (1 | ID) + (0 + t\_point\_sc | trial\_series) + (1 | trial\_series) |

### Variance and repeatability estimates

``` r
# Bootstrap confidence intervals
res_bootOut <- 
  bootMer(m_bkgrd_rpt_final,  # final rpt model 
          use.u = FALSE, 
          FUN = varEst_bkgrd,  # custom variance estimates function
          seed = 1,
          nsim = 10,  
          type = "parametric",
          parallel = "multicore",
          .progress = "txt")
```

    ## ====================================================================================================================

``` r
# Compile variance and repeatability estimates 
res_varEst <- 
  as_tibble(res_bootOut$t) %>%
  rename("Individual Intercept" = "V1",
         "Individual Slope" = "V2",
         "Series Intercept" = "V3",
         "Series Slope" = "V4",
         "Residual" = "V5",
         "Reaction Norm Intercept Repeatability" = "V6",
         "Reaction Norm Slope Repeatability" = "V7",
         "Long-term Repeatability" = "V8") %>%
  map_df(~ quantile(.x, c((1 - 0.95)/2, 1 - (1 - 0.95)/2),  # calculate 95CI from bootstrap
                  na.rm = TRUE), .id = "scale") %>%
  rename(Variable = scale, Lower95CI = `2.5%`, Upper95CI = `97.5%`) %>%
  mutate(Estimate = pluck(res_bootOut$t0)) %>%  # add R estimates
  relocate(Variable, Estimate) %>%  # reorder columns
  mutate_at(c(2:4), round, 2)  # # round estimates and errors to 2 digits

kable(res_varEst)
```

| Variable                              | Estimate | Lower95CI | Upper95CI |
| :------------------------------------ | -------: | --------: | --------: |
| Individual Intercept                  |     3.55 |      1.77 |      7.57 |
| Individual Slope                      |     0.80 |      0.29 |      2.82 |
| Series Intercept                      |     1.15 |      0.70 |      1.38 |
| Series Slope                          |     0.79 |      0.17 |      1.45 |
| Residual                              |     4.46 |      4.01 |      5.37 |
| Reaction Norm Intercept Repeatability |     0.75 |      0.62 |      0.86 |
| Reaction Norm Slope Repeatability     |     0.50 |      0.25 |      0.91 |
| Long-term Repeatability               |     0.33 |      0.18 |      0.53 |

## Activity

  - Number of lines crossed

**Create data set with NA values removed**

``` r
data_act_rpt <- 
  remove_dat_na(data_act, c("line_cross", "t_point_sc", "trial",
                            "sex", "ID", "trial_series", "day",
                            "obsID"))
```

### Random-effect selection

**Test significance of individual intercept (1|ID)**

``` r
# Fit base model
m_act_base <- 
  glmmTMB(line_cross ~ t_point_sc * trial + sex +
            (1|ID) + (t_point_sc|trial_series),
          data = data_act_rpt,
          contrasts = list(trial = c(-1,1), sex = c(-1,1)),
          ziformula= ~ 1,
          family="nbinom2")

# Fit model without individual ID
m_act_reduc <- refitModelRanEf(m_act_base,
                                 "(t_point_sc|trial_series)")

# LRT for significance act individual ID
rpt_SigTest(m_act_base, m_act_reduc)
```

| Df | AIC     | BIC     | logLik    | deviance | Chisq | Chi Df | p\_value |
| :- | :------ | :------ | :-------- | :------- | :---- | :----- | :------- |
| 10 | 6567.80 | 6625.75 | \-3273.90 | 6547.80  | NA    | NA     | NA       |
| 11 | 6547.79 | 6611.53 | \-3262.89 | 6525.79  | 22.01 | 1      | \<0.001  |

Likelihood ratio test for significance of random effect

**Test significance of individual intercept and time slope correlation
(t\_point\_sc|ID)**

  - Run before testing random slope because result act slope intercept
    will determine which slope model to use

<!-- end list -->

``` r
# Fit model with intercept-slope correlation
m_act_full <- refitModelRanEf(m_act_base,
                                 "(t_point_sc|ID) + (t_point_sc|trial_series)")  

# Fit model without intercept-slope m_act_base
m_act_reduc <- refitModelRanEf(m_act_base,
                                 "(0 + t_point_sc|ID) + (1|ID) +
                                  (t_point_sc|trial_series)")

# LRT for significance act intercept-slope correlation
rpt_SigTest(m_act_full, m_act_reduc)
```

| Df | AIC     | BIC     | logLik    | deviance | Chisq | Chi Df | p\_value |
| :- | :------ | :------ | :-------- | :------- | :---- | :----- | :------- |
| 12 | 6549.79 | 6619.32 | \-3262.89 | 6525.79  | NA    | NA     | NA       |
| 13 | 6551.08 | 6626.41 | \-3262.54 | 6525.08  | 0.71  | 1      | 0.2      |

Likelihood ratio test for significance of random effect

``` r
# print non-significant intercept-slope correlation
kable(as_tibble_col(attr(VarCorr(m_act_full)$cond$ID,"correlation")[1,2],
              column_name = "Random intercept-slope variance estimate"),
      align = "l")
```

| Random intercept-slope variance estimate |
| :--------------------------------------- |
| 0.9999215                                |

**Test significance of time slope (t\_point\_sc)**

``` r
# Fit model with random slope
m_act_full <- refitModelRanEf(m_act_base,
                                 "(0 + t_point_sc|ID) + (1|ID) +
                                  (t_point_sc|trial_series)") 

# Fit model without random slope
m_act_reduc <- refitModelRanEf(m_act_base,
                                 "(1|ID) + (t_point_sc|trial_series)")

# LRT for significance act random slope 
rpt_SigTest(m_act_full, m_act_reduc)
```

| Df | AIC     | BIC     | logLik    | deviance | Chisq | Chi Df | p\_value |
| :- | :------ | :------ | :-------- | :------- | :---- | :----- | :------- |
| 11 | 6547.79 | 6611.53 | \-3262.89 | 6525.79  | NA    | NA     | NA       |
| 12 | 6549.79 | 6619.32 | \-3262.89 | 6525.79  | 0     | 1      | 0.5      |

Likelihood ratio test for significance of random effect

``` r
# print non-significant random slope
kable(as_tibble_col(VarCorr(m_act_full)$cond$ID,
              column_name = "Random slope variance estimate"),
      align = "l")
```

| Random slope variance estimate |
| :----------------------------- |
| 3.10617e-07                    |

**Test significance of (reaction norm) trial series intercept**

``` r
# Fit model with random intercept
m_act_full <- refitModelRanEf(m_act_base,
                                 "(1|ID) + (1|trial_series)")

# Fit model without random intercept
m_act_reduc <- refitModelRanEf(m_act_base,
                                 "(1|ID)")

# LRT for significance act intercept
rpt_SigTest(m_act_full, m_act_reduc)
```

| Df | AIC     | BIC     | logLik    | deviance | Chisq | Chi Df | p\_value |
| :- | :------ | :------ | :-------- | :------- | :---- | :----- | :------- |
| 8  | 6659.68 | 6706.04 | \-3321.84 | 6643.68  | NA    | NA     | NA       |
| 9  | 6637.31 | 6689.46 | \-3309.66 | 6619.31  | 24.37 | 1      | \<0.001  |

Likelihood ratio test for significance of random effect

**Test significance of (reaction norm) trial series intercept and time
slope correlation**

  - Run before testing random slope b/c result act slope intercept will
    determine which slope model to use

<!-- end list -->

``` r
# Fit model with random intercept-slope correlation
m_act_full <- refitModelRanEf(m_act_base,
                                 "(1|ID) + (t_point_sc|trial_series)")  

# Fit model without random intercept-slope correlation
m_act_reduc <- refitModelRanEf(m_act_base,
                                 "(1|ID) + (0 + t_point_sc|trial_series) +
                                  (1|trial_series)")

# LRT for significance act intercept-slope correlation
rpt_SigTest(m_act_full, m_act_reduc)
```

| Df | AIC     | BIC     | logLik    | deviance | Chisq | Chi Df | p\_value |
| :- | :------ | :------ | :-------- | :------- | :---- | :----- | :------- |
| 10 | 6564.30 | 6622.25 | \-3272.15 | 6544.30  | NA    | NA     | NA       |
| 11 | 6547.79 | 6611.53 | \-3262.89 | 6525.79  | 18.52 | 1      | \<0.001  |

Likelihood ratio test for significance of random effect

**Test significance of (reaction norm) time slope**

``` r
# Fit model with random slope
m_act_full <- refitModelRanEf(m_act_base,
                                 "(1|ID) + (t_point_sc|trial_series)")  

# Fit model without random slope
m_act_reduc <- refitModelRanEf(m_act_base,
                                 "(1|ID) + (1|trial_series)")

# LRT for significance act random slope
rpt_SigTest(m_act_full, m_act_reduc)
```

| Df | AIC     | BIC     | logLik    | deviance | Chisq | Chi Df | p\_value |
| :- | :------ | :------ | :-------- | :------- | :---- | :----- | :------- |
| 9  | 6637.31 | 6689.46 | \-3309.66 | 6619.31  | NA    | NA     | NA       |
| 11 | 6547.79 | 6611.53 | \-3262.89 | 6525.79  | 93.53 | 2      | \<0.001  |

Likelihood ratio test for significance of random effect

**Final random effect model structure**

  - Removed trial series intercept and time slope correlation due to
    model convergence errors during bootstrapping

<!-- end list -->

``` r
m_act_rpt_final <- 
  glmmTMB(line_cross ~ t_point_sc * trial + sex +
            (1|ID) + (0 + t_point_sc|trial_series) + (1|trial_series),
          data = data_act_rpt,
          contrasts = list(trial = c(-1,1), sex = c(-1,1)),
          ziformula= ~ 1,
          family="nbinom2")

kable(as_tibble_col(as.character(c(formula(m_act_rpt_final))),
              column_name = "Final formula"))
```

| Final formula                                                                                                    |
| :--------------------------------------------------------------------------------------------------------------- |
| line\_cross \~ t\_point\_sc \* trial + sex + (1 | ID) + (0 + t\_point\_sc | trial\_series) + (1 | trial\_series) |

### Variance and repeatability estimates

``` r
# Bootstrap confidence intervals
res_bootOut <- 
  bootMer(m_act_rpt_final,  # final rpt model 
          use.u = FALSE, 
          FUN = varEst_act,  # custom variance estimates function
          seed = 1,
          nsim = 10,  
          type = "parametric",
          parallel = "multicore",
          .progress = "txt")
```

    ## ====================================================================================================================

``` r
# Compile variance and repeatability estimates 
res_varEst <- 
  as_tibble(res_bootOut$t) %>%
  rename("Individual Intercept" = "V1",
         "Series Intercept" = "V2",
         "Series Slope" = "V3",
         "Residual" = "V4",
         "Reaction Norm Intercept Repeatability" = "V5",
         "Long-term Repeatability" = "V6") %>%
  map_df(~ quantile(.x, c((1 - 0.95)/2, 1 - (1 - 0.95)/2),  # calculate 95CI from bootstrap
                  na.rm = TRUE), .id = "scale") %>%
  rename(Variable = scale, Lower95CI = `2.5%`, Upper95CI = `97.5%`) %>%
  mutate(Estimate = pluck(res_bootOut$t0)) %>%  # add R estimates
  relocate(Variable, Estimate) %>%  # reorder columns
  mutate_at(c(2:4), round, 2)  # # round estimates and errors to 2 digits

kable(res_varEst)
```

| Variable                              | Estimate | Lower95CI | Upper95CI |
| :------------------------------------ | -------: | --------: | --------: |
| Individual Intercept                  |     4.30 |      2.75 |      6.06 |
| Series Intercept                      |     0.80 |      0.23 |      0.63 |
| Series Slope                          |     1.38 |      0.71 |      1.74 |
| Residual                              |     4.36 |      3.03 |      6.85 |
| Reaction Norm Intercept Repeatability |     0.84 |      0.85 |      0.96 |
| Long-term Repeatability               |     0.40 |      0.34 |      0.51 |