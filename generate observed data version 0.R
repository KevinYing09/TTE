#in this version, the function generate.obs.data() does not determine eligibility. We leave this to be determined after the data is generated.
# Load required libraries
library(dplyr)

# Set seed
set.seed(123)

generate.obs.data <- function(n, J, beta0.trt0,beta.adherence, beta0.outcome, beta.trt.effect, beta0.censoring) {
  
  # Initialize baseline data
  mean.age=35; sd.age=12
  baseline <- data.frame(
    id = 1:n,
    X1_0 = rbinom(n, 1, 0.5),
    X2_0 = rnorm(n, 0, 1),
    age_0 = rnorm(n, mean.age, sd.age),
    X3 = rbinom(n, 1, 0.5),
    X4 = rnorm(n, 0, 1)
  )%>%rowwise()%>%mutate( logit_pi0 =beta0.trt0,
                          A0 = rbinom(1, 1, plogis(logit_pi0)))%>%select(-logit_pi0)%>%ungroup()
  
  # Expand to long format
  long_data <- baseline %>%
    slice(rep(1:n(), each = J)) %>%
    mutate(time = rep(0:(J - 1), times = n)) %>%
    group_by(id) %>%
    mutate(
      age = age_0 + time,
      ages = (age - mean.age) / sd.age,
      A_prev = NA,
      Y_prev = NA,
      C_prev = NA,
      A=if_else(time==0, A0, NA),
      Y = if_else(time == 0, 0, NA),
      C = if_else(time == 0, 0, NA),
      #eligible = 1,
      X1 = if_else(time == 0, X1_0, NA),
      X2 = if_else(time == 0, X2_0, NA)
    ) %>%
    ungroup()
  
  # Simulate time-varying covariates, treatment, outcome, censoring, and eligibility
  long_data <- long_data %>%
    group_by(id) %>%
    group_modify(~ {
      df <- .x
      for (j in 2:nrow(df)) {
        # if (j > 1) {
        df$X1[j] <- rbinom(1, 1, plogis(-df$A[j - 1]))
        df$X2[j] <- rnorm(1, mean = -0.3 * df$A[j - 1], sd = 1)
        df$A_prev[j] <- df$A[j - 1]
        df$Y_prev[j] <- df$Y[j - 1]
        df$C_prev[j] <- df$C[j - 1]
        # Treatment
        logit_pi <-beta.adherence* (df$A_prev[j]-0.5 )+ 0.5 * df$X1[j] + 0.5 * df$X2[j] - 0.2 * df$X3[j] + df$X4[j] - 0.3 * df$ages[j]
        df$A[j] <- rbinom(1, 1, plogis(logit_pi))
        #}
        #df$eligible[j] <- as.integer(df$age[j] >= 18 & all(df$A[1:j - 1] == 0) & all(df$Y[1:j - 1] == 0))
        # Outcome
        if (df$Y_prev[j] == 0) {
          logit_lambda <- beta0.outcome +beta.trt.effect * df$A[j] + 0.5 * df$X1[j] + 0.5 * df$X2[j] + df$X3[j] + df$X4[j] + 0.5 * df$ages[j]
          df$Y[j] <- rbinom(1, 1, plogis(logit_lambda))
        } else {
          df$Y[j] <- 1  # once outcome is positive, remain positive
        }
        
        # Censoring
        if (df$C_prev[j] == 0){
          logit_q <- beta0.censoring - df$A_prev[j] - 0.5 * df$X1[j] + 0.5 * df$X2[j] - 0.2 * df$X3[j] + 0.2 * df$X4[j] - df$ages[j]
          df$C[j] <- rbinom(1, 1, plogis(logit_q))
        } else {
          df$C[j] <- 1  # once censored, remain censored
        }
        
      }
      df
    }) %>%
    ungroup()
  #drop data after censoring or event
  t0=long_data%>%filter(C==1 | Y==1)%>% group_by(id)%>%summarise(t0=min(time))
  long_data=long_data%>%left_join(t0)%>%filter(is.na(t0)| time<=t0)%>%select(-t0)
  
  # #drop data before first eligible
  # long_data=long_data %>%
  #   group_by(id) %>%
  #   mutate(first_eligible = min(time[eligible == 1], na.rm = TRUE)) %>%
  #   filter(time >= first_eligible) %>%
  #   ungroup() %>%
  #   select(-first_eligible)
  # Ensure all variables are numeric
  long_data <- long_data %>%
    mutate(across(c(id, time, A, Y, C,  X1, X2, X3, X4, age), as.numeric))%>%
    select(id,time,X1,X2,X3,X4, age,A, Y, C) 
  long_data
  
}

#checking
#simdata_censored=generate.obs.data(n=10000,J=10,beta0.trt0=0,beta.adherence=10, beta0.outcome=-5, beta.trt.effect=-1.2, beta0.censoring=-1)
