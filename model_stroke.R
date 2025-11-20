model {
  # likelihood
  for (i in 1:N){
    Y[i] ~ dbern(p[i])
    logit(p[i]) <- beta0
                  + beta_site[Site[i]]
                  + beta_time[Time[i]]
                  + beta_age * Age[i]
                  + beta_gender[Gender[i]]
                  + beta_race[Race[i]]
                  + beta_transport * Transport[i]
                  + beta_notify * Notify[i]
                  + beta_TPA * hadTPA[i]
                  + beta_Thromb * hadThrombectomy[i]
                  + beta_tpaComp * tpaComplic[i]
                  + beta_thrComp * thrComplic[i]
  }
  
  # Notes for missing data models:
  # 1. All models should condition on site
  # 2. Notify must depend on Transport
  # 3. State conditional independence assumptions & justifications
  
  # Race | Site
  for (i in 1:N){
    Race[i] ~ dcat(pi_race[Site[i],1:3])
  }
  for (s in 1:9){
    pi_race[s,1:3] ~ ddirch(alpha[1:3])
  }

  alpha[1] <- 0.5  # Jeffreys prior
  alpha[2] <- 0.5
  alpha[3] <- 0.5
  
  # Transport | Race, Site
  for (i in 1:N){
    Transport[i] ~ dbern(pi_T[i])
    logit(pi_T[i]) <- gamma0
                      + gamma_site[Site[i]]
                      + gamma_race[Race[i]]
  }
  gamma0 ~ dnorm(0, 0.001)
  for (s in 1:9){ gamma_site[s] ~ dnorm(0, 0.001) }
  for (r in 1:3){ gamma_race[r] ~ dnorm(0, 0.001) }
  
  # Notify | Transport, Race, Site
  # Transport[i] = 1 if EMS, 0 if car
  for (i in 1:N){
    Notify[i] ~ dbern(pi_N[i])
    pi_N[i] <- pi_EMS[i] * Transport[i]
    logit(pi_EMS[i]) <- delta0
                        + delta_site[Site[i]]
                        + delta_race[Race[i]]
  }
  delta0 ~ dnorm(0, 0.001)
  for (s in 1:9){ delta_site[s] ~ dnorm(0, 0.001) }
  for (r in 1:3){ delta_race[r] ~ dnorm(0, 0.001) }
  
  # Age | Notify, Transport, Race, Site
  for (i in 1:N){
    Age[i] ~ dnorm(mu_A[i], tau_A)
    mu_A[i] <- eta0
              + eta_site[Site[i]]
              + eta_race[Race[i]]
              + eta_transport * Transport[i]
              + eta_notify * Notify[i]
  }
  eta0 ~ dnorm(0, 0.001)
  for (s in 1:9){ eta_site[s] ~ dnorm(0, 0.001) }
  for (r in 1:3){ eta_race[r] ~ dnorm(0, 0.001) }
  eta_transport ~ dnorm(0, 0.001)
  eta_notify ~ dnorm(0, 0.001)
  
  tau_A ~ dgamma(0.001, 0.001)
  sigma_A <- 1/sqrt(tau_A)
  
  
  # Priors for outcome model
  beta0 ~ dnorm(0, 0.001)
  for (s in 1:9){ beta_site[s] ~ dnorm(0, 0.001) }
  for (q in 1:8){ beta_time[q] ~ dnorm(0, 0.001) }
  beta_age ~ dnorm(0, 0.001)
  for (g in 1:2){ beta_gender[g] ~ dnorm(0, 0.001) }
  for (r in 1:3){ beta_race[r] ~ dnorm(0, 0.001) }
  
  # Binary covariates / treatments
  beta_transport ~ dnorm(0, 0.001)
  beta_notify    ~ dnorm(0, 0.001)
  beta_TPA       ~ dnorm(0, 0.001)
  beta_Thromb    ~ dnorm(0, 0.001)
  beta_tpaComp   ~ dnorm(0, 0.001)
  beta_thrComp   ~ dnorm(0, 0.001)
}