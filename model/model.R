
# this script contains the model function

model_function = function(t, pop, param, commu_FOI, NCD_rate) {
  
  with(as.list(c(param, pop)), {
    
    N=S+E+Ia+P+Is+R+S_c+E_c+Ia_c+P_c+Is_c+R_c
    
    # Teleworking activation
    # Teleworking is only activated if we are past the activation time
    alpha = ifelse(t >= t_alpha, alpha, 0)
    
    # Community force of infection
    # Extrapolate at time t in the solver from the provided interpolating function commu_FOI()
    lambda_v = commu_FOI(t)
    
    # Total force of infection
    # workplace infections + weekend infections + telework infections
    # For a frequency-dependent model, R0 = beta*alpha/gamma_a
    # For a density-dependent model, R0 = beta*alpha*N/gamma_a
    
    if (alpha == 1) {
      beta = 0
    } else {
      beta = R0 * rho * gamma_a / ((1-prop_a) * gamma_a + rho * nu * prop_a) 
    }
    
    lambda = 5/7 * (1-alpha) * beta * (nu*(Ia+Ia_c) + (P+P_c))/(N-Is)  + 
      5/7 * alpha * epsilon * lambda_v + 
      2/7 * lambda_v
    
    # Rate of work-related chronic disease
    # Extrapolate for the proportion of telework alpha from the provided function chronic_rate()
    # Note this approach allows us to flexibly change alpha during the simulation, to
    # reflect scenarios where telework might be introduced and then lifted
    omega = baseline_NCD*NCD_rate(alpha)
    
    #force omega = 0 if alpha = 0
    # if(alpha == 0) omega = 0
    
    # Individuals without work-related chronic disease ####
    # Susceptible
    # - workplace infections - weekend infections - telework infections - chronic disease incidence
    dS = -lambda*S - omega*S
    
    # Exposed (non-infectious)
    # + infections - progressions to infectious - chronic disease incidence
    dE = lambda*S - sigma*E - omega*E
    
    # Infected (asymptomatic)
    # + progressions from exposed - recoveries - chronic disease incidence
    dIa = sigma*E*prop_a - gamma_a*Ia - omega*Ia
    
    # Pre-symptomatic (infectious)
    # + progressions from exposed - progressions to symptomatic - chronic disease incidence
    dP = sigma*E*(1-prop_a) - rho*P - omega*P
    
    # Infected (symptomatic)
    # Note: no chronic disease incidence, as we assume these individuals are not working
    # + progressions from pre-symptomatic - recoveries
    dIs = rho*P - gamma_s*Is
    
    # Recovered
    # + recoveries - chronic disease incidence
    dR = Ia*gamma_a + Is*gamma_s - omega*R
    
    
    
    # Individuals without work-related chronic disease ####
    # Susceptible
    # - workplace infections - weekend infections - telework infections + chronic disease incidence
    dS_c = -lambda*S_c + omega*S
    
    # Exposed (non-infectious)
    # + infections - progressions to infectious + chronic disease incidence
    dE_c = lambda*S_c - sigma*E_c + omega*E
    
    # Infected (asymptomatic)
    # + progressions from exposed - recoveries + chronic disease incidence
    dIa_c = sigma*E_c*prop_a - gamma_a*Ia_c + omega*Ia
    
    # Pre-symptomatic (infectious)
    # + progressions from exposed - progressions to symptomatic + chronic disease incidence
    dP_c = sigma*E_c*(1-prop_a) - rho*P_c + omega*P
    
    # Infected (symptomatic)
    # + progressions from pre-symptomatic - recoveries
    dIs_c = rho*P_c - gamma_s*Is_c
    
    # Recovered
    # + recoveries + chronic disease incidence
    dR_c = Ia_c*gamma_a + Is_c*gamma_s + omega*R
    
    list(c(dS, dE, dIa, dP, dIs, dR, dS_c, dE_c, dIa_c, dP_c, dIs_c, dR_c))
    
  })
  
}
