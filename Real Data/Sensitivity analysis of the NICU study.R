expit <- function(x) return(exp(x)/(1+exp(x)))
library(boot)
#dt = read.csv('matches_reduced_strongerIV2005.csv')
dt = read.csv('allyears_matched_originalIV.csv')
dt_strengthen = read.csv('allyears_matched_strongerIV.csv')

# Impute
dt_die.e = dt$death.e + dt$death_fatel.e
dt_die.u = dt$death.u + dt$death_fatel.u
survisor_losI.e = dt$losI.e[which(dt_die.e == 0)]
survisor_losI.u = dt$losI.u[which(dt_die.u == 0)]
imputed_value = quantile(c(survisor_losI.e, survisor_losI.u), 0.99)

dt$losI.e[which(dt_die.e == 1)] = imputed_value
dt$losI.u[which(dt_die.u == 1)] = imputed_value

dt = dt[complete.cases(dt),]

dt_strengthen_die.e = dt_strengthen$death.e + dt_strengthen$death_fatel.e
dt_strengthen_die.u = dt_strengthen$death.u + dt_strengthen$death_fatel.u


dt_strengthen$losI.e[which(dt_strengthen_die.e == 1)] = imputed_value
dt_strengthen$losI.u[which(dt_strengthen_die.u == 1)] = imputed_value

dt_strengthen = dt_strengthen[complete.cases(dt_strengthen),]


# Fit an overall population model using all data
library(dplyr)
dt_1 = dt %>%
  dplyr::select(dif.travel.time.e, bthwght.e, gestage.e, gestdiabetes.e, singlebirth.e, parity.e, ageM.e,
                medu.e, white.e, race.mis.e, below_poverty.e, high_level_NICU.e, losI.e)
colnames(dt_1) <- c('dif.travel.time', 'bthwght', 'gestage', 'gestdiabetes',
                    'singlebirth', 'parity','ageM', 'medu', 'white', 'race.mis', 'below_poverty', 'high_level_NICU', 'LOSI')
dt_2 = dt %>%
  dplyr::select(dif.travel.time.u, bthwght.u, gestage.u, gestdiabetes.u, singlebirth.u, parity.u, ageM.u,
                medu.u, white.u, race.mis.u, below_poverty.u, high_level_NICU.u, losI.u)

colnames(dt_2) <- c('dif.travel.time', 'bthwght', 'gestage', 'gestdiabetes',
                    'singlebirth', 'parity','ageM', 'medu', 'white', 'race.mis', 'below_poverty', 'high_level_NICU', 'LOSI')

dt_all = rbind(dt_1, dt_2)
lm_model = lm(dif.travel.time ~ bthwght + gestage + gestdiabetes + singlebirth + parity + ageM +
                medu + white + race.mis + below_poverty, data = dt_all)

p = predict(lm_model, newdata = dt_all)
z_tilde = scale(dt_all$dif.travel.time - p)
dt_all$z_tilde = z_tilde



p_1 = predict(lm_model, newdata = dt_1)
z_tilde_e = scale(dt_1$dif.travel.time - p_1)

p_2 = predict(lm_model, newdata = dt_2)
z_tilde_u = scale(dt_2$dif.travel.time - p_2)

dt$z_tt.e = z_tilde_e
dt$z_tt.u = z_tilde_u


dt_1 = dt_strengthen %>%
  dplyr::select(dif.travel.time.e, bthwght.e, gestage.e, gestdiabetes.e, singlebirth.e, parity.e, ageM.e,
                medu.e, white.e, race.mis.e, below_poverty.e, high_level_NICU.e, losI.e)
colnames(dt_1) <- c('dif.travel.time', 'bthwght', 'gestage', 'gestdiabetes',
                    'singlebirth', 'parity','ageM', 'medu', 'white', 'race.mis', 'below_poverty', 'high_level_NICU', 'LOSI')
dt_2 = dt_strengthen %>%
  dplyr::select(dif.travel.time.u, bthwght.u, gestage.u, gestdiabetes.u, singlebirth.u, parity.u, ageM.u,
                medu.u, white.u, race.mis.u, below_poverty.u, high_level_NICU.u, losI.u)

colnames(dt_2) <- c('dif.travel.time', 'bthwght', 'gestage', 'gestdiabetes',
                    'singlebirth', 'parity','ageM', 'medu', 'white', 'race.mis', 'below_poverty', 'high_level_NICU', 'LOSI')

dt_all = rbind(dt_1, dt_2)
lm_model = lm(dif.travel.time ~ bthwght + gestage + gestdiabetes + singlebirth + parity + ageM +
                medu + white + race.mis + below_poverty, data = dt_all)

p = predict(lm_model, newdata = dt_all)
z_tilde = scale(dt_all$dif.travel.time - p)
dt_all$z_tilde = z_tilde


p_1 = predict(lm_model, newdata = dt_1)
z_tilde_e = scale(dt_1$dif.travel.time - p_1)

p_2 = predict(lm_model, newdata = dt_2)
z_tilde_u = scale(dt_2$dif.travel.time - p_2)

dt_strengthen$z_tt.e = z_tilde_e
dt_strengthen$z_tt.u = z_tilde_u



# Estimate sigma
estimate_sigma <- function(lambda_0, lambda_1, delta, dt_all){
  prob = expit(lambda_0 + lambda_1 * dt_all$z_tilde)

  sigma_est_vec = numeric(10)
  for (i in 1:10){
    dt_all$u = rbinom(length(prob), 1, prob)
    transformed_lost = dt_all$LOSI - delta * dt_all$u
    model_2 = lm(transformed_lost ~ bthwght + gestage + gestdiabetes + singlebirth + parity + ageM +
                   medu + white + race.mis + below_poverty + high_level_NICU, data = dt_all)
    sigma_est_vec[i] = sigma(model_2)
  }

  return(mean(sigma_est_vec))
}


sensitivity_interval <- function(dt, delta, lambda_0, lambda_1, sigma_est){
  I = dim(dt)[1]
  E_e = rbinom(I, 1, prob = expit(lambda_0 + lambda_1*dt$z_tt.e))
  E_u = rbinom(I, 1, prob = expit(lambda_0 + lambda_1*dt$z_tt.u))


  compliance = abs(mean(dt$high_level_NICU.e, na.rm = TRUE) - mean(dt$high_level_NICU.u, na.rm = TRUE))
 
  estimate_bias = delta*abs(mean(E_e) - mean(E_u))/compliance
  
  wald = (mean(dt$losI.e) - mean(dt$losI.u))/compliance
  
  point_est = wald - estimate_bias
  sd_point_est = sqrt(2)*sigma_est/(sqrt(I)*compliance)
  
  return(c(point_est, sd_point_est^2))

}

# Solve for lambda_0 for a given tau and lambda_1
solve_lambda_0 <- function(tau, lambda_1){
  func <- function(x, lambda_0, lambda_1) expit(lambda_0 + lambda_1*x)*(1/sqrt(2*pi))*exp(-x^2)
  func_2 <- function(lambda_0, lambda_1, tau) integrate(func, lower = 0, upper = 10, lambda_0 = lambda_0, lambda_1 = lambda_1)$value -
    integrate(func, lower = -10, upper = 0, lambda_0 = lambda_0, lambda_1 = lambda_1)$value - 2*tau

  tryCatch({uniroot(func_2, c(0, 10), extendInt = c('yes'), lambda_1 = lambda_1, tau = tau)$root},
           error=function(error_message) {
             message(" Error in integrate(func, lower = 0, upper = 10, lambda_0 = lambda_0, lambda_1 = lambda_1) :
                     non-finite function value ")
             return(NA)
           })
}


sensitivity_interval_mi <- function(dt, dt_all, tau, lambda_1, delta){
  point_est_vec = numeric(5)
  var_vec = numeric(5)
  for (i in 1:5){
    lambda_0 = solve_lambda_0(tau, lambda_1)
    sigma_est = estimate_sigma(lambda_0, lambda_1, delta, dt_all)
    res_imputation = sensitivity_interval(dt, delta, lambda_0, lambda_1, sigma_est)
    point_est_vec[i] = res_imputation[1]
    var_vec[i] = res_imputation[2]
  }
  point_est = mean(point_est_vec)
  total_var = mean(var_vec) + (1 + 1/5)*var(point_est_vec)
  cat(mean(var_vec), (1 + 1/5)*var(point_est_vec), '\n')
  return(c(point_est, sqrt(total_var)))
}



find_delta <- function(dt, dt_all, tau, lambda_1, delta_0, tol = 1){
  delta = delta_0
  ASI_low = 100
  while (delta < 1000){
    ASI_low = sensitivity_interval_mi(dt, dt_all, tau, lambda_1, delta)[1]
  if (ASI_low < 0) return(delta)
    delta = delta + tol
  }
  return(delta)
}

###############################################################
# Not strengthen

delta_0_0_1 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_1 = c(delta_0_0_1, find_delta(dt, dt_all, 0.01, lambda_1, 20))
}

delta_0_0_1_5 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_1_5 = c(delta_0_0_1_5, find_delta(dt, dt_all, 0.015, lambda_1, 20))
}

delta_0_0_2 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_2 = c(delta_0_0_2, find_delta(dt, dt_all, 0.02, lambda_1, 30))
}

delta_0_0_2_5 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_2_5 = c(delta_0_0_2_5, find_delta(dt, dt_all, 0.025, lambda_1, 40))
}

delta_0_0_3 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_3 = c(delta_0_0_3, find_delta(dt, dt_all, 0.03, lambda_1, 40))
}

delta_0_0_3_5 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_3_5 = c(delta_0_0_3_5, find_delta(dt, dt_all, 0.035, lambda_1, 30))
}

delta_0_0_4 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_4 = c(delta_0_0_4, find_delta(dt, dt_all, 0.04, lambda_1, 20))
}


tau = rep(c(0.01, 0.015, 0.02,0.025,0.03,0.035,0.04), each = 11)
lambda_1 = rep(seq(1,3,0.2), 7)
delta_max = c(delta_0_0_1, delta_0_0_1_5, delta_0_0_2, delta_0_0_2_5, delta_0_0_3, delta_0_0_3_5, delta_0_0_4)
delta_max = delta_max/15.8
mat_plot = data.frame(tau, lambda_1, delta_max)


v <- ggplot(mat_plot, aes(lambda_1, tau, z = delta_max)) + geom_raster(aes(fill = delta_max)) +
  geom_contour(colour = "white") + scale_fill_gradient(name = expression(delta[sup]), low = "grey", high = "black") +
  xlab(expression(lambda[1])) + ylab(expression(tau)) +
  theme(legend.title = element_text(size=14),
        axis.title = element_text(size=16))




###############################################################
# Strengthen

# Not strengthen

dt = dt_strengthen
delta_0_0_1 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_1 = c(delta_0_0_1, find_delta(dt, dt_all, 0.01, lambda_1, 20))
}

delta_0_0_1_5 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_1_5 = c(delta_0_0_1_5, find_delta(dt, dt_all, 0.015, lambda_1, 20))
}

delta_0_0_2 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_2 = c(delta_0_0_2, find_delta(dt, dt_all, 0.02, lambda_1, 20))
}

delta_0_0_2_5 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_2_5 = c(delta_0_0_2_5, find_delta(dt, dt_all, 0.025, lambda_1, 20))
}

delta_0_0_3 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_3 = c(delta_0_0_3, find_delta(dt, dt_all, 0.03, lambda_1, 20))
}

delta_0_0_3_5 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_3_5 = c(delta_0_0_3_5, find_delta(dt, dt_all, 0.035, lambda_1, 20))
}

delta_0_0_4 = NULL
for (lambda_1 in seq(1, 3, 0.2)){
  delta_0_0_4 = c(delta_0_0_4, find_delta(dt, dt_all, 0.04, lambda_1, 20))
}


tau = rep(c(0.01, 0.015, 0.02,0.025,0.03,0.035,0.04), each = 11)
lambda_1 = rep(seq(1,3,0.2), 7)
delta_max = c(delta_0_0_1, delta_0_0_1_5, delta_0_0_2, delta_0_0_2_5, delta_0_0_3, delta_0_0_3_5, delta_0_0_4)
delta_max = delta_max/15.8
mat_plot = data.frame(tau, lambda_1, delta_max)


v <- ggplot(mat_plot, aes(lambda_1, tau, z = delta_max)) + geom_raster(aes(fill = delta_max)) +
  geom_contour(colour = "white") + scale_fill_gradient(name = expression(delta[sup]), low = "grey", high = "black") +
  xlab(expression(lambda[1])) + ylab(expression(tau)) +
  theme(legend.title = element_text(size=14),
        axis.title = element_text(size=16))











