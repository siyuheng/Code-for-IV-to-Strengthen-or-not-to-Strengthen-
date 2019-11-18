library(MASS)

bisection <- function(f, a, b, n = 20, tol = 1e-3) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  if (f(a)>0 || f(b)<0) {
    return(99)
  } 
  
  for (j in 1:n) {
    print(j)
    c <- as.integer((a+b)/2) # Calculate midpoint
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
    # function and return the root.
    if ((f(c) == 0) || ((b - a) / 2) < tol) {
      return(c)
    }
    
    # If another iteration is required, 
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(f(c)) == sign(f(a)), 
           a <- c,
           b <- c)
  }
  # If the max number of iterations is reached and no root has been found, 
  # return message and end function.
  print('Too many iterations')
}

power_wilcoxon_compare<-function(I, alpha, beta, beta_0, iA, iC, iN, S){
  
  gamma=1
  
  EW_Null=(gamma/(gamma+1))*I*(I+1)/2
  VW_Null=(gamma/((1+gamma)*(1+gamma)))*I*(I+1)*(2*I+1)/6
  CV_Null=EW_Null+qnorm(1-alpha)*sqrt(VW_Null)
  
  power_W=0
  counting_W=0
  T_W=rep(0,S)
  
  prop_1=(iA+iC)*iA+iN*(iN+iC)
  prop_2=iN*iA
  prop_3=iN*iA+iC
  for(i in 1:S){
    randomsample_W_false=rep(0,I)
    indicator_W=rep(0,I)
    # Sample from normal mixture
    components <- sample(1:3,prob=c(prop_1, prop_2, prop_3), size = I, replace=TRUE)
    mus <- c(0, beta_0 - beta, beta - beta_0)
    sds <- sqrt(c(1, 1, 1))
    randomsample_W_false <- rnorm(n = I, mean = mus[components], sd=sds[components])
    rank_W=rank(abs(randomsample_W_false))
    for(j in 1:I){
      if(randomsample_W_false[j]>0){
        indicator_W[j]=1
      }
    }
    T_W[i]=sum(rank_W*indicator_W)
    if(T_W[i]>CV_Null){
      counting_W=counting_W+1
    }
  }
  power_W=counting_W/S
  return(power_W)
}


power_wilcoxon_compare_laplace<-function(I, alpha, beta, beta_0, iA, iC, iN, S){
  
  gamma=1
  
  EW_Null=(gamma/(gamma+1))*I*(I+1)/2
  VW_Null=(gamma/((1+gamma)*(1+gamma)))*I*(I+1)*(2*I+1)/6
  CV_Null=EW_Null+qnorm(1-alpha)*sqrt(VW_Null)
  
  power_W=0
  counting_W=0
  T_W=rep(0,S)
  
  prop_1=(iA+iC)*iA+iN*(iN+iC)
  prop_2=iN*iA
  prop_3=iN*iA+iC
  for(i in 1:S){
    randomsample_W_false=rep(0,I)
    indicator_W=rep(0,I)
    
    # Sample from Laplace mixture
    components <- sample(1:3,prob=c(prop_1, prop_2, prop_3), size = I, replace=TRUE)
    mus <- c(0, beta_0 - beta, beta - beta_0)
    sds <- sqrt(c(1, 1, 1))
    randomsample_W_false <- rlaplace(I, mus[components], sqrt(2)/2) 
    rank_W=rank(abs(randomsample_W_false))
    
    rank_W=rank(abs(randomsample_W_false))
    for(j in 1:I){
      if(randomsample_W_false[j]>0){
        indicator_W[j]=1
      }
    }
    T_W[i]=sum(rank_W*indicator_W)
    if(T_W[i]>CV_Null){
      counting_W=counting_W+1
    }
  }
  power_W=counting_W/S
  return(power_W)
}


power_sign_compare<-function(I, alpha, beta, beta_0, iA, iC, iN, S){
  power_W=0
  counting_W=0
  T_W=rep(0,S)
  
  prop_1=(iA+iC)*iA+iN*(iN+iC)
  prop_2=iN*iA
  prop_3=iN*iA+iC
  for(i in 1:S){
    randomsample_W_false=rep(0,I)
    indicator_W=rep(0,I)
    
    # Sample from normal mixture
    components <- sample(1:3,prob=c(prop_1, prop_2, prop_3), size = I, replace=TRUE)
    mus <- c(0, beta_0 - beta, beta - beta_0)
    sds <- sqrt(c(1, 1, 1))
    randomsample_W_false <- rnorm(n = I, mean = mus[components], sd=sds[components])
    
    T_W[i]=(sum(randomsample_W_false>0)-I/2)/sqrt(I/4)
    if(T_W[i]>qnorm(1-alpha)){
      counting_W=counting_W+1
    }
  }
  power_W=counting_W/S
  return(power_W)
}


power_sign_compare_laplace<-function(I, alpha, beta, beta_0, iA, iC, iN, S){
  power_W=0
  counting_W=0
  T_W=rep(0,S)
  
  prop_1=(iA+iC)*iA+iN*(iN+iC)
  prop_2=iN*iA
  prop_3=iN*iA+iC
  for(i in 1:S){
    randomsample_W_false=rep(0,I)
    indicator_W=rep(0,I)
    for(j in 1:I){
      randomnumber=runif(1)
      prop_1=(iA+iC)*iA+iN*(iN+iC)
      prop_2=iN*iA
      prop_3=iN*iA+iC
      indic_1=as.numeric(randomnumber<=prop_1)
      indic_2=as.numeric(randomnumber>prop_1)*as.numeric(randomnumber<(1-prop_3))
      indic_3=as.numeric(randomnumber>=(1-prop_3))
      compo_1=rlaplace(1, location = 0, scale = sqrt(2)/2)
      compo_2=rlaplace(1, location = beta_0 - beta, scale = sqrt(2)/2)
      compo_3=rlaplace(1, location = beta - beta_0, scale = sqrt(2)/2)
      randomsample_W_false[j]=indic_1*compo_1 + indic_2*compo_2 + indic_3*compo_3
    }
    T_W[i]=(sum(randomsample_W_false>0)-I/2)/sqrt(I/4)
    if(T_W[i]>qnorm(1-alpha)){
      counting_W=counting_W+1
    }
  }
  power_W=counting_W/S
  return(power_W)
}

samplesize_wilc<-function(alpha, pow, beta, beta_0, iA, iC, iN, S, sample_size_lower, sample_size_upper){
  equation<-function(I){
    y=power_wilcoxon_compare(I, alpha, beta, beta_0, iA, iC, iN, S)-pow
    return(y)
  }
  ds<-bisection(equation, a=sample_size_lower, b=sample_size_upper, n=15, tol = 5)
  return(ds)
}


samplesize_laplace_wilc<-function(alpha, pow, beta, beta_0, iA, iC, iN, S, sample_size_lower, sample_size_upper){
  equation<-function(I){
    y=power_wilcoxon_compare_laplace(I, alpha, beta, beta_0, iA, iC, iN, S)-pow
    return(y)
  }
  ds<-bisection(equation, a=sample_size_lower, b=sample_size_upper, n=15, tol=1)
  return(ds)
}


samplesize_sign<-function(alpha, pow, beta, beta_0, iA, iC, iN, S, sample_size_lower, sample_size_upper){
  equation<-function(I){
    y=power_sign_compare(I, alpha, beta, beta_0, iA, iC, iN, S)-pow
    return(y)
  }
  ds<-bisection(equation, a=sample_size_lower, b=sample_size_upper, n=15, tol=1)
  return(ds)
}


samplesize_laplace_sign<-function(alpha, pow, beta, beta_0, iA, iC, iN, S, sample_size_lower, sample_size_upper){
  equation<-function(I){
    y=power_sign_compare_laplace(I, alpha, beta, beta_0, iA, iC, iN, S)-pow
    return(y)
  }
  ds<-bisection(equation, a=sample_size_lower, b=sample_size_upper, n=15, tol=1)
  return(ds)
}

# Plot sample size against power
get_power_curve <- function(iota_c, lower, upper, S){
  sample_vec = seq(lower, upper, 10)
  power_vec = NULL
  for (I in sample_vec){
    power_temp = power_wilcoxon_compare(I, 0.05, 0.1, 0, 0, iota_c, 1 - iota_c, S)
    power_vec = c(power_vec, power_temp)
  }
  return(power_vec)
}

###########################################################################
# Normal error

# pair I
res_0_5 = get_power_curve(0.5, 100, 4000, 10000)
res_0_6 = get_power_curve(0.6, 100, 4000, 10000)

# pair II
res_0_4 = get_power_curve(0.4, 100, 6000, 10000)
res_0_7 = get_power_curve(0.7, 100, 6000, 10000)

# pair III
res_0_3 = get_power_curve(0.3, 100, 9000, 10000)
res_0_8 = get_power_curve(0.8, 100, 9000, 10000)


# Make plots
dt_pair_1 = data.frame(sample_size = rep(seq(100, 4000, 10), 2), power = c(res_0_5, res_0_6), 
                       IV = c(rep('Compliance = 0.5', length(seq(100, 4000,10))),
                              rep('Compliance = 0.6', length(seq(100, 4000,10)))))

p <- ggplot(data = dt_pair_1, aes(x = sample_size, y = power, colour = IV, fill = IV)) + 
  geom_smooth(method = "loess", se = FALSE) +
  xlim(c(0, 4000)) + ylim(c(0, 1)) +
  labs(x = 'sample size', y = 'power') + 
  theme(legend.position = c(0.85, 0.2), 
        legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        axis.title = element_text(size=16))
ggsave('plot_power_pair_1.png', device = 'png')



dt_pair_2 = data.frame(sample_size = rep(seq(100, 6000, 10), 2), power = c(res_0_4, res_0_7), 
                       IV = c(rep('Compliance = 0.4', length(seq(100, 6000,10))),
                              rep('Compliance = 0.7', length(seq(100, 6000,10)))))

p <- ggplot(data = dt_pair_2, aes(x = sample_size, y = power, colour = IV, fill = IV)) + 
  geom_smooth(method = "loess", se = FALSE) +
  xlim(c(0, 6000)) + ylim(c(0, 1)) +
  labs(x = 'sample size', y = 'power') + 
  theme(legend.position = c(0.85, 0.2), 
        legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        axis.title = element_text(size=16))
ggsave('plot_power_pair_2.png', device = 'png')


dt_pair_3 = data.frame(sample_size = rep(seq(100, 9000, 10), 2), power = c(res_0_3, res_0_8), 
                       IV = c(rep('Compliance = 0.3', length(seq(100, 9000,10))),
                              rep('Compliance = 0.8', length(seq(100, 9000,10)))))

p <- ggplot(data = dt_pair_3, aes(x = sample_size, y = power, colour = IV, fill = IV)) + 
  geom_smooth(method = "loess", se = FALSE) +
  xlim(c(0, 9000)) + ylim(c(0, 1)) +
  labs(x = 'sample size', y = 'power') + 
  theme(legend.position = c(0.85, 0.2), 
        legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        axis.title = element_text(size=16))
ggsave('plot_power_pair_3.png', device = 'png')

############################################################################
# Lapalce error
get_power_curve_laplace <- function(iota_c, lower, upper, S){
  sample_vec = seq(lower, upper, 10)
  power_vec = NULL
  for (I in sample_vec){
    power_temp = power_wilcoxon_compare_laplace(I, 0.05, 0.1, 0, 0, iota_c, 1 - iota_c, S)
    power_vec = c(power_vec, power_temp)
  }
  return(power_vec)
}

# pair I
res_0_5 = get_power_curve_laplace(0.5, 100, 4000, 10000)
res_0_6 = get_power_curve_laplace(0.6, 100, 4000, 10000)

# pair II
res_0_4 = get_power_curve_laplace(0.4, 100, 6000, 10000)
res_0_7 = get_power_curve_laplace(0.7, 100, 6000, 10000)

# pair III
res_0_3 = get_power_curve_laplace(0.3, 100, 9000, 10000)
res_0_8 = get_power_curve_laplace(0.8, 100, 9000, 10000)


# Match sample size
num_ind = 35

sp_0_3 = NULL
sp_0_8 = NULL

for (power in seq(min(res_0_3) + 0.15, max(res_0_3), 0.02)){
  sp_1 = which.min(abs(res_0_3 - power))
  sp_2 = which.min(abs(res_0_8 - power))
  sp_0_3 = c(sp_0_3, 100 + (sp_1 - 1)*10)
  sp_0_8 = c(sp_0_8, 100 + (sp_2 - 1)*10)
}
ind_3_8 = sample(seq(1, length(sp_0_3)), num_ind, replace = FALSE)
sp_0_3 = sp_0_3[ind_3_8]
sp_0_8 = sp_0_8[ind_3_8]

sp_0_4 = NULL
sp_0_7 = NULL

for (power in seq(min(res_0_4) + 0.1, max(res_0_4), 0.02)){
  sp_1 = which.min(abs(res_0_4 - power))
  sp_2 = which.min(abs(res_0_7 - power))
  sp_0_4 = c(sp_0_4, 100 + (sp_1 - 1)*10)
  sp_0_7 = c(sp_0_7, 100 + (sp_2 - 1)*10)
}
ind_4_7 = sample(seq(1, length(sp_0_4)), num_ind, replace = FALSE)
sp_0_4 = sp_0_4[ind_4_7]
sp_0_7 = sp_0_7[ind_4_7]

sp_0_5 = NULL
sp_0_6 = NULL

for (power in seq(min(res_0_5) + 0.1, max(res_0_5), 0.02)){
  sp_1 = which.min(abs(res_0_5 - power))
  sp_2 = which.min(abs(res_0_6 - power))
  sp_0_5 = c(sp_0_5, 100 + (sp_1 - 1)*10)
  sp_0_6 = c(sp_0_6, 100 + (sp_2 - 1)*10)
}
ind_5_6 = sample(seq(1, length(sp_0_5)), num_ind, replace = FALSE)
sp_0_5 = sp_0_5[ind_5_6]
sp_0_6 = sp_0_6[ind_5_6]

dt = data.frame(stronger_instrument = c(sp_0_8, sp_0_7, sp_0_6), 
                weaker_instrument = c(sp_0_3, sp_0_4, sp_0_5), 
                pair = c(rep('Pair 3', length(sp_0_8)), rep('Pair 2', length(sp_0_7)), rep('Pair 1', length(sp_0_6))))
g2 = ggplot(dt, aes(stronger_instrument, weaker_instrument, colour = pair)) + 
  geom_point() + 
  scale_color_manual(breaks = c("Pair 3", "Pair 2", "Pair 1"),
                     values=c("dark red", "dark blue", "dark green")) +
  labs(x = 'stronger instrument', y = 'weaker instrument') + 
  theme(axis.title = element_text(size=16)) + 
  geom_abline(intercept = 0, slope = 64/9, size = 1.5, colour = 'light green', alpha = 0.3) +
  geom_abline(intercept = 0, slope = 49/16, size = 1.5, colour = 'light blue', alpha = 0.3) +
  geom_abline(intercept = 0, slope = 36/25, size = 1.5, colour = 'pink', alpha = 0.3) + xlim(c(200, 6000)) + ylim(c(200, 6000))





