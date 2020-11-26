## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# loading the Rdata
load("data.RData")

# create the variable I(t)

I = covid[,2][9:220]
T = length(I)
r = rep(0,T)
for (t in 2:length(I)){
  r[t] = log(I[t]/I[t-1])
}

# plot r(t)
plot(r[2:T],type = "l",main="plot of r(t)",
     xlab = "time t",ylab = "r(t)")



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# create the difference variable 

del_r = rep(0,T)
for(t in 3:T){
  del_r[t] = r[t] - r[t-1]
}

# plot del_r(t)
plot(del_r[3:T],type = "l",main="plot of delta_r(t)",
     xlab = "time t",ylab = "delta_r(t)", ylim=c(-0.4, 0.4))



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# functions for estimating the mean and covariance functions

mean_estimator = function(stoc.process,timepoints){
  return(sum(stoc.process)/timepoints)
}

cov_estimator = function(stoc.process, timepoints, h){
  if (h >= T){
    print("Give value of h less then total number of time points.")
  }else{
    total_sum = 0
    for (t in 1:(timepoints - h)){
      within_sum_expression = (stoc.process[t]-mean_estimator(stoc.process, timepoints)) *
                              (stoc.process[t + h] - 
                               mean_estimator(stoc.process, timepoints))
      total_sum = total_sum + within_sum_expression
    }
    return(total_sum / timepoints)
  }
}

print(paste0("Estimated mean: ", mean_estimator(del_r[3:T], T - 2), "."))

# plot of the data for h from 1 to 30
cov_lag_h = rep(NA,T - 3)
for (h in 1:(T - 3)){
  cov_lag_h[h] = cov_estimator(del_r[3:T], T - 2, h)
}
plot(cov_lag_h[1:30],type = "l", main="plot of covariance estimators with lag",
     xlab = "lag h", ylab = "covariance function")



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# confidence interval 

estimated_mean = mean_estimator(del_r[3:T], T - 2)
estimated_variance = ((T - 2) * cov_estimator(del_r[3:T], T - 2, 0) + 
                      (2 * sum(cov_lag_h))) / (T - 2)^2
factor = qnorm(0.975, 0, 1) * sqrt(estimated_variance)
print(paste0("The 95% confidence interval is W = (",
            (estimated_mean - factor), ", ", (estimated_mean + factor),")."))



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)

# load 'mnormt' library
library(mnormt)

n_observations = 200

mean_vector = matrix(0, nrow = 1, ncol = n_observations)

covariance_matrix = matrix(0, nrow = n_observations, ncol = n_observations)

for(i in 1:n_observations) {
    for(j in 1:n_observations) {
        covariance_matrix[i, j] = exp(-abs(i - j))
    }
}

# Generate a random random vector with multivariate Normal distribution
generated_process = rmnorm(n = 1, mean = mean_vector, varcov = covariance_matrix)

# Plot the generated process
par(mar = c(5.1, 4.1, 3.1, 2.1))
plot(generated_process, type = 'l', main = 'Simulated Process', xlab = 'Time', 
     ylab = 'Observed Values')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
mean_estimator = function(sequence) {
    return(mean(sequence))
}

covariance_estimator = function(sequence, lag, est_mean) {
    sum = 0
    for(t in 1:(length(sequence) - lag)) {
        sum = sum + ((sequence[t] - est_mean) * (sequence[t + lag] - est_mean))
    }
    return(sum / length(sequence))
}


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
est_mean = mean_estimator(generated_process)
print(paste0('The true mean is \'0\' and the estimated mean is \'', est_mean, '\'.'))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
tru_covariance = matrix(0, nrow = 1, ncol = n_observations)
est_covariance = matrix(0, nrow = 1, ncol = n_observations)

for(i in 0:(n_observations - 1)) {
    tru_covariance[1, i + 1] = exp(-abs(i))
    est_covariance[1, i + 1] = covariance_estimator(generated_process, i, est_mean)
}


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot the real and estimated covariance functions
par(mar = c(5.1, 4.1, 3.1, 2.1))
plot(tru_covariance[1, 0:10], type = 'l', main = 'REAL vs ESTIMATED covariance functions',
     xlab = 'Lag', ylab = 'Covariance', ylim = c(-0.20, 1))
lines(est_covariance[1, 0:10], col = 'red')
legend('topright', legend = c('Real covariace', 'Estimated covariace'), 
       col = c('black', 'red'), pch = 19)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
mean_vector2 = matrix(0, nrow = 1, ncol = n_observations)

for(i in 1:n_observations) {
    mean_vector2[1, i] = sqrt(i)
}

generated_process2 = rmnorm(n = 1, mean = mean_vector2, varcov = covariance_matrix)

# Plot the generated process
par(mar = c(5.1, 4.1, 3.1, 2.1))
plot(generated_process2, type = 'l', main = 'Simulated Process 2', 
     xlab = 'Time', ylab = 'Observed Values')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
est_mean2 = mean_estimator(generated_process2)
est_mean2 = matrix(est_mean2, nrow = 1, ncol = n_observations)

par(mar = c(5.1, 4.1, 3.1, 2.1))
plot(mean_vector2[1, ], type = 'l', main = 'REAL vs ESTIMATED mean functions (2)',
     xlab = 'Time', ylab = 'Means', ylim = c(0, 14))
lines(est_mean2[1, ], col = 'red')
legend('topleft', legend = c('Real mean', 'Estimated mean'), 
       col = c('black', 'red'), pch = 19)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
tru_covariance2 = matrix(0, nrow = 1, ncol = n_observations)
est_covariance2 = matrix(0, nrow = 1, ncol = n_observations)

for(i in 0:(n_observations - 1)) {
    tru_covariance2[1, i + 1] = exp(-abs(i))
    est_covariance2[1, i + 1] = covariance_estimator(generated_process2, i, 
                                                     est_mean2[1, 1])
}

par(mar = c(5.1, 4.1, 3.1, 2.1))
plot(tru_covariance2[1, ], type = 'l', 
     main = 'REAL vs ESTIMATED covariance functions (2)',
     xlab = 'Lag', ylab = 'Covariance', ylim = c(-4, 4))
lines(est_covariance2[1, ], col = 'red')
legend('topright', legend = c('Real covariace', 'Estimated covariace'), 
       col = c('black', 'red'), pch = 19)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
log_return = matrix(NA, nrow = 1, ncol = length(sp))

for(i in 2:length(sp)) {
    log_return[1, i] = log(sp[i] / sp[i - 1])
}

# Plot the original data
par(mar = c(5.6, 4.1, 3.1, 2.1))
plot(sp, type = 'l', main = 'S&P 500 stock market index', 
     xlab = 'Day', ylab = 'Index value')


# Plot the transformed data
plot(log_return[1, ], type = 'l', 
     main = 'Log-returns of S&P 500 stock market index', 
     xlab = 'Day', ylab = '"Log-return" values', col = 'red')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Predicting the three next values of our process
new_length = length(sp) + 2

mean_vector_sp = matrix(0, nrow = 1, ncol = new_length)
covariance_matrix_sp = matrix(0, nrow = new_length, ncol = new_length)

sigma_2 = (0.02**2)
a = 0.4

for(i in 1:new_length) {
    for(j in 1:new_length) {
        covariance_matrix_sp[i, j] = (sigma_2 / (1 - a**2)) * ((-a)**(abs(i - j)))
    }
}

# Compute partial quantities of interest
x_a = log_return[1, 2:length(log_return[1, ])]

mean_a = mean_vector_sp[1, 1:(length(sp) - 1)]
mean_b = mean_vector_sp[1, length(sp):length(mean_vector_sp[1, ])]

sigma_aa = covariance_matrix_sp[1:(length(sp) - 1), 1:(length(sp) - 1)]
sigma_ab = covariance_matrix_sp[1:(length(sp) - 1), 
                                length(sp):length(covariance_matrix_sp[1, ])]
sigma_ba = covariance_matrix_sp[length(sp):length(covariance_matrix_sp[, 1]), 
                                1:(length(sp) - 1)]
sigma_bb = covariance_matrix_sp[length(sp):length(covariance_matrix_sp[, 1]), 
                                length(sp):length(covariance_matrix_sp[1, ])]


# Compute the conditional mean convariance values
mean_b_given_a  = mean_b + ((sigma_ba %*% solve(sigma_aa)) %*% (x_a - mean_a))
sigma_b_given_a = sigma_bb - (sigma_ba %*% solve(sigma_aa) %*% sigma_ab)

# Print final answer
for(i in 1:(new_length - (length(sp) - 1))) {
    print(paste0('R_', length(sp) + i, ' has mean \'',
                 format(mean_b_given_a[i], width = 10, digits = 5), 
                 '\' and variance \'', format(sigma_b_given_a[i, i], 
                                              width = 9, digits = 4), '\'.'))
}


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
compute_fourier_transform = function(process, freq) {
  return(sum(c(process) * exp(-1 * complex(real = 0, imaginary = 1) * 2 * pi * freq 
                              * c(0:(length(c(process)) - 1)))))
}

w = seq(from = -0.49, to = 0.50, by = 0.01)
periodogram = spectral_df = matrix(0, nrow = 1, ncol = length(w))

for(i in 1:(length(w))) {
  periodogram[, i] = (1 / length(generated_process)) * 
                     (Mod(compute_fourier_transform(generated_process, w[i])) ** 2)
  spectral_df[, i] = 2 / (1 + (2 * pi * w[i]) ** 2)
}

par(mar = c(5.1, 4.1, 3.1, 2.1))
plot(w, periodogram[1, ], type = 'l', main = 'PERIODOGRAM vs SPECTRAL DENSITY FUNCTION',
     xlab = 'Frequency', ylab = 'Power', xlim = c(-0.5, 0.5), ylim = c(0, 6))
lines(w, spectral_df[1, ], col = 'red')
legend('topright', legend = c('Periodogram', 'Spectral Density Function'), 
       col = c('black', 'red'), pch = 19)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
n_observations = 1000

mean_vector = matrix(0, nrow = 1, ncol = n_observations)
covariance_matrix = matrix(0, nrow = n_observations, ncol = n_observations)

for(i in 1:n_observations) {
  for(j in 1:n_observations) {
    covariance_matrix[i, j] = exp(-abs(i - j))
  }
}
generated_process3 = rmnorm(n = 1, mean = mean_vector, varcov = covariance_matrix)

periodogram3 = matrix(0, nrow = 1, ncol = length(w))

for(i in 1:(length(w))) {
  periodogram3[, i] = (1 / length(generated_process3)) * 
                      (Mod(compute_fourier_transform(generated_process3, w[i])) ** 2)
}

par(mar = c(5.1, 4.1, 3.1, 2.1))
plot(w, periodogram3[1, ], type = 'l', main = 'PERIODOGRAM vs SPECTRAL DENSITY FUNCTION (2)',
     xlab = 'Frequency', ylab = 'Power', xlim = c(-0.5, 0.5), ylim = c(0, 6))
lines(w, spectral_df[1, ], col = 'red')
legend('topright', legend = c('Periodogram', 'Spectral Density Function'), 
       col = c('black', 'red'), pch = 19)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------

compute_fourier_transform_v2 = function(process, freq, d) {
  return(sum(c(process) * exp(-1 * complex(real = 0, imaginary = 1) * 2 * pi * freq 
                              * (d * c(0:(length(c(process)) - 1))))))
}

w2 = seq(from = 0, to = 49.99, by = 0.01)

periodogram_eeg = matrix(0, nrow = 1, ncol = length(w2))

for(i in 1:(length(w2))) {
  periodogram_eeg[, i] = (1/100) * (1 / length(eeg)) * 
                      (Mod(compute_fourier_transform_v2(eeg, w2[i], 1/100)) ** 2)
}

par(mar = c(5.1, 4.1, 3.1, 2.1))
plot(w2, periodogram_eeg[1, ], type = 'l', main = 'EEG PERIODOGRAM',
     xlab = 'Frequency', ylab = 'Power', xlim = c(0, 50))



## ----message=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
library(bspec)
welsh = welchPSD(ts(eeg, frequency = 100), seglength = 4)
plot(welsh$frequency, welsh$power, type = 'l', main = 'EEG WELSH\'S PERIODOGRAM',
     xlab = 'Frequency', ylab = 'Power', xlim = c(0, 50))

compute_approx_area = function(low, upper) {
  low_position   = match(low,   welsh$frequency)
  upper_position = match(upper, welsh$frequency)
  area = 0
  for(i in (low_position + 1):upper_position) {
    area = area + ((welsh$frequency[i] - welsh$frequency[i - 1]) * 
                   welsh$power[i])
  }
  return(area)
}  
print(paste0('The average band power of the \'delta\' band was ', 
             round(compute_approx_area(0.5, 4), 3), '.'))
print(paste0('The average band power of the \'theta\' band was ', 
             round(compute_approx_area(4,   8), 3), '.'))
print(paste0('The average band power of the \'alpha\' band was ', 
             round(compute_approx_area(8,  12), 3), '.'))
print(paste0('The average band power of the \'beta\'  band was ', 
             round(compute_approx_area(12, 30), 3), '.'))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
impulse_response = function(a, b, t) {
  return((sin(2 * pi * t * b) - sin(2 * pi * t * a)) / (pi * t))
}

frequency_filter = function(process, d, a, b) {
  new_process = rep(0, length(process))
  
  for(t in 1:(length(process))) {
    sum = 0
    for(i in 0:(length(process) - 1)) {
      if((t - i + 1) > 0) {
        if(i == 0) {
          sum = sum + (2 * (b - a) * process[t - i + 1] * d)
        } else {
          sum = sum + impulse_response(a, b, d * i) * process[t - i + 1] * d
        }
      }
    }
    new_process[t] = sum
  }
  
  return(new_process)
}

delta = frequency_filter(eeg, 1/100, 0.5, 4)
theta = frequency_filter(eeg, 1/100, 4,   8)
alpha = frequency_filter(eeg, 1/100, 8,  12)
beta  = frequency_filter(eeg, 1/100, 12, 30)

x_axis = 1:(length(eeg)) * (1/100)
plot(x_axis, delta, type = 'l', main = 'Delta Band',
   xlab = 'Time', ylab = 'Value')

plot(x_axis, theta, type = 'l', main = 'Theta Band',
   xlab = 'Time', ylab = 'Value')

plot(x_axis, alpha, type = 'l', main = 'Alpha Band',
   xlab = 'Time', ylab = 'Value')

plot(x_axis, beta,  type = 'l', main = 'Beta  Band',
   xlab = 'Time', ylab = 'Value')




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
customer = samba$customer
processing_time = samba$processing_time
T = length(customer)
interarrival_time = rep(NA,T)
interarrival_time[1] = customer[1]

for (i in 2:T){
  interarrival_time[i] = customer[i]-customer[i-1]
}

# we know inter arrival time follows exponential 

est_lambda = 1 / (mean(interarrival_time))

# processing time also follows exponential 
est_mu = 1 / (mean(processing_time))

# estimates of lambda and mu

print(paste0("The estimates of mu and lambda are: ",
            est_lambda, " and ", est_mu, ", respectively."))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# simulate 5000 points 

n = 5000
set.seed(233784)
interarrival_time_simulated = rexp(n , rate = est_lambda)
processing_time_simulated = rexp(n , rate = est_mu)

customer_simulated = rep(NA,n)
customer_simulated[1] = interarrival_time_simulated[1]

for(i in 2:n){
  customer_simulated[i] = customer_simulated[i-1] + interarrival_time_simulated[i]
}
waiting_time = rep(NA,n)
waiting_time[1] = 0

for(i in 2:n){
  waiting_time[i] = max(0, 
                        (customer_simulated[i-1]+
                           waiting_time[i-1] +
                          processing_time_simulated[i-1]-
                          customer_simulated[i]))
}
average_waiting_time = mean(waiting_time)
print(paste("Average waiting time from simulated data is : ", average_waiting_time))
h <- hist(waiting_time, breaks = 20 , col="brown", probability = TRUE, 
          xlab="waiting time of simlulated data",
          main="Histogram with a smooth curve")
lines(h$mids,h$density,col="red",lwd=4)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Waiting time for a window of 4
waiting_time_revised = rep(NA,n)
waiting_time_revised[1] = 0
waiting_time_revised[2] = 0
waiting_time_revised[3] = 0
waiting_time_revised[4] = 0

index = c(1,2,3,4)

for(i in 5:n){
  t1 = index[1]
  t2 = index[2]
  t3 = index[3]
  t4 = index[4]
  k1 = customer_simulated[t1] + processing_time_simulated[t1] + waiting_time_revised[t1]
  k2 = customer_simulated[t2] + processing_time_simulated[t2] + waiting_time_revised[t2]
  k3 = customer_simulated[t3] + processing_time_simulated[t3] + waiting_time_revised[t3]
  k4 = customer_simulated[t4] + processing_time_simulated[t4] + waiting_time_revised[t4]
  f = c(k1,k2,k3,k4)
  m = which.min(f)
  waiting_time_revised[i] = max(0,(min(f) - customer_simulated[i]))
  index[m] = i
}

average_waiting_time_revised = mean(waiting_time_revised)
print(paste("Average waiting time from simulated data is : ", 
            average_waiting_time_revised))
h <- hist(waiting_time_revised, breaks = 20 , col="brown", probability = TRUE,
          xlab="waiting time of simlulated data",
          main="Histogram of waiting times with window of 4 with a smooth curve")
lines(h$mids,h$density,col="red",lwd=4)

