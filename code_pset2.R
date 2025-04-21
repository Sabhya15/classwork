rm(list = ls())

library(readxl)
library(lubridate)
library(dplyr)
library(lattice)
library(ggplot2)
library(ggthemes)
library(zoo)
library(stargazer)

# Set the working directory to the location of your Excel file
root <- "//Users//sg137//BOSTON UNIVERSITY Dropbox//Sabhya Gupta//Apps//GoodNotes 5//Goodnotes//PhD program//Second Year//Spring 2025//IO//Assignments//Assignment_2"
setwd(root)
dt <- read_excel(paste0(root, "//dvddataEC732.xls"))

# cleaning the data
dt$coax <- ifelse(is.na(dt$coax), 0, dt$coax)
dt$mp3 <- ifelse(is.na(dt$mp3), 0, dt$mp3)

# market share variables
market_size = 100000000
dt$s_jt = (dt$sales)/market_size
dt <- dt %>%
  group_by(month) %>%
  mutate(s_t = sum(s_jt), s_0t = 1-sum(s_t))

dt$ln_s_jt = log(dt$s_jt)
dt$ln_st = log(dt$s_t)

dt$s_jt_s_t = dt$s_jt/dt$s_t
dt$ln_s_jt_s_t = dt$ln_s_jt - dt$ln_st
dt$ln_price = log(dt$price)

# Part a

# estimate logit to get the beta values
logit_model <- lm(ln_s_jt_s_t ~ coax + mp3 + ln_price + factor(month) -1  , data = dt)
summary(logit_model)

# export regression results
stargazer(logit_model, title = "Logit Regression Results", 
          dep.var.labels = c("DVD Player Sales"), keep = c("coax", "mp3", "ln_price"),
          covariate.labels = c("Coax Cable", "MP3 Output", "ln(price)"), 
          omit.stat = c("f", "ser", "adj.rsq") , column.sep.width = "1pt", 
          out = file.path("logit_regression_results.tex"), label = "reg",
           digits = 3,header = FALSE, no.space = TRUE, omit.table.layout = "n")



# Part b
state_spaces <- seq(from = -1, to = 3, length.out = 20)
r_values <- tibble(-summary(logit_model)$coefficients[4:22,1])
colnames(r_values) <- c("r_t")
r_values$r_t_1 <- lead(r_values$r_t)

r_regression <- lm(r_t_1 ~ r_t, data=r_values)
summary(r_regression)

alpha_0 <- summary(r_regression)$coefficients[1]
alpha_1 <- summary(r_regression)$coefficients[2]
mu_sigma <- summary(r_regression)$sigma

stargazer(r_regression, title = "AR1 Process for value of purchase", 
          dep.var.labels = c("Month Fixed Effect"), keep = c("r_t"),
          covariate.labels = c("Lagged Month Fixed Effect"), 
          omit.stat = c("f", "ser", "adj.rsq") , column.sep.width = "1pt", 
          out = file.path("ar1.tex"), label = "ar1",
          digits = 3,header = FALSE, no.space = TRUE, omit.table.layout = "n")




# transition matrix
rtauchen_func = function(grid, alpha_0, alpha_1, mu_sigma){
  # grid: vector of state space
  # alpha_0: intercept of the regression
  # alpha_1: slope of the regression
  # mu_sigma: standard deviation of the regression
  
  n <- length(grid)
  delta <- grid[2] - grid[1]
  
  transition_matrix <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (j == 1) {
        transition_matrix[i, j] <- pnorm((grid[j] - alpha_0 - alpha_1 * grid[i] + delta / 2) / mu_sigma)
      } else if (j == n) {
        transition_matrix[i, j] <- pnorm((grid[j] - alpha_0 - alpha_1 * grid[i] - delta / 2) / mu_sigma, lower.tail = FALSE)
      } else {
        transition_matrix[i, j] <- pnorm((grid[j] - alpha_0 - alpha_1 * grid[i] + delta / 2) / mu_sigma) -
          pnorm((grid[j] - alpha_0 - alpha_1 * grid[i] - delta / 2) / mu_sigma)
      }
    }
  }
  
  return(transition_matrix)

} # tauchen method

transition_matrix <- rtauchen_func(state_spaces, alpha_0, alpha_1, mu_sigma)

# graph matrix 
rgb.palette <- colorRampPalette(c("white", "red"), space = "rgb")

png(filename="transition_matrix.png", type="cairo",
    units="in", 
    width=5, 
    height=4, 
    pointsize=12, 
    res=96)
levelplot(transition_matrix, xlab="Bin of this period r_t", ylab="Bin of next period r_t",  col.regions=rgb.palette(120))

dev.off()


# Part c
beta = 0.99
c = 0
n = length(state_spaces)

r_midpoint = rep(0, n)
for (i in 1:n-1){
  r_midpoint[i] = (state_spaces[i+1] + state_spaces[i])/2
}
r_midpoint[n] = 3

# value function iteration
v_guess = rep(0, length(state_spaces))
V = rep(0, length(state_spaces))
tol = 1e-6
max_iter = 1000

vold = v_guess
distance = 10
iter = 0

VIF <- function(c, V, V_new, distance, iter, beta){
  while ((distance > tol) | (iter < max_iter)) {
    V = vold
    V_new = rep(0,n)
    
    EV = transition_matrix %*% V # Expected value at different time periods 
    max_possible_value <- pmax(r_midpoint, c + beta * EV, na.rm = TRUE) # best choice at each r value
    
    V_new = max_possible_value + log(exp(r_midpoint - max_possible_value) + exp(c + beta * EV - max_possible_value))
    
    distance = max(abs(V_new - V))
    vold = V_new
    iter = iter + 1
  }
  
  V = vold
  return(V)
}

V_dataset <- as.data.frame(r_midpoint)
V_dataset$V = VIF(c, V, V_new, distance, iter, beta)

# plot 

ggplot(V_dataset, aes(x = r_midpoint, y = V)) + geom_line() + geom_point() +
  labs(x = "r_t", y = "Value function at c = 0") +
  theme_light() +
  scale_x_continuous(breaks = seq(-1, 3, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 6, by = 0.25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("V_function_c0.png", width = 8, height = 6, dpi = 300)

# Part d

# s_T+1  = people who never purchased (or everyone else buys in T+1)
s_T_1 <- 1 - sum(unique(dt$s_t))

s_t = dt %>%
  group_by(month) %>%
  summarise(s_t = sum(s_jt)) %>%
  select(s_t)

actual_s_t <- c(s_t$s_t, s_T_1)

s_predicted_t <- rep(0, length(actual_s_t))

# predicted s_t function
predicted_s_t <- function(c, V, V_new, distance, iter, beta){
 
  V_dataset_c <- as.data.frame(r_midpoint)
  V_dataset_c$V_c = VIF(c, V, V_new, distance, iter, beta)
  
  remainingshare = 1
  s_predicted_t <- rep(0, (length(s_t$s_t) + 1))
  # compute shares
  for (i in 1:length(s_t$s_t)){
    r_value = r_values$r_t[i]
    # closest midpoint
    r_midpoint_index = which(abs(r_midpoint-r_value ) == min(abs(r_midpoint-r_value )))

    W_t = c + beta * sum(transition_matrix[r_midpoint_index,] * V_dataset_c$V_c)
    max_possible_values <- max(r_midpoint[r_midpoint_index], W_t, na.rm = TRUE) 
    
    
    logit_prob_buy = exp(r_midpoint[r_midpoint_index] - max_possible_values)/(exp(W_t - max_possible_values) + exp(r_value - max_possible_values))
    s_predicted_t[i] = remainingshare * logit_prob_buy
    remainingshare =  remainingshare - logit_prob_buy
  }
  
  s_predicted_t[length(s_t$s_t) + 1] = remainingshare
  

  # return the predicted s_t
  return(s_predicted_t)
}

s_t_hat <- predicted_s_t(c, V, V_new, distance, iter, beta)

max_likelihood <- function(c, actual_s_t){
  
  # calculate the predicted shares
  
  s_t_hat <- predicted_s_t(c, V, V_new, distance, iter, beta)
  
  # calculate the log-likelihood
  log_likelihood <- sum(log(s_t_hat) * actual_s_t)
  
  # return the negative log-likelihood
  return(-log_likelihood)
  
}

# graph likelihood
c_values <- as.list(seq(0, 0.2, length.out = 100))
likelihood_c <- -unlist(lapply(c_values, FUN = max_likelihood,actual_s_t  = actual_s_t))
c_values = unlist(c_values)

ggplot() + geom_line(aes(x=c_values, y = likelihood_c), color = "blue") +
  labs(x = "c", y = "Log-likelihood") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("log_likelihood_c.png", width = 8, height = 6, dpi = 300)


# optimize the log-likelihood function
optim_result <- optim(par = 0, fn = max_likelihood, actual_s_t = actual_s_t, method = "L-BFGS-B", lower = 0, upper = 1)

optimal_c = optim_result$par
likelihood_optimal <- max_likelihood(optimal_c, actual_s_t)

# graph the value function at optimal c
V_dataset_optimal <- as.data.frame(r_midpoint)
V_dataset_optimal$V_c_optimal = VIF(optimal_c, V, V_new, distance, iter, beta)

ggplot(V_dataset_optimal, aes(x = r_midpoint, y = V_c_optimal)) + geom_line() + geom_point() +
  labs(x = "r_t", y = "Value function at optimal c value") +
  theme_light() +
  scale_x_continuous(breaks = seq(-1, 3, by = 0.5)) +
  scale_y_continuous(breaks = seq(9, 10, by = .005)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("V_function_coptimal.png", width = 8, height = 6, dpi = 300)

# graph the predicted and market shares
s_t_graph <- as.data.frame(actual_s_t)
s_t_graph$predicted_s_t_optimal <- predicted_s_t(optimal_c, V, V_new, distance, iter, beta)
s_t_graph$month <-c(unique(dt$month), 465)
s_t_graph <- s_t_graph[1:19,]
s_t_graph$month_text <- seq(as.Date("1997/03/1"), by = "month", length.out = 19)

ggplot() +
  geom_line(data = s_t_graph, aes(x = month_text, y = actual_s_t, color = "Actual")) + 
  geom_line(data = s_t_graph, aes(x = month_text, y = predicted_s_t_optimal, color = "predicted")) + 
  labs(x = "Month", y = "Market Share", color = "") +  scale_y_continuous(breaks = seq(0, 0.001, by = .00005))+
   scale_x_date(breaks = "1 month", date_labels = "%b %Y") + theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = c(Actual = "blue", predicted = "red"),
                    labels = c(Actual = "Actual market share", predicted = "Predicted market share")) 
                    

ggsave("market_shares_predicted_actual.png", width = 12, height = 6)
