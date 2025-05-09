rm(list=ls())

library(tidyverse)
library(readxl)
library(stargazer)
library(data.table)
library(fastDummies)
library(ggplot2)

root = "//Users//sg137//BOSTON UNIVERSITY Dropbox//Sabhya Gupta//Apps//GoodNotes 5//Goodnotes//PhD program//Second Year//Spring 2025//IO//Assignments//Final"
data = paste0(root, "//data")
code = paste0(root, "//code")
table = paste0(root, "//tables")
fig = paste0(root, "//figures")

setwd(root)


dt <- read_excel(paste0(data, "//airlineEntry(2).xlsx"))

######## Question 1 #####

probit = glm(entry ~ lnPax + cfNumEntrant, family = binomial(link = "probit"), 
             data = dt)
summary(probit)

ols = lm(entry ~ lnPax + cfNumEntrant, data = dt)
summary(ols)

stargazer(probit, ols, title = "Probit and OLS of entry on covariates", 
          dep.var.labels = c("Dummy on Entry"), keep = c("lnPax", "cfNumEntrant"),
          covariate.labels = c("Log(Passengers)", "Counterfactual Entrants"), 
          column.labels = c("Probit", "OLS"),omit.stat = c("f", "ser", "adj.rsq", "aic"), 
          column.sep.width = "1pt", out = file.path(paste0(table, "//IO_q1.tex")), 
          label = "probit",digits = 3,header = FALSE, no.space = TRUE, 
          omit.table.layout = "n", model.names = FALSE, model.numbers = FALSE)

######## Question 2 #####

beta_0 = -2
beta_1 = -0.1
beta_2 = 0.1
beta_3 = -0.5
beta_4 = -0.5
beta_5 = -0.5

beta = c(beta_0, beta_1, beta_2, beta_3, beta_4, beta_5)

# calculate P_im and max likelihood using 2 methods 

dt <- cbind(dt, dummy_cols(dt$cr)[2:ncol(dummy_cols(dt$cr))])
dt <- dt %>%
  rename(AA = .data_AA, 
         DL = .data_DL,
         UA = .data_UA, 
         WN = .data_WN)


# method 1
dt$p=  pnorm(beta_0 + beta_1*(dt$cfNumEntrant) + beta_2 * (dt$lnPax) + 
  beta_3 * (dt$AA) + beta_4 * (dt$DL) + beta_5 * dt$UA, mean = 0, sd = 1)

mean(dt$p)
sd(dt$p)


likelihood1 <- function(beta, data){
  
  # calculate p_im
  data$p=  pnorm(beta[1] + beta[2] *(data$cfNumEntrant) + beta[3]  * (data$lnPax) + 
                 beta[4]  * (data$AA) + beta[5]  * (data$DL) + beta[6]  * data$UA, mean = 0, sd = 1)
  
  # likelihood function
  likelihood_value <- sum(data$entry * log(data$p)  + (1 - data$entry) * log(1 - data$p))
  
  return(-likelihood_value)
}

optim_result <- optim(par = beta, fn = likelihood1, data = dt)
params_method1 = optim_result$par

# method 2
tol = 1e-6
max_iter = 100

dt = dt %>%
  group_by(apCode1, apCode2) %>%
  mutate(expected_n_im_all = sum(p))

dt$expected_n_im <- dt$expected_n_im_all - dt$p

p_new = pnorm(beta_0 + beta_1*(dt$expected_n_im) + beta_2 * (dt$lnPax) + 
                beta_3 * (dt$AA) + beta_4 * (dt$DL) + beta_5 * dt$UA, mean = 0, sd = 1)

distance = 100
iter = 0

solve_p_im = function(beta, tol, max_iter, data){
  while ((distance > tol) & (iter < max_iter)) {

    data$p_new = pnorm(beta[1] + beta[2]*(data$expected_n_im) + beta[3] * (data$lnPax) + 
                  beta[4] * (data$AA) + beta[5] * (data$DL) + beta[6] * data$UA, mean = 0, sd = 1)
    
    data = data %>%
      group_by(apCode1, apCode2) %>%
      mutate(expected_n_im_all_new = sum(p_new))
    
    
    data$expected_n_im_new <- data$expected_n_im_all_new - data$p_new
    
    distance = max(abs(data$p_new - data$p))
    
    data$p = data$p_new
    data$expected_n_im_all = data$expected_n_im_all_new
    data$expected_n_im = data$expected_n_im_new
    
    iter = iter + 1
    print(iter)
    print(distance)
    }

  data$p = data$p_new
  return(data)

}

tol = 1e-6
max_iter = 100
dt = solve_p_im(beta, tol, max_iter, dt)

mean(dt$p)
sd(dt$p)

# solve the likelihood 

likelihood2 <- function(beta,  tol, max_iter, data){

  data <- solve_p_im(beta, tol, max_iter, data)
 
  
  likelihood_value <- sum(data$entry * log(data$p)  + (1 - data$entry) * log(1 - data$p)) # likelihood function
  return(-likelihood_value)
  
}


optim_result <- optim(par = beta, fn = likelihood2, tol = tol, max_iter = max_iter, data = dt)
params_method2 = optim_result$par

######## Question 3 #####
# 3.a moment inq 
b1_support = seq(-2, 2, length.out = 500)
b0_more_than = rep(0, length(b1_support))
b0_less_than = rep(0, length(b1_support))


for (i in 1:length(b1_support)){
  b0_more_than[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 1)] + dt$lnPax[which(dt$entry == 1)])) / length(dt$entry[which(dt$entry == 1)])
  b0_less_than[i] = - sum((b1_support[i]  * dt$cfNumEntrant[which(dt$entry == 0)] + dt$lnPax[which(dt$entry == 0)])) / length(dt$entry[which(dt$entry == 0)])
  
}

moments = data.frame(b1_support = b1_support, b0_less_than = b0_less_than, b0_more_than =b0_more_than )

ggplot(data = moments, aes(x = b1_support)) + 
  geom_line(aes(y = b0_less_than, color = 'No entry')) +
  geom_line(aes(y = b0_more_than, color = 'Entry')) +
  geom_ribbon(aes(ymin = b0_more_than, ymax = b0_less_than), fill="purple", alpha=0.25) +
  scale_color_manual(values = c('No entry' = 'blue','Entry' = 'red')) +
  ylab("Beta_0") + xlab("Beta_1") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(filename = paste0(fig, "//moments_IO_q3_a.png"))

# 3.b instruments 

b0_more_than_z1 = rep(0, length(b1_support))
b0_less_than_z0 = rep(0, length(b1_support))
b0_more_than_lnpax1 =  rep(0, length(b1_support))
b0_less_than_lnpax0 = rep(0, length(b1_support))
b0_more_than_numentrant1 =  rep(0, length(b1_support))
b0_less_than_numentrant0 =  rep(0, length(b1_support))


for (i in 1:length(b1_support)){
  b0_more_than_z1[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 1)] + dt$lnPax[which(dt$entry == 1)])) / length(dt$entry[which(dt$entry == 1)])
  b0_less_than_z0[i] = - sum((b1_support[i]  * dt$cfNumEntrant[which(dt$entry == 0)] + dt$lnPax[which(dt$entry == 0)])) / length(dt$entry[which(dt$entry == 0)])
  
  b0_more_than_lnpax1[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 1)] * dt$lnPax2015[which(dt$entry == 1)] + dt$lnPax[which(dt$entry == 1)]  * dt$lnPax2015[which(dt$entry == 1)])) / sum(dt$lnPax2015[which(dt$entry == 1)])
  b0_less_than_lnpax0[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 0)] * dt$lnPax2015[which(dt$entry == 0)] + dt$lnPax[which(dt$entry == 0)]  * dt$lnPax2015[which(dt$entry == 0)])) / sum(dt$lnPax2015[which(dt$entry == 0)])
  
  b0_more_than_numentrant1[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 1)] * dt$numEntrant2015[which(dt$entry == 1)] + dt$lnPax[which(dt$entry == 1)]  * dt$numEntrant2015[which(dt$entry == 1)])) / sum(dt$numEntrant2015[which(dt$entry == 1)])
  b0_less_than_numentrant0[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 0)] * dt$numEntrant2015[which(dt$entry == 0)] + dt$lnPax[which(dt$entry == 0)]  * dt$numEntrant2015[which(dt$entry == 0)])) / sum(dt$numEntrant2015[which(dt$entry == 0)])
  
}


moments_inst = data.frame(b1_support = b1_support, 
                          b0_more_than_z1 = b0_more_than_z1, 
                          b0_less_than_z0 = b0_less_than_z0, 
                          b0_more_than_lnpax1 = b0_more_than_lnpax1, 
                          b0_less_than_lnpax0 = b0_less_than_lnpax0, 
                          b0_more_than_numentrant1 = b0_more_than_numentrant1, 
                          b0_less_than_numentrant0 = b0_less_than_numentrant0)



ggplot(data = moments_inst, aes(x = b1_support)) + 
  geom_line(aes(y = b0_more_than_z1, color = 'Entry')) +
  geom_line(aes(y = b0_less_than_z0, color = 'No Entry')) +
  geom_line(aes(y = b0_more_than_lnpax1, color = 'Entry')) +
  geom_line(aes(y = b0_less_than_lnpax0, color = 'No Entry')) +
  geom_line(aes(y = b0_more_than_numentrant1, color = 'Entry')) +
  geom_line(aes(y = b0_less_than_numentrant0, color = 'No Entry'))+ 
  geom_ribbon(aes(ymin = b0_more_than_z1, ymax = b0_less_than_z0), fill="yellow", alpha=0.3) + 
  geom_ribbon(aes(ymin = b0_more_than_lnpax1, ymax = b0_less_than_lnpax0), fill="red", alpha=0.3) + 
  geom_ribbon(aes(ymin = b0_more_than_numentrant1, ymax = b0_less_than_numentrant0), fill="red", alpha=0.3) +
  scale_color_manual(values = c('No Entry' = 'blue','Entry' = 'red')) +
  ylab("Beta_0") + xlab("Beta_1") +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave(filename = paste0(fig, "//moments_inst_IO_q3_b.png"))

# 3.c moment inequality

b0_more_than_numentrant1_1 = rep(0, length(b1_support))
b0_less_than_numentrant0_1 = rep(0, length(b1_support))

b0_more_than_numentrant1_2 = rep(0, length(b1_support))
b0_less_than_numentrant0_2 = rep(0, length(b1_support))

b0_more_than_numentrant1_3 = rep(0, length(b1_support))
b0_less_than_numentrant0_3 = rep(0, length(b1_support))

b0_more_than_numentrant1_4 = rep(0, length(b1_support))
b0_less_than_numentrant0_4 = rep(0, length(b1_support))

for (i in 1:length(b1_support)){
  b0_more_than_numentrant1_1[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 1 & dt$numEntrant2015 == 1)] * dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 1)] + dt$lnPax[which(dt$entry == 1 & dt$numEntrant2015 == 1)]  * dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 1)])) / sum(dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 1)])
  b0_less_than_numentrant0_1[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 0 & dt$numEntrant2015 == 1)] * dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 1)] + dt$lnPax[which(dt$entry == 0 & dt$numEntrant2015 == 1)]  * dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 1)])) / sum(dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 1)])
  
  b0_more_than_numentrant1_2[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 1 & dt$numEntrant2015 == 2)] * dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 2)] + dt$lnPax[which(dt$entry == 1 & dt$numEntrant2015 == 2)]  * dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 2)])) / sum(dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 2)])
  b0_less_than_numentrant0_2[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 0 & dt$numEntrant2015 == 2)] * dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 2)] + dt$lnPax[which(dt$entry == 0 & dt$numEntrant2015 == 2)]  * dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 2)])) / sum(dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 2)])
  
  b0_more_than_numentrant1_3[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 1 & dt$numEntrant2015 == 3)] * dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 3)] + dt$lnPax[which(dt$entry == 1 & dt$numEntrant2015 == 3)]  * dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 3)])) / sum(dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 3)])
  b0_less_than_numentrant0_3[i] = - sum((b1_support[i] *dt$cfNumEntrant[which(dt$entry == 0 & dt$numEntrant2015 == 3)] * dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 3)] + dt$lnPax[which(dt$entry == 0 & dt$numEntrant2015 == 3)]  * dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 3)])) / sum(dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 3)])
  
  b0_more_than_numentrant1_4[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 1 & dt$numEntrant2015 == 4)] * dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 4)] + dt$lnPax[which(dt$entry == 1 & dt$numEntrant2015 == 4)]  * dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 4)])) / sum(dt$numEntrant2015[which(dt$entry == 1 & dt$numEntrant2015 == 4)])
  b0_less_than_numentrant0_4[i] = - sum((b1_support[i] * dt$cfNumEntrant[which(dt$entry == 0 & dt$numEntrant2015 == 4)] * dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 4)] + dt$lnPax[which(dt$entry == 0 & dt$numEntrant2015 == 4)]  * dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 4)])) / sum(dt$numEntrant2015[which(dt$entry == 0 & dt$numEntrant2015 == 4)])
  
}

moments_inst_2 = data.frame(b1_support = b1_support, 
                          b0_more_than_z1 = b0_more_than_z1, 
                          b0_less_than_z0 = b0_less_than_z0, 
                          b0_more_than_lnpax1 = b0_more_than_lnpax1, 
                          b0_less_than_lnpax0 = b0_less_than_lnpax0, 
                          b0_more_than_numentrant1_1 = b0_more_than_numentrant1_1,
                          b0_less_than_numentrant0_1 = b0_less_than_numentrant0_1,
                          b0_more_than_numentrant1_2 = b0_more_than_numentrant1_2,
                          b0_less_than_numentrant0_2 = b0_less_than_numentrant0_2,
                          b0_more_than_numentrant1_3 = b0_more_than_numentrant1_3,
                          b0_less_than_numentrant0_3 = b0_less_than_numentrant0_3,
                          b0_more_than_numentrant1_4 = b0_more_than_numentrant1_4,
                          b0_less_than_numentrant0_4 = b0_less_than_numentrant0_4)

ggplot(data = moments_inst_2, aes(x = b1_support)) + 
  geom_line(aes(y = b0_more_than_z1, color = 'Entry')) +
  geom_line(aes(y = b0_less_than_z0, color = 'No Entry')) +
  geom_line(aes(y = b0_more_than_lnpax1, color = 'Entry')) +
  geom_line(aes(y = b0_less_than_lnpax0, color = 'No Entry')) +
  geom_line(aes(y = b0_more_than_numentrant1_1, color = 'Entry')) +
  geom_line(aes(y = b0_less_than_numentrant0_1, color = 'No Entry')) +
  geom_line(aes(y = b0_more_than_numentrant1_2, color = 'Entry')) +
  geom_line(aes(y = b0_less_than_numentrant0_2, color = 'No Entry'))+ 
  geom_line(aes(y = b0_more_than_numentrant1_3, color = 'Entry')) +
  geom_line(aes(y = b0_less_than_numentrant0_3, color = 'No Entry'))+ 
  geom_line(aes(y = b0_more_than_numentrant1_4, color = 'Entry')) +
  geom_line(aes(y = b0_less_than_numentrant0_4, color = 'No Entry'))  +
  annotate(geom = "rect", xmin=-1.5, xmax=-.5, ymin= -26, ymax = -23.5, alpha=0.3, fill="red") + 
  ylab("Beta_0") + xlab("Beta_1") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(filename = paste0(fig, "//moments_inst_IO_q3_c.png"))
  

