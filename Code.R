  # 1. Import dataset
load("2536706.RData")
ls()
View(dat)


  # 2. Exploratory Analysis
library(ggplot2)
library(GGally)

print(head(dat))

treatment_colors <- c("red", "blue")

#pairs(~ rest_pulse + stimulated_pulse + bmi + treatment, data=dat, col = treatment_colors[as.numeric(as.character(dat$treatment)) + 1])
ggpairs(dat, mapping = ggplot2::aes(colour=treatment))

print(groups_count <- table(dat$treatment))

  
  # 3. Model by the clinician
# Using R functions
fit0 <- lm(stimulated_pulse ~ rest_pulse + treatment, data=dat)

summary(fit0)
par(mfrow=c(2,2))
plot(fit0)


  # 4. Model by the statistician
      fit1 <- glm((stimulated_pulse - rest_pulse) ~ bmi + treatment, data=dat, family = Gamma(link = "inverse"))
      summary(fit1)
      par(mfrow=c(2,2))
      plot(fit1)
      par(mfrow=c(1,1))
      
      fit2 <- glm((stimulated_pulse - rest_pulse) ~ bmi * treatment + bmi*treatment, data=dat, family = Gamma(link = "inverse"))
      summary(fit2)
      par(mfrow=c(2,2))
      plot(fit2)
      par(mfrow=c(1,1))
      
      
# Build X matrix
data <- dat[, c(3,4)]
X_bt <- data.matrix(data)

# Initialise beta vector parameters
epsilon = 1e-9
beta_4 <- rep(epsilon, 3)

# stimulated_pulse - rest_pulse
y = dat[, 2] - dat[, 1]

# IWLS
for (i in 1:250){
  eta_4 <- cbind(1,X_bt)%*%beta_4 # Estimated linear predictor
  mu_4 <- 1 / eta_4 # Estimated mean response
  z_4 <-  eta_4 + (mu_4 - y)/mu_4^2 # Form the adjusted dependent variable in terms of y and mu
  w_4 <-  mu_4^2 # Compute weights in terms of mu
  lmod <- lm(z_4 ~ X_bt, weights = w_4) # Regress z on x with weights w.
  beta_4 <- as.numeric(lmod$coeff) # New beta
}

print(beta_4)


  # 5.
names(which.max(table(dat$bmi))) # Get most recurrent bmi

eta_4 <- cbind(1,X_bt)%*%beta_4
mu_4 <- 1 / eta_4
w_4 <- as.vector(mu_4^2)
W_4 <-  diag(w_4)

print(mu_4[92] - mu_4[38]) # Difference between two individuals with same bmi but different group

R_4_0 = 1 / (eta_4[92] - 1.96*sqrt(t(X_bt[92, ]) %*% solve(t(X_bt) %*% W_4 %*% X_bt) %*% X_bt[92, ]))
L_4_0 = 1 / (eta_4[92] + 1.96*sqrt(t(X_bt[92, ]) %*% solve(t(X_bt) %*% W_4 %*% X_bt) %*% X_bt[92, ]))

R_4_1 = 1 / (eta_4[38] - 1.96*sqrt(t(X_bt[38, ]) %*% solve(t(X_bt) %*% W_4 %*% X_bt) %*% X_bt[38, ]))
L_4_1 = 1 / (eta_4[38] + 1.96*sqrt(t(X_bt[38, ]) %*% solve(t(X_bt) %*% W_4 %*% X_bt) %*% X_bt[38, ]))

print(paste('(',L_4_0 - L_4_1,',',R_4_0 - R_4_1,')')) # Confidence interval


# Quick computations of results
mean_placebo = mean(mu_4[dat$treatment == 1])
mean_treatment = mean(mu_4[dat$treatment == 0])
mean_difference = mean_treatment - mean_placebo


  # 6.
data <- dat[, c(3,4)]
# Add predictors
data$prod <- data$bmi * as.numeric(as.character(data$treatment))
X_aug <- data.matrix(data)

# Initialise beta, vector parameters
epsilon = 1e-9
beta_5 <- rep(epsilon, 4)

# stimulated_pulse - rest_pulse
y = dat[, 2] - dat[, 1]

# IWLS
for (i in 1:250){
  eta_5 <- cbind(1,X_aug)%*%beta_5 # Estimated linear predictor
  mu_5 <- 1 / eta_5 # Estimated mean response
  z_5 <-  eta_5 + (mu_5 - y)/mu_5^2 # Form the adjusted dependent variable in terms of y and mu
  w_5 <-  mu_5^2 # Compute weights in terms of mu
  lmod <- lm(z_5 ~ X_aug, weights = w_5) # Regress z on x with weights w.
  beta_5 <- as.numeric(lmod$coeff) # New beta
}

print(beta_5)

eta_5 <- cbind(1,X_aug)%*%beta_5
mu_5 <- 1 / eta_5
z_5 <-  eta_5 + (mu_5 - y)/mu_5^2
w_5 <- as.vector(mu_5^2)
W_5 <-  diag(w_5)

print(mu_5[92] - mu_5[38]) # Difference between two individuals with same bmi but different group

R_5_0 = 1 / (eta_5[92] - 1.96*sqrt(t(X_aug[92, ]) %*% solve(t(X_aug) %*% W_5 %*% X_aug) %*% X_aug[92, ]))
L_5_0 = 1 / (eta_5[92] + 1.96*sqrt(t(X_aug[92, ]) %*% solve(t(X_aug) %*% W_5 %*% X_aug) %*% X_aug[92, ]))

R_5_1 = 1 / (eta_5[38] - 1.96*sqrt(t(X_aug[38, ]) %*% solve(t(X_aug) %*% W_5 %*% X_aug) %*% X_aug[38, ]))
L_5_1 = 1 / (eta_5[38] + 1.96*sqrt(t(X_aug[38, ]) %*% solve(t(X_aug) %*% W_5 %*% X_aug) %*% X_aug[38, ]))

print(paste('(',L_5_0 - L_5_1,',',R_5_0 - R_5_1,')')) # Confidence interval


# Deviance
deviance <- function(mu_hat, y){
  dev <- sum(1/mu_hat * y - log(mu_hat))
}

dev_4 <- 2 * (deviance(y, y) - deviance(mu_4, y))
dev_5 <- 2 * (deviance(y, y) - deviance(mu_5, y))

dev <- 2 * (deviance(mu_5, y) - deviance(mu_4, y))


models_comparison <- data.frame(
  Model = c("fit1", "fit2"),
  Deviance = c(dev_4, dev_5),
  AIC = c(AIC(fit1), AIC(fit2))
)

models_comparison


  # 7. For the clinician's model
plot(dat$bmi, mu_4, xlab = "BMI", ylab = "Expected increase in heart rate", 
     col = treatment_colors[as.numeric(as.character(dat$treatment)) + 1], pch=16)

# Add legend for treatment groups
legend("topright", legend = levels(dat$treatment), col = treatment_colors, 
       pch = 16, cex=0.7, title = "Treatment")

# Add guide lines
library(KernSmooth)
LR.fit0 <- locpoly(dat$bmi[dat$treatment == 0], mu_4[dat$treatment == 0], degree = 3,
                   bandwidth = 0.5,
                   kernel = "normal")
LR.fit1 <- locpoly(dat$bmi[dat$treatment == 1], mu_4[dat$treatment == 1], degree = 3,
                   bandwidth = 0.5,
                   kernel = "normal")

lines(LR.fit0$x, LR.fit0$y, type = "s", col = "red", lwd = 1, lty = 3)
lines(LR.fit1$x, LR.fit1$y, type = "s", col = "blue", lwd = 1, lty = 3)


  # 8. Simulation study
n_sim <- 1000
treat_effs <- seq(0, 10, length.out = 20)
power <- matrix(0, nrow = length(treat_effs), ncol = 1)
j <- 1

for(treat_eff in treat_effs){
  p_store <- matrix(0, nrow = n_sim, ncol = 1) 
  for(i in 1:n_sim){
    dat$diff_sim <- sample(unlist(simulate(fit1)))
    wh_treat <- dat$treatment == 0
    dat$diff_sim[wh_treat] <- dat$diff_sim[wh_treat] + treat_eff
    fit_sim1 <- glm(diff_sim ~ bmi + treatment, data=dat, family = Gamma(link = "inverse"))
    sum1 <- summary(fit_sim1)
    p_store[i,1] <- sum1$coefficients[3,4]
  }
 
  power[j] <- colMeans(p_store < 0.05)
  j <- j + 1
}


power_long <- data.frame(
  effect = treat_effs,
  power = power
)

power_plot <- ggplot(data = power_long, mapping = aes(x = effect,
                                                      y = power) ) + 
                    geom_line(col = 'red') +
                    geom_hline(yintercept = 0.05, lty = 2) 

power_plot




# Create summary table
summary_table <- dat %>%
  group_by(as.numeric(as.character(treatment))) %>%
  summarise(
    n = n(),
    Mean_rest_pulse = mean(rest_pulse),
    Sd_rest_pulse = sd(rest_pulse),
    Mean_stimulated_pulse = mean(stimulated_pulse),
    Sd_stimulated_pulse = sd(stimulated_pulse),
    Mean_bmi = mean(bmi),
    Sd_bmi = sd(bmi)
  )

# Convert the summary table to a dataframe
summary_table_df <- as.data.frame(summary_table)
summary_table_df <- round(summary_table_df, 2)

# Print the formatted table with lines between columns and rows
summary_table_html <- kable(summary_table_df, format = "html", 
                            caption = "Summary Statistics by Treatment Group",
                            align = "c",
                            col.names = c("Treatment Group", "Count", "Mean Rest Pulse", "SD Rest Pulse", "Mean Stimulated Pulse", "SD Stimulated Pulse", "Mean BMI", "SD BMI"),
                            escape = FALSE) %>%
  as.character()  # Convert to character to manipulate HTML

# Add CSS classes and additional attributes to the table
summary_table_html <- gsub("<table", "<table class='table table-bordered table-striped' style='border-collapse: separate; border-spacing: 0px;'", summary_table_html)
summary_table_html <- gsub("<th", "<th style='border: 1px solid black; padding: 5px;'", summary_table_html)
summary_table_html <- gsub("<td", "<td style='border: 1px solid black; padding: 5px;'", summary_table_html)

# Print the HTML table
cat(summary_table_html)
