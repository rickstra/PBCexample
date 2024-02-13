library(tidyverse)
library(rms)
library(table1)

set.seed(8675309)

d <- subset(survival::pbc, !is.na(trt))  # Exclude subjects not randomized

{
  d$trt     <- factor(d$trt, levels=1:2, labels=c("D-penicillamine", "Placebo"))
  d$sex     <- factor(d$sex, levels=c("m", "f"), labels=c("Male", "Female"))
  d$stage   <- factor(d$stage, levels=1:4, labels=paste("Stage", 1:4))
  d$edema   <- factor(d$edema, levels=c(0, 0.5, 1),
                      labels=c("No edema",
                               "Untreated or successfully treated",
                               "Edema despite diuretic therapy"))
  d$spiders <- factor(d$spiders, levels=c(0, 1), labels=c("No", "Yes"))
  d$hepato <- factor(d$hepato, levels=c(0, 1), labels=c("No", "Yes"))
  d$ascites <- factor(d$ascites, levels=c(0, 1), labels=c("No", "Yes"))
  d$stage4 <- factor(d$stage == "Stage 4", levels = c(FALSE, TRUE), 
                     labels = c("Fibrosis stage 1-3", "Fibrosis stage 4"))
}
{
  label(d$age)      <- "Age (y)"
  label(d$sex)      <- "Sex"
  label(d$stage)    <- "Histologic stage of disease"
  label(d$edema)    <- "Edema status"
  label(d$spiders)  <- "Presence of spiders"
  label(d$hepato)   <- "Presence of hepatomegaly"
  label(d$ascites)  <- "Presence of ascites"
  label(d$platelet) <- "Platelet count (10^9 per liter)"
  label(d$protime)  <- "Standardised blood clotting time"
  label(d$albumin)  <- "Serum albumin (g/dL)"
  label(d$alk.phos) <- "Alkaline phosphotase (U/L)"
  label(d$ast)      <- "Aspartate aminotransferase (U/mL)"
  label(d$bili)     <- "Serum bilirubin (mg/dL)"
  label(d$chol)     <- "Serum cholesterol (mg/dL)"
  label(d$copper)   <- "Urine copper (&mu;g/day)"
  label(d$trig)     <- "Triglycerides (mg/dL)"
}

#-------------------------------------------------------------------------------

d <- d %>% filter(ascites == "No") %>% # Exclude ascites patients
  rename(platelets = platelet,         # Rename some predictors
         cholesterol = chol, 
         bilirubin = bili, 
         triglycerides = trig)


#-------------------------------------------------------------------------------
# Table 1
tab1 <- table1(~ age + sex + 
                 edema + spiders + hepato + 
                 platelets + protime + albumin + alk.phos + ast + bilirubin + cholesterol +
                 copper + triglycerides | as.factor(stage4), data = d, 
               render = function(x, ...) render.default(x, render.continuous = function (x, ...){
                 with(stats.apply.rounding(stats.default(x, ...), ...), 
                      c("", `Mean (SD)` = sprintf("%s (%s)", MEAN, SD), 
                        `Median [IQR]` = sprintf("%s [%s, %s]", MEDIAN, Q1, Q3)
                      ))
               }))

tab1
#-------------------------------------------------------------------------------
# Single impute the missing values
# It is however recommended to use multiple imputation

d <- complete(mice::mice(d, m = 1, printFlag = FALSE), 1)

#-------------------------------------------------------------------------------
# Create new outcome variables
d$stage4 <- as.numeric(d$stage == "Stage 4")

# Rename the levels of some predictors for later
levels(d$hepato) <- c("no hepatomegaly", "hepatomegaly")
levels(d$spiders) <- c("no spiders", "spiders")

n <- nrow(d)

#-------------------------------------------------------------------------------
# Define the formula for the maximal model, 
# which includes 3-knot splines for the continuous predictors
(f_full <- formula(stage4 ~ sex + rcs(age, 3) + rcs(bilirubin, 3) + 
                     rcs(albumin, 3) + rcs(ast, 3) + rcs(platelets, 3) + 
                     rcs(cholesterol, 3) + rcs(alk.phos, 3) + 
                     rcs(triglycerides, 3) + spiders + hepato + edema))
# Fit the maximal model (logistic regression)
m_full <- lrm(f_full, d)

# Inspect the fitted model
print(m_full, coefs = FALSE)

#-------------------------------------------------------------------------------

step1 <- function(y # The response variable
                  ) {
  p <- mean(y)
  n_required <- p * (1 - p) * (1.96 / 0.05)^2
  delta_achieved <- 1.96 * sqrt(p * (1 - p) / length(y))
  c("Sample size" = round(n_required), Value = round(delta_achieved, 3), "df" = NA)
}

step2 <- function(y, # The response variable
                  k  # The desired degrees of freedom
                  ) {
  n <- length(y)
  p <- mean(y)
  mape <- exp(-0.508 - 0.544 * log(n) + 
                0.259 * log(p) + 0.504 * log(k))
  n_required <- exp((-0.508 - log(0.05) + 
                       0.259 * log(p) + 0.504 * log(k)) / 0.544)
  k_allowed <- exp((-0.508 - log(0.05) + 
                      0.259 * log(p) - 0.544 * log(n)) / -0.504)
  c("Sample size" = round(n_required), Value = round(mape, 3), "df" = round(k_allowed, 1))
}

step3 <- function(n, # acquired sample size
                  m_chisq, # the model chi-squared for the full model
                  k, # desired degrees of freedom
                  opt = 0.8 # Which proportion of the full model's R2 do you expect the k df model to have?
                  ) {
  S <- 0.9
  R2cs <- (1 - exp(-(as.numeric(m_chisq) / n))) * opt
  n_required <- k / ((S - 1) * log(1 - R2cs / S))
  k_allowed <- n * ((S - 1) * log(1 - R2cs / S))
  S <- uniroot(function(s) n * ((s - 1) * log(1 - R2cs / s)) - k, c(0.5, 1))$root
  c("Sample size" = round(n_required), Value = round(S, 3), "df" = round(k_allowed, 1))
}

step4 <- function(y, 
                  m_chisq, # the model chi-squared for the full model
                  k, 
                  opt = 0.8 # Which proportion of the full model's R2 do you expect the k df model to have? 
                  ) {
  p <- mean(y)
  n <- length(y)
  logLnull <- p * log(p) + (1 - p) * log(1 - p)
  maxR2 <- 1 - exp(2 * logLnull)
  R2cs <- (1 - exp(-(as.numeric(m_chisq) / n))) * opt
  
  S <- R2cs / (R2cs + 0.05 * maxR2)
  
  n_required <- k / ((S - 1) * log(1 - R2cs / S))
  k_allowed <- n * ((S - 1) * log(1 - R2cs / S))
  
  opt <- uniroot(function(d) {
    S <- R2cs / (R2cs + d * maxR2)
    k / ((S - 1) * log(1 - R2cs / S)) - n
  }, c(0.01, 0.7))
  
  c("Sample size" = round(n_required), Value = round(opt$root, 3), "d.f." = round(k_allowed, 1))  
}
# Calculate based on 14 degrees of freedom
k_cand <- 14
R2frac <- 0.8 #k_cand/m_full$stats[4]
sample_size <- data.frame(Criterion = 1:4, 
                          Description = c("Precisely estimate the overall risk", 
                                          "Small mean absolute prediction error", 
                                          'Small estimated "shrinkage factor"', 
                                          "Small expected optimism of the model fit"),
                          Target = c("<0.05", "<0.05", ">0.90", "<0.05"), 
                          rbind(step1(d$stage4),
                                step2(d$stage4, k_cand),
                                step3(length(d$stage4), m_full$stats[3], k_cand, R2frac),
                                step4(d$stage4, m_full$stats[3], k_cand, R2frac))
)
sample_size
#-------------------------------------------------------------------------------
# Rank the predictors based on their predictive power to determine where to 
# spend more than 1 degree of freedom
plot(anova(m_full), margin = c("d.f.", "chisq"))

# Combine categories for edema -> -1 df
d$edema01 <- factor(as.numeric(d$edema != "No edema"), levels = 0:1, 
                    labels = c("No", "Yes"))
label(d$edema01) <- label(d$edema)

# Define a new model with all predictors where we model 
# bilirubin, albumin, and cholesterol using 3-knot splines 
# using 2 degrees of freedom each, which allows for non-linearity
(f <- formula(stage4 ~ ast  + triglycerides + sex + spiders + edema01 + 
                platelets + age + alk.phos + rcs(bilirubin, 3) + 
                rcs(albumin, 3) + rcs(cholesterol, 3) + hepato))
m <- lrm(f, d, x = T, y = T)

# Inspect fitted model
print(m, coefs = FALSE)

# We want the AIC to be lower for the reduced model than for the maximal model
# This indicates that the degrees of freedom we cut didn't lose us too much 
# predictive power
c(AIC_full = AIC(m_full), AIC_reduced = AIC(m))


# The model still has 15 degrees of freedom, which is too many for our sample size
# Use penalized maximum likelihood estimation to shrink the estimates and 
# compensate for the overfitting

# First find the optimal shrinkage parameter with a grid search
pen <- pentrace(m, list(simple = seq(0, 10, by = 1), 
                        nonlinear = seq(0, 50, by = 1)))
opt.pen <- pen$penalty
# One more iteration
pen <- pentrace(m, 
         list(simple = do.call(seq, as.list(c(opt.pen$simple + c(-0.5, 0.5), 0.1))), 
              nonlinear = do.call(seq, as.list(c(opt.pen$nonlinear + c(-0.5, 0.5), 0.1)))))
(opt.pen <- pen$penalty)

# Fit the model with penalized regression
m_pen <- update(m, penalty = list(simple = opt.pen$simple, 
                                  nonlinear = opt.pen$nonlinear))
print(m_pen, coefs = F)

# Predict from both the penalized and unpenalized model
d$lp0 <- predict(m, d, type = "lp")
d$p <- predict(m_pen, d, type = "fitted")
d$lp <- predict(m_pen, d, type = "lp")

# Find a simplified approximate model through ordinary linear regression on the 
# predicted values
(ols_f <- formula(paste("lp~", as.character(f))[3]))
apx <- ols(ols_f, data = d, 
           sigma = 1, x = TRUE)
fastbw(apx, aics = 1e10)
# age, spiders, platelets, cholesterol, albumin, and hepatomegaly are enough 
# to explain the predictions from the full model with R2>=0.95

# Fit the final approximated model
m_apx <- ols(lp ~ hepato + rcs(albumin, 3) + platelets + 
      rcs(cholesterol, 3) + spiders + age, 
    data = d, x = TRUE)

# Generate predicted values
d$lp2 <- predict(m_apx, d, type = "lp")
d$p2 <- plogis(d$lp2)

#-------------------------------------------------------------------------------
# Optionally force the approximated model into a logistic regression object
m_apx <- lrm(stage4 ~ hepato + rcs(albumin, 3) + platelets + 
                      rcs(cholesterol, 3) + spiders + age, 
             data = d, x = TRUE, 
             penalty = list(simple = opt.pen$simple, 
                            nonlinear = opt.pen$nonlinear))

# Replace the estimates with the appropriate ones based on the penalized model 
X <- cbind(Intercept = 1, m_pen$x) # full model design
Z <- cbind(Intercept = 1, m_apx$x) # approx. model design 
W <- solve(t(Z) %*% Z, t(Z)) %*% X # contrast matrix 
V <- vcov(m_pen) # var(full model)
V_apx <- W %*% V %*% t(W)

m_apx$coefficients <- W %*% m_pen$coefficients
m_apx$var <- V_apx

# Genereate predictions from the approximated model as a logistic model
d$lp2 <- predict(m_apx, d, type = "lp")
d$p2 <- predict(m_apx, d, type = "fitted")

latex(m_apx, digits = 3)

#-------------------------------------------------------------------------------
# Define the performance measure functions
#-------------------------------------------------------------------------------
# Discrimination
EstimateDiscrimination <- function(p, y, plot = FALSE, ...) {
  y <- y[order(-p)]
  sens <- cumsum(y) / sum(y)
  spec <- 1-cumsum(!y) / sum(!y)
  auc <- sum((matrix(-sort(-p)[!y], sum(y), sum(!y), byrow = TRUE) - 
                (-sort(-p)[!!y])) < 0) / sum(y) / sum(!y)
  if(plot) {
    plot(1-spec, sens, type = 'l', ...)
    abline(a = 0, b=1, lty = 3)
    text(0.7, 0.3, paste("AUC =", round(auc, 3)))
  }
  data.frame(sens, spec, auc)
}

#-------------------------------------------------------------------------------
# Utility (Net Benefit)
EstimateNetBenefit <- function(p, y, pt, standardized = FALSE, plot = FALSE) {
  n <- length(y)
  s <- ifelse(standardized, mean(y), 1)
  
  p_pos <- p[y == 1]
  p_neg <- p[y == 0]
  
  TP <- rowSums((matrix(p_pos, length(pt), length(p_pos), byrow = T) - pt) >= 0)
  FP <- rowSums((matrix(p_neg, length(pt), length(p_neg), byrow = T) - pt) >= 0)
  
  NB <- (TP - pt / (1 - pt) * FP) / n / s
  NB0 <- (length(p_pos) - pt / (1 - pt) * length(p_neg)) / n / s
  NB_def <- pmax(NB0, 0)
  TT <- 1 / (NB - NB_def)
  
  out <- data.frame(pt = pt, NB = NB, NB0 = NB0, NB_def = NB_def, TT = TT)
  attr(out, "standardized") <- standardized
  
  if(plot) with(out, {
    plot(pt, NB, col = "gray", type = 'l', xlim = c(0, 1), ylim = c(-0.01, max(NB0)), 
         xlab = "cut-off", ylab = "Net benefit")
    lines(lowess(pt, NB, 0.2))
    lines(pt, NB0, lty = 2)
    abline(h = 0, lty = 3)
  })
  out
} 
#-------------------------------------------------------------------------------
# Perform Bootstrap validation
B <- 500
set.seed(8675309)
bs_out <- list()
d <- d[order(d$p2), ]
for (i in 1:B) {
  # Create new BS sample and fit model
  bd <- d[sample(1:nrow(d), nrow(d), replace = TRUE), ]
  bm <- update(m_pen, data = bd)
  
  term <- c(ast = "ast", triglycerides = "triglycerides", 
    sex = "sex", spiders = "spiders", edema01 = "edema01", 
    platelets = "platelets", age = "age", alk.phos = "alk.phos", 
    bilirubin = "rcs(bilirubin, 3)" , albumin = "rcs(albumin, 3)", 
    cholesterol = "rcs(cholesterol, 3)", hepato = "hepato")
  
  (ols_f <- formula(paste("blp~", as.character(f))[3]))
  apx <- ols(ols_f, data = bd, 
             sigma = 1, x = TRUE)
  
  # Generate prediction in both the BS sample and the original sample 
  bd$blp <- predict(bm, bd, type = "lp")
  d$vlp <- predict(bm, d, type = "lp")
  bp <- plogis(bd$blp)
  vp <- plogis(bd$vlp)
  
  # Calibration error
  cal_b <- lowess(bp, bd$stage4, iter = 0)
  cal_v <- lowess(vp, d$stage4, iter = 0)
  obs_b <- approx(cal_b, xout=d$p2, ties=function(x)x[1])$y
  obs_v <- approx(cal_v, xout=d$p2, ties=function(x)x[1])$y
  err_b <- obs_b - d$p2
  err_v <- obs_v - d$p2
  cal_opt <- err_b - err_v 
  
  # O:E Ratio
  oe_b <- sum(bd$stage4) / sum(bp)
  oe_v <- sum(d$stage4) / sum(vp)
  
  # Calibration Intercept + Slope
  cal_weak_b <- lrm(stage4 ~ blp, data = bd)$coefficients
  cal_weak_v <- lrm(stage4 ~ vlp, data = d)$coefficients
  
  # AUC
  auc_b <- EstimateDiscrimination(bp, bd$stage4)$auc[1]
  auc_v <- EstimateDiscrimination(vp, d$stage4)$auc[1]
  auc_opt <- auc_b - auc_v
  
  # Decision curves
  nb_b <- EstimateNetBenefit(bp, bd$stage4, seq(0, 1, by = 0.01))
  nb_v <- EstimateNetBenefit(vp, d$stage4, seq(0, 1, by = 0.01))
  nb_opt <- nb_b$NB - nb_v$NB
  
  
  bs_out[[i]] <- list(bp = bp, vp = vp, cal_weak_b = cal_weak_b, 
                      cal_weak_v = cal_weak_v, oe_b = oe_b, oe_v = oe_v, 
                      err_b = err_b, cal_opt = cal_opt, 
                      auc_b = auc_b, auc_opt = auc_opt, 
                      nb_b = nb_b, nb_opt = nb_opt)
}
auc_app <- EstimateDiscrimination(d$p2, d$stage4)$auc[1]
auc_opt <- mean(sapply(bs_out, function(lst) {lst$auc_opt}))

auc_95ci <- quantile(sapply(bs_out, function(lst) {lst$auc_b}), 
                     c(0.025, 0.975))
auc_95ci_corr <- auc_95ci - auc_opt

err_b <- do.call(rbind, lapply(bs_out, function(lst) {lst$err_b}))
l <- lowess(d$p2, d$stage4, iter = 0)
cal_opt <- do.call(rbind, lapply(bs_out, function(lst) {lst$cal_opt})) %>% colMeans(na.rm = T)
cal_df <- data.frame(pred = l$x, obs = c(l$y, l$y - cal_opt), 
                     t(apply(err_b, 2, function(col) quantile(col, c(0.025, 0.975),
                                                              na.rm = TRUE))) + d$p2,
                     model = rep(c("Apparent", "Corrected"), each = length(l$x)))
cal_df[-(1:n), 3:4] <- NA
colnames(cal_df)[3:4] <- c("lower", "upper")

tab_ci <- apply(abs(err_b), 1, function(row) c(ICI = mean(row, na.rm = T), 
                                               E50 = as.numeric(quantile(row, c(0.50), na.rm =T)), 
                                               E90 = as.numeric(quantile(row, c(0.90), na.rm = T)))) %>% 
  apply(1, function(row) quantile(row, c(0.025, 0.975))) %>% t %>% rbind(AUC = auc_95ci)

weak_bs <- sapply(bs_out, function(lst) lst$cal_weak_b) 

weak_opt <- (weak_bs - sapply(bs_out, function(lst) lst$cal_weak_v)) %>% rowMeans

# weak_bs %>% apply(1, function(row) quantile(row, c(0.025, 0.975))) %>% t

weak_app <- lrm(stage4 ~ lp2, d)$coefficients

cal_error <- abs(l$y - l$x)
cal_error_corr <- abs(l$y - cal_opt - l$x)

oe_opt <- mean(sapply(bs_out, function(lst) lst$oe_b) - sapply(bs_out, function(lst) lst$oe_v)) 


meas_tab <- cbind(Apparent = 
                    c(AUC = EstimateDiscrimination(d$p2, d$stage4)$auc[1],
                      ICI = mean(cal_error), 
                      E50 = as.numeric(quantile(cal_error, 0.5)), 
                      E90 = as.numeric(quantile(cal_error, 0.9)), 
                      weak_app),
                  rbind(tab_ci[c(4, 1:3), ], (weak_bs %>% apply(1, function(row) quantile(row, c(0.025, 0.975))) %>% t)),
                  Corrected = 
                    c(AUC = EstimateDiscrimination(d$p2, d$stage4)$auc[1] - auc_opt, 
                      ICI = mean(cal_error_corr), 
                      E50 = as.numeric(quantile(cal_error_corr, 0.5)), 
                      E90 = as.numeric(quantile(cal_error_corr, 0.9)), 
                      weak_app -weak_opt)
) %>% round(digits = 3)

rownames(meas_tab)[6] <- "Slope"

nb <- EstimateNetBenefit(d$p2, d$stage4, seq(0, 1, by = 0.01))
nb_opt <- rowMeans(sapply(bs_out, function(lst) {lst$nb_opt}), na.rm = TRUE)
nb_b <- do.call(rbind, lapply(bs_out, function(lst) {lst$nb_b$NB}))
nb_b0 <- do.call(rbind, lapply(bs_out, function(lst) {lst$nb_b$NB0}))
nb$NBc <- nb$NB - nb_opt
nb$NBs <- lowess(nb$pt, nb$NB, 0.2)$y
nb$NBcs <- lowess(nb$pt, nb$NBc, 0.2)$y

nb_df <- data.frame(pt = nb$pt, NB = c(nb$NBs, nb$NBcs), 
                    t(apply(nb_b, 2, function(col) quantile(col, c(0.025, 0.975), na.rm = TRUE))),
                    model = rep(c("Apparent", "Corrected"), each = nrow(nb)))
nb_df[-(1:nrow(nb)), 3:4] <- NA
colnames(nb_df)[3:4] <- c("lower", "upper")
nb_df$lower[1:nrow(nb)] <- lowess(nb$pt, nb_df$lower[1:nrow(nb)], 0.2)$y
nb_df$upper[1:nrow(nb)] <- lowess(nb$pt, nb_df$upper[1:nrow(nb)], 0.2)$y

nb0_df <- data.frame(pt = nb$pt, NB = nb$NB0, 
                     t(apply(nb_b0, 2, function(col) quantile(col, c(0.025, 0.975), na.rm = TRUE))))
nb0_df[-(1:nrow(nb)), 3:4] <- NA
colnames(nb0_df)[3:4] <- c("lower", "upper")
nb0_df$lower <- pmax(nb0_df$lower, -0.05)

tt_df <- data.frame(pt = nb$pt, tt = nb$TT,  
                    t(sapply(bs_out, function(lst) {lst$nb_b$TT})) %>% 
                      apply(MARGIN = 2, FUN = function(col) quantile(col, c(0.025, 0.975), na.rm = TRUE)) %>% t
)
colnames(tt_df)[3:4] <- c("lower", "upper")


# Table 1
tab1

# Table 2
sample_size

# Table 3
meas_tab

# Figures

ggplot(fit_df, aes(X, Y)) + geom_point(alpha = 0.5) + 
  geom_line(aes(y = Y_hat, color = fit)) + 
  scale_color_discrete(guide = FALSE) + 
  xlab("Predictor") + ylab("Outcome") + 
  facet_grid(cols = vars(fit_df$fit))

ggplot(d, aes(x = p2, fill = as.factor(stage4), color = as.factor(stage4))) + 
  geom_density(alpha = 0.4) + xlab("Predicted risk") + ylab("Density") + 
  scale_fill_discrete(name = "Cirrhosis", labels = c("No", "Yes")) + 
  scale_color_discrete(guide = FALSE) + 
  theme(legend.justification = c(1, 1), legend.position = c(0.9, 0.9))

ggplot(cal_df, aes(pred, obs, color = model, fill = model)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, color = NULL), alpha = 0.2) + 
  geom_line() + 
  geom_abline(intercept = 0, slope = 1, linetype = 3)  + 
  xlab("Predicted risk") + ylab("Observed risk") + 
  scale_fill_discrete(name = "Model Performance", 
                      labels = c("Apparent", "Corrected")) + 
  scale_color_discrete(guide = FALSE) + 
  theme(legend.justification = c(1, 1), legend.position = c(0.9, 0.3))

ggplot(nb, aes(x = pt)) + 
  geom_line(aes(y = NB), color = "pink", alpha = 0.5)  + 
  geom_line(aes(y = NBc), color = "lightblue", alpha = 0.5) + 
  geom_line(data = nb_df, mapping = aes(y = NB, color = model)) + 
  geom_ribbon(data = nb_df, 
              mapping = aes(ymin = lower, ymax = upper, 
                            fill = model, color = NULL), alpha = 0.1) + 
  geom_hline(yintercept = 0, linetype = 3) + 
  geom_line(aes(y = NB0), linetype = 2)  + 
  geom_ribbon(data = nb0_df, 
              mapping = aes(ymin = lower, ymax = upper,  color = NULL), 
              alpha = 0.05) + 
  scale_color_discrete(name = "Model Performance", 
                       labels = c("Apparent", "Corrected")) + 
  scale_fill_discrete(guide = FALSE) + 
  theme(legend.justification = c(1, 1), legend.position = c(0.9, 0.9)) + 
  ylim(-0.05, max(nb_df$upper)) + 
  xlab("Threshold probability") + ylab("Net Benefit")

plot(nomogram(m_apx, lp = FALSE, fun = function(x) 1/(1 + exp(-x)), 
              funlabel = "Predicted Risk", vnames = "labels"))



                                                                                                                                                                                                     
