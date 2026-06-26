##############################################################
#Run all thre models (Temp,Nitrate and Biomass - MC) 
#using random-effects, multilevel meta-analyses using metafor
##############################################################

##############################################################
# Random-effects multilevel meta-analysis
# Temperature → Microcystin
##############################################################

# Load required packages 
library(dplyr)
library(metafor)

## Ensure StudyID is a factor (important for multilevel models)
temp_data$StudyID <- factor(temp_data$StudyID)

## Fit multilevel random-effects model
res_temp <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = temp_data,
  method = "REML"
)

## Model summary
summary(res_temp)

##################################################
## Effect size interpretation (PRIMARY RESULT)
## Back-transform Fisher’s Z → Pearson’s r
##################################################

# Pooled effect (Pearson's r)
temp_r <- tanh(res_temp$b[1])

# 95% confidence interval for effect size
temp_ci <- tanh(c(res_temp$ci.lb, res_temp$ci.ub))

temp_r
temp_ci

##################################################
## Heterogeneity diagnostics (SECONDARY)
##################################################

# Between-study variance (tau^2)
res_temp$sigma2

# Confidence interval for tau^2
confint(res_temp)






############################################################
### VERIFY HETEROGENEITY WITH MODERATOR
### Temperature → Microcystin
### Baseline = Intracellular microcystin
############################################################

library(dplyr)
library(metafor)

##----------------------------------------------------------
## 1. Sanity checks on moderator
##----------------------------------------------------------

# Check the moderator exists
stopifnot("Type.of.Peptide.Measured" %in% names(temp_data))

# Check categories and counts
table(temp_data$Type.of.Peptide.Measured)

##----------------------------------------------------------
## 2. FORCE Intracellular to be the reference level
##    (THIS IS THE KEY CORRECTION)
##----------------------------------------------------------

temp_data$Type.of.Peptide.Measured <- factor(
  temp_data$Type.of.Peptide.Measured,
  levels = c("Intracellular", "Extracellular", "Total")
)

# Confirm ordering
levels(temp_data$Type.of.Peptide.Measured)
table(temp_data$Type.of.Peptide.Measured)

##----------------------------------------------------------
## 3. Baseline model (no moderator)
##----------------------------------------------------------

res_temp_base <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = temp_data,
  method = "REML"
)

summary(res_temp_base)

##----------------------------------------------------------
## 4. Moderator model (Intracellular as baseline)
##----------------------------------------------------------

res_temp_mod <- rma.mv(
  yi, vi,
  mods   = ~ Type.of.Peptide.Measured,
  random = ~ 1 | StudyID,
  data   = temp_data,
  method = "REML"
)

summary(res_temp_mod)

##----------------------------------------------------------
## 5. Compare heterogeneity
##----------------------------------------------------------

# Tau^2 before and after
res_temp_base$sigma2
res_temp_mod$sigma2

# Likelihood ratio test
anova(res_temp_base, res_temp_mod, refit = TRUE)

##----------------------------------------------------------
## 6. Back-transform subgroup effects
##----------------------------------------------------------

# Extract coefficients (Fisher’s Z)
coef_z <- coef(res_temp_mod)

# Back-transform to Pearson’s r
coef_r <- tanh(coef_z)

coef_z
coef_r

##----------------------------------------------------------
## 7. Confidence intervals for subgroup effects
##----------------------------------------------------------

# Confidence intervals on Z scale
ci_z <- cbind(
  ci.lb = res_temp_mod$ci.lb,
  ci.ub = res_temp_mod$ci.ub
)

# Back-transform CIs to r
ci_r <- tanh(ci_z)
ci_z
ci_r

summary(res_temp_mod)
coef_r
ci_r



#####################################################
##Extract coefficients and variance–covariance matrix
#####################################################
# Fixed-effect estimates (Fisher's Z)
coef_z <- coef(res_temp_mod)

# Variance–covariance matrix of fixed effects
V <- vcov(res_temp_mod)

coef_z
V


####Compute Fisher’s Z for each peptide type
# Intracellular (baseline)
z_intra <- coef_z["intrcpt"]

# Extracellular = intercept + extracellular difference
z_extra <- coef_z["intrcpt"] +
  coef_z["Type.of.Peptide.MeasuredExtracellular"]

# Total = intercept + total difference
z_total <- coef_z["intrcpt"] +
  coef_z["Type.of.Peptide.MeasuredTotal"]


##Compute standard errors correctly
# Standard error for intracellular
se_intra <- sqrt(V["intrcpt", "intrcpt"])

# Standard error for extracellular
se_extra <- sqrt(
  V["intrcpt", "intrcpt"] +
    V["Type.of.Peptide.MeasuredExtracellular",
      "Type.of.Peptide.MeasuredExtracellular"] +
    2 * V["intrcpt", "Type.of.Peptide.MeasuredExtracellular"]
)

# Standard error for total
se_total <- sqrt(
  V["intrcpt", "intrcpt"] +
    V["Type.of.Peptide.MeasuredTotal",
      "Type.of.Peptide.MeasuredTotal"] +
    2 * V["intrcpt", "Type.of.Peptide.MeasuredTotal"]
)


###Build 95% confidence intervals on the Z scale
z_crit <- qnorm(0.975)

ci_intra_z <- c(
  z_intra - z_crit * se_intra,
  z_intra + z_crit * se_intra
)

ci_extra_z <- c(
  z_extra - z_crit * se_extra,
  z_extra + z_crit * se_extra
)

ci_total_z <- c(
  z_total - z_crit * se_total,
  z_total + z_crit * se_total
)


####Back-transform to Pearson’s r
# Point estimates
r_intra <- tanh(z_intra)
r_extra <- tanh(z_extra)
r_total <- tanh(z_total)

# Confidence intervals
ci_intra_r <- tanh(ci_intra_z)
ci_extra_r <- tanh(ci_extra_z)
ci_total_r <- tanh(ci_total_z)

# Combine into a clean table
results <- data.frame(
  Peptide_Type = c("Intracellular", "Extracellular", "Total"),
  r = c(r_intra, r_extra, r_total),
  CI_lower = c(ci_intra_r[1], ci_extra_r[1], ci_total_r[1]),
  CI_upper = c(ci_intra_r[2], ci_extra_r[2], ci_total_r[2])
)

results





############################################
# TEMPERATURE AND BIOMASS
############################################
##############################################
# STEP A: Temperature → Biomass
# Random-effects multilevel meta-analysis
##############################################
library(dplyr)
library(metafor)

##------------------------------------------------------------
## 1. Prepare the Temperature–Biomass dataset
##------------------------------------------------------------

# Ensure StudyID is a factor
biomass_temp_data <- biomass_data %>%
  filter(Predictor == "Temperature")

biomass_temp_data$StudyID <- factor(biomass_temp_data$StudyID)

# Sanity checks
nrow(biomass_temp_data)
summary(biomass_temp_data$vi)
stopifnot(all(is.finite(biomass_temp_data$vi)))

##------------------------------------------------------------
## 2. Baseline multilevel model
##------------------------------------------------------------

res_temp_biomass <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = biomass_temp_data,
  method = "REML"
)

summary(res_temp_biomass)

##------------------------------------------------------------
## 3. Back-transform effect size (PRIMARY RESULT)
##------------------------------------------------------------

# Pooled effect (Pearson's r)
temp_biomass_r <- tanh(res_temp_biomass$b[1])

# 95% CI
temp_biomass_ci <- tanh(c(
  res_temp_biomass$ci.lb,
  res_temp_biomass$ci.ub
))

temp_biomass_r
temp_biomass_ci

##------------------------------------------------------------
## 4. Heterogeneity diagnostics (SECONDARY)
##------------------------------------------------------------

# Between-study variance (tau^2)
res_temp_biomass$sigma2

# CI for tau^2
confint(res_temp_biomass)


summary(res_temp_biomass)
temp_biomass_r
temp_biomass_ci
res_temp_biomass$sigma2