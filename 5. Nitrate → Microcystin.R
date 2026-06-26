##############################################################
#Run all thre models (Temp,Nitrate and Biomass - MC) 
#using random-effects, multilevel meta-analyses using metafor
##############################################################



##############################################################
# Random-effects multilevel meta-analysis
# Nitrate → Microcystin
##############################################################

library(dplyr)
library(metafor)

##------------------------------------------------------------
## 1. Ensure StudyID is a factor
##------------------------------------------------------------
nit_data$StudyID <- factor(nit_data$StudyID)

##------------------------------------------------------------
## 2. BASELINE MODEL (no moderator)
##------------------------------------------------------------

res_nit_base <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = nit_data,
  method = "REML"
)

summary(res_nit_base)

##############################################################
## Effect size interpretation (baseline)
##############################################################

# Back-transform pooled effect
nit_r  <- tanh(res_nit_base$b[1])
nit_ci <- tanh(c(res_nit_base$ci.lb, res_nit_base$ci.ub))

nit_r
nit_ci

##############################################################
## Heterogeneity diagnostics (baseline)
##############################################################

res_nit_base$sigma2
confint(res_nit_base)


##############################################################
### VERIFY HETEROGENEITY WITH MODERATOR(PEPTIDE TYPE)
### Nitrate → Microcystin
### Baseline = Intracellular microcystin
##############################################################

##------------------------------------------------------------
## 3. Sanity checks on moderator
##------------------------------------------------------------

stopifnot("Type.of.Peptide.Measured" %in% names(nit_data))

table(nit_data$Type.of.Peptide.Measured)

##------------------------------------------------------------
## 4. FORCE Intracellular as the reference level
##------------------------------------------------------------

nit_data$Type.of.Peptide.Measured <- factor(
  nit_data$Type.of.Peptide.Measured,
  levels = c("Intracellular", "Extracellular", "Total")
)

levels(nit_data$Type.of.Peptide.Measured)
table(nit_data$Type.of.Peptide.Measured)

##------------------------------------------------------------
## 5. MODERATOR MODEL
##------------------------------------------------------------

res_nit_mod <- rma.mv(
  yi, vi,
  mods   = ~ Type.of.Peptide.Measured,
  random = ~ 1 | StudyID,
  data   = nit_data,
  method = "REML"
)

summary(res_nit_mod)

##------------------------------------------------------------
## 6. Compare heterogeneity
##------------------------------------------------------------

res_nit_base$sigma2
res_nit_mod$sigma2

anova(res_nit_base, res_nit_mod, refit = TRUE)

##############################################################
## 7. Extract absolute r and CI for each peptide type
##############################################################

## Fixed-effect estimates (Fisher's Z)
coef_z <- coef(res_nit_mod)

## Variance–covariance matrix
V <- vcov(res_nit_mod)

coef_z
V

##------------------------------------------------------------
## 8. Fisher’s Z per peptide type
##------------------------------------------------------------

# Intracellular (baseline)
z_intra <- coef_z["intrcpt"]

# Extracellular
z_extra <- coef_z["intrcpt"] +
  coef_z["Type.of.Peptide.MeasuredExtracellular"]

# Total
z_total <- coef_z["intrcpt"] +
  coef_z["Type.of.Peptide.MeasuredTotal"]

##------------------------------------------------------------
## 9. Standard errors (accounting for covariance)
##------------------------------------------------------------

se_intra <- sqrt(V["intrcpt", "intrcpt"])

se_extra <- sqrt(
  V["intrcpt", "intrcpt"] +
    V["Type.of.Peptide.MeasuredExtracellular",
      "Type.of.Peptide.MeasuredExtracellular"] +
    2 * V["intrcpt", "Type.of.Peptide.MeasuredExtracellular"]
)

se_total <- sqrt(
  V["intrcpt", "intrcpt"] +
    V["Type.of.Peptide.MeasuredTotal",
      "Type.of.Peptide.MeasuredTotal"] +
    2 * V["intrcpt", "Type.of.Peptide.MeasuredTotal"]
)

##------------------------------------------------------------
## 10. 95% CI on Z scale
##------------------------------------------------------------

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

##------------------------------------------------------------
## 11. Back-transform to Pearson’s r
##------------------------------------------------------------

r_intra <- tanh(z_intra)
r_extra <- tanh(z_extra)
r_total <- tanh(z_total)

ci_intra_r <- tanh(ci_intra_z)
ci_extra_r <- tanh(ci_extra_z)
ci_total_r <- tanh(ci_total_z)

##------------------------------------------------------------
## 12. Final working results table (NOT for paper yet)
##------------------------------------------------------------

results_nitrate <- data.frame(
  Peptide_Type = c("Intracellular", "Extracellular", "Total"),
  r = c(r_intra, r_extra, r_total),
  CI_lower = c(ci_intra_r[1], ci_extra_r[1], ci_total_r[1]),
  CI_upper = c(ci_intra_r[2], ci_extra_r[2], ci_total_r[2])
)

results_nitrate

summary(res_nit_base)

summary(res_nit_mod)

anova(res_nit_base, res_nit_mod)




########################################
#NITRATE -BIOMASS
#######################################
##############################################################
# STEP B: Nitrate → Biomass
# Random-effects multilevel meta-analysis
##############################################################

library(dplyr)
library(metafor)

##------------------------------------------------------------
## 1. Prepare the Nitrate–Biomass dataset
##------------------------------------------------------------

# Subset biomass data to Nitrate predictor
biomass_nit_data <- biomass_data %>%
  filter(Predictor == "Nitrate")

# Ensure StudyID is a factor
biomass_nit_data$StudyID <- factor(biomass_nit_data$StudyID)

# Sanity checks
nrow(biomass_nit_data)
summary(biomass_nit_data$vi)
stopifnot(all(is.finite(biomass_nit_data$vi)))

##------------------------------------------------------------
## 2. Baseline multilevel model
##------------------------------------------------------------

res_nit_biomass <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = biomass_nit_data,
  method = "REML"
)

summary(res_nit_biomass)

##------------------------------------------------------------
## 3. Back-transform effect size (PRIMARY RESULT)
##------------------------------------------------------------

# Pooled effect (Pearson's r)
nit_biomass_r <- tanh(res_nit_biomass$b[1])

# 95% CI
nit_biomass_ci <- tanh(c(
  res_nit_biomass$ci.lb,
  res_nit_biomass$ci.ub
))

nit_biomass_r
nit_biomass_ci

##------------------------------------------------------------
## 4. Heterogeneity diagnostics (SECONDARY)
##------------------------------------------------------------

# Between-study variance (tau^2)
res_nit_biomass$sigma2

# CI for tau^2
confint(res_nit_biomass)


summary(res_nit_biomass)
nit_biomass_r
nit_biomass_ci
res_nit_biomass$sigma2