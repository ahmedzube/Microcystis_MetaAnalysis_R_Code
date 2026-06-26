########################
#BIOMASS - MICROCYSTIS
########################
##############################################################
# STEP C1: Biomass → Microcystin
# Baseline random-effects multilevel meta-analysis
##############################################################

library(dplyr)
library(metafor)

##------------------------------------------------------------
## 1. Prepare the Biomass–Microcystin dataset
##------------------------------------------------------------



biomass_mc_data <- micro_data %>%
  filter(!is.na(yi), !is.na(vi))

# Ensure StudyID is a factor
biomass_mc_data$StudyID <- factor(biomass_mc_data$StudyID)

# Sanity checks
nrow(biomass_mc_data)
summary(biomass_mc_data$vi)
stopifnot(all(is.finite(biomass_mc_data$vi)))

##------------------------------------------------------------
## 2. Baseline multilevel model
##------------------------------------------------------------

res_biomass_mc <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = biomass_mc_data,
  method = "REML"
)

summary(res_biomass_mc)

##------------------------------------------------------------
## 3. Back-transform effect size (PRIMARY RESULT)
##------------------------------------------------------------

# Pooled effect (Pearson's r)
biomass_mc_r <- tanh(res_biomass_mc$b[1])

# 95% confidence interval
biomass_mc_ci <- tanh(c(
  res_biomass_mc$ci.lb,
  res_biomass_mc$ci.ub
))

biomass_mc_r
biomass_mc_ci

##------------------------------------------------------------
## 4. Heterogeneity diagnostics (SECONDARY)
##------------------------------------------------------------

# Between-study variance (tau^2)
res_biomass_mc$sigma2

# CI for tau^2
confint(res_biomass_mc)

summary(res_biomass_mc)
biomass_mc_r
biomass_mc_ci
res_biomass_mc$sigma2



##############################################################
# STEP C2: Biomass → Microcystin × Peptide compartment
# Baseline = Intracellular microcystin
##############################################################

library(dplyr)
library(metafor)

##------------------------------------------------------------
## 1. Sanity checks on moderator
##------------------------------------------------------------

stopifnot("Type.of.Peptide.Measured" %in% names(biomass_mc_data))

table(biomass_mc_data$Type.of.Peptide.Measured)

##------------------------------------------------------------
## 2. Force Intracellular as reference level
##------------------------------------------------------------

biomass_mc_data$Type.of.Peptide.Measured <- factor(
  biomass_mc_data$Type.of.Peptide.Measured,
  levels = c("Intracellular", "Extracellular", "Total")
)

levels(biomass_mc_data$Type.of.Peptide.Measured)

##------------------------------------------------------------
## 3. Baseline model (no moderator)
##------------------------------------------------------------

res_biomass_mc_base <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = biomass_mc_data,
  method = "REML"
)

summary(res_biomass_mc_base)

##------------------------------------------------------------
## 4. Moderator model
##------------------------------------------------------------

res_biomass_mc_mod <- rma.mv(
  yi, vi,
  mods   = ~ Type.of.Peptide.Measured,
  random = ~ 1 | StudyID,
  data   = biomass_mc_data,
  method = "REML"
)

summary(res_biomass_mc_mod)

##------------------------------------------------------------
## 5. Compare heterogeneity
##------------------------------------------------------------

res_biomass_mc_base$sigma2
res_biomass_mc_mod$sigma2

anova(res_biomass_mc_base, res_biomass_mc_mod, refit = TRUE)


##############################################################
# STEP C2: Biomass → Microcystin × Peptide compartment
# Baseline = Intracellular microcystin
##############################################################

library(dplyr)
library(metafor)

##------------------------------------------------------------
## 1. Sanity checks on moderator
##------------------------------------------------------------

stopifnot("Type.of.Peptide.Measured" %in% names(biomass_mc_data))

table(biomass_mc_data$Type.of.Peptide.Measured)

##------------------------------------------------------------
## 2. Force Intracellular as reference level
##------------------------------------------------------------

biomass_mc_data$Type.of.Peptide.Measured <- factor(
  biomass_mc_data$Type.of.Peptide.Measured,
  levels = c("Intracellular", "Extracellular", "Total")
)

levels(biomass_mc_data$Type.of.Peptide.Measured)

##------------------------------------------------------------
## 3. Baseline model (no moderator)
##------------------------------------------------------------

res_biomass_mc_base <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = biomass_mc_data,
  method = "REML"
)

summary(res_biomass_mc_base)

##------------------------------------------------------------
## 4. Moderator model
##------------------------------------------------------------

res_biomass_mc_mod <- rma.mv(
  yi, vi,
  mods   = ~ Type.of.Peptide.Measured,
  random = ~ 1 | StudyID,
  data   = biomass_mc_data,
  method = "REML"
)

summary(res_biomass_mc_mod)

##------------------------------------------------------------
## 5. Compare heterogeneity
##------------------------------------------------------------

res_biomass_mc_base$sigma2
res_biomass_mc_mod$sigma2

anova(res_biomass_mc_base, res_biomass_mc_mod, refit = TRUE)

summary(res_biomass_mc_mod)
anova(res_biomass_mc_base, res_biomass_mc_mod, refit = TRUE)
res_biomass_mc_base$sigma2
res_biomass_mc_mod$sigma2




##############################################################
# STEP C2b: Extract compartment-specific effects (Biomass → MC)
##############################################################

# Fixed-effect estimates (Fisher's Z)
coef_z <- coef(res_biomass_mc_mod)

# Variance–covariance matrix
V <- vcov(res_biomass_mc_mod)

##------------------------------------------------------------
## Compute Z estimates per compartment
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
## Correct standard errors
##------------------------------------------------------------

se_intra <- sqrt(V["intrcpt", "intrcpt"])

se_extra <- sqrt(
  V["intrcpt", "intrcpt"] +
    V["Type.of.Peptide.MeasuredExtracellular",
      "Type.of.Peptide.MeasuredExtracellular"] +
    2 * V["intrcpt",
          "Type.of.Peptide.MeasuredExtracellular"]
)

se_total <- sqrt(
  V["intrcpt", "intrcpt"] +
    V["Type.of.Peptide.MeasuredTotal",
      "Type.of.Peptide.MeasuredTotal"] +
    2 * V["intrcpt",
          "Type.of.Peptide.MeasuredTotal"]
)

##------------------------------------------------------------
## 95% confidence intervals (Z scale)
##------------------------------------------------------------

z_crit <- qnorm(0.975)

ci_intra_z <- c(z_intra - z_crit * se_intra,
                z_intra + z_crit * se_intra)

ci_extra_z <- c(z_extra - z_crit * se_extra,
                z_extra + z_crit * se_extra)

ci_total_z <- c(z_total - z_crit * se_total,
                z_total + z_crit * se_total)

##------------------------------------------------------------
## Back-transform to Pearson's r
##------------------------------------------------------------

results_biomass_mc <- data.frame(
  Peptide_Type = c("Intracellular", "Extracellular", "Total"),
  r = tanh(c(z_intra, z_extra, z_total)),
  CI_lower = tanh(c(ci_intra_z[1], ci_extra_z[1], ci_total_z[1])),
  CI_upper = tanh(c(ci_intra_z[2], ci_extra_z[2], ci_total_z[2]))
)

results_biomass_mc



########Biomass - MC
res_bio<- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = biomass_mc_data,
  method = "REML"
)

summary(res_bio)

##################################################
##Back-transform to Pearson’s r
##################################################
bio_r <- tanh(res_bio$b[1])
bio_ci <- c(tanh(res_bio$ci.lb), tanh(res_bio$ci.ub))

bio_r
bio_ci

##Calculate Heterogeneity in %
I2_ml <- function(model, vi) {
  # Extract variance components (study-level and within-study)
  vc <- as.numeric(model$sigma2)
  total_het <- sum(vc)
  
  # Mean sampling variance
  mean_vi <- mean(vi, na.rm = TRUE)
  
  # Total variability
  total_var <- total_het + mean_vi
  
  # I2 in proportion
  I2 <- total_het / total_var
  
  # Return I2 as %
  return(I2 * 100)
}
I2_bio <- I2_ml(res_bio, biomass_mc_data$vi)
I2_bio

summary(res_bio)
bio_r
bio_ci
I2_bio




version

