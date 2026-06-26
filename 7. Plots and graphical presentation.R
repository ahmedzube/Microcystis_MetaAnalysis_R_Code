list.files(recursive = TRUE)
files <- list.files("Code", full.names = TRUE)
file.edit(files)
##Representation of Result###
##(Context Table)##
table1 <- data.frame(
  Metric = c("Number of studies",
             "Total effect sizes",
             "Experimental setting",
             "Drivers tested",
             "Response variables",
             "Microcystin compartments"),
  Value = c(18,
            143,
            "Laboratory",
            "Temperature, Nitrate",
            "Biomass, Microcystin",
            "Intracellular, Extracellular, Total")
)
res_temp_biomass$k
res_nit_biomass$k
res_temp$k
res_nit_base$k
res_biomass_mc$k
cat(
  "Temp - Biomass:", res_temp_biomass$k, "\n",
  "Nitrate - Biomass:", res_nit_biomass$k, "\n",
  "Temp - MC:", res_temp$k, "\n",
  "Nitrate - MC:", res_nit_base$k, "\n",
  "Biomass - MC:", res_biomass_mc$k, "\n"
)
##Summary Meta-analytic Results (CORE TABLE)##
table2 <- data.frame(
  Relationship = c("Temperature → Biomass",
                   "Nitrate → Biomass",
                   "Temperature → Microcystin",
                   "Nitrate → Microcystin",
                   "Biomass → Microcystin"),
  k = c(37, 58, 83, 58, 141),
  r = c(0.81, 0.95, 0.39, 0.68, 0.61),
  CI = c("0.02–0.98",
         "0.89–0.98",
         "−0.34–0.82",
         "0.21–0.89",
         "0.20–0.84"),
  tau2 = c(1.87, 0.55, 1.33, 1.15, 1.26)
)


##FIGURE 2 — Drivers → Biomass → MC (SYNTHESIS PLOT)##
effects <- data.frame(
  Relationship = c("Temp → Biomass",
                   "Nitrate → Biomass",
                   "Temp → MC",
                   "Nitrate → MC",
                   "Biomass → MC"),
  r = c(0.81, 0.95, 0.39, 0.68, 0.61),
  lower = c(0.02, 0.89, -0.34, 0.21, 0.20),
  upper = c(0.98, 0.98, 0.82, 0.89, 0.84)
)

library(ggplot2)
ggplot(effects, aes(x = r, y = Relationship)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  xlab("Correlation coefficient (Pearson’s r)")
ggsave("Figure2_driver_biomass_MC.png", width = 6, height = 4, dpi = 300)

##FIGURE 3 — Compartment-specific Biomass → MC Effects (KEY FIGURE)##
compartment <- data.frame(
  Compartment = c("Intracellular", "Extracellular", "Total"),
  r = c(0.38, 0.66, 0.77),
  lower = c(-0.14, 0.24, 0.37),
  upper = c(0.73, 0.87, 0.93)
)

ggplot(compartment, aes(x = r, y = Compartment)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  xlab("Correlation coefficient (Pearson’s r)")
ggsave("Figure3_compartment_biomass_MC.png", width = 6, height = 4, dpi = 300)

library(metafor)

##Forest Plot Nitrate-Microcystin##
nit_data <- nit_data[order(nit_data$StudyID), ]
res_nit_base <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = nit_data,
  method = "REML"
)
nit_data$StudyLabel <- gsub(",", "", nit_data$Study)
study_breaks <- cumsum(table(nit_data$StudyID))
tiff(
  "Figure3_Nitrate_MC.tiff",
  width = 10,
  height = 16,
  units = "in",
  res = 600,
  compression = "lzw"
)
############plot 
nit_data$StudyLabel <- gsub(",", "", nit_data$Study)

pdf(
  "Figure3_Nitrate_MC.pdf",
  width = 10,
  height = 16
)

forest(
  res_nit_base,
  transf = tanh,
  refline = 0,
  xlab = "Pearson correlation coefficient (r)",
  slab = nit_data$StudyLabel,
  annotate = FALSE,
  cex = 0.8,
  psize = 1.2,
  header = "Study"
)

dev.off()

tiff(
  "Figure3_Nitrate_MC.tiff",
  width = 10,
  height = 16,
  units = "in",
  res = 600,
  compression = "lzw"
)

forest(
  res_nit_base,
  transf = tanh,
  refline = 0,
  xlab = "Pearson correlation coefficient (r)",
  slab = nit_data$StudyLabel,
  annotate = FALSE,
  cex = 0.8,
  psize = 1.2,
  header = "Study"
)

dev.off()

png(
  "Figure3_Nitrate_MC.png",
  width = 4000,
  height = 6500,
  res = 600
)

forest(
  res_nit_base,
  transf = tanh,
  refline = 0,
  xlab = "Pearson correlation coefficient (r)",
  slab = nit_data$StudyLabel,
  annotate = FALSE,
  cex = 0.8,
  psize = 1.2,
  header = "Study"
)

dev.off()
##Forest plot — Biomass → Microcystin##
# 1. Order data
micro_data <- micro_data[order(micro_data$StudyID), ]

# 2. Fit model
res_all_mc <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = micro_data,
  method = "REML"
)

# 3. Clean labels
micro_data$StudyLabel <- gsub(",", "", micro_data$Study)

# 4. Slab labels


# 5. Plot
micro_data$StudyLabel <- gsub(",", "", micro_data$Study)

pdf(
  "Figure4_Biomass_MC.pdf",
  width = 10,
  height = 24
)

forest(
  res_all_mc,
  transf = tanh,
  refline = 0,
  xlab = "Pearson correlation coefficient (r)",
  slab = micro_data$StudyLabel,
  annotate = FALSE,
  cex = 0.8,
  psize = 1.2,
  header = "Study"
)

dev.off()

tiff(
  "Figure4_Biomass_MC.tiff",
  width = 10,
  height = 24,
  units = "in",
  res = 600,
  compression = "lzw"
)

forest(
  res_all_mc,
  transf = tanh,
  refline = 0,
  xlab = "Pearson correlation coefficient (r)",
  slab = micro_data$StudyLabel,
  annotate = FALSE,
  cex = 0.8,
  psize = 1.2,
  header = "Study"
)

dev.off()

png(
  "FigureS1_Biomass_MC_FullForest.png",
  width = 5000,
  height = 12000,
  res = 600
)

forest(
  res_all_mc,
  transf = tanh,
  refline = 0,
  xlab = "Pearson correlation coefficient (r)",
  slab = micro_data$StudyLabel,
  annotate = FALSE,
  cex = 0.6,
  psize = 1,
  header = "Study"
)

dev.off()

######Biomass Pool effect sizes within each study
library(metafor)
library(dplyr)

study_pooled <- micro_data %>%
  group_by(StudyID, Study) %>%
  summarise(
    yi = weighted.mean(yi, 1/vi),
    vi = mean(vi),
    .groups = "drop"
  )
nrow(study_pooled)
res_study <- rma(
  yi,
  vi,
  data = study_pooled,
  method = "REML"
)
pdf(
  "Figure4_Biomass_MC_StudyLevel.pdf",
  width = 8,
  height = 8
)

forest(
  res_study,
  transf = tanh,
  slab = study_pooled$Study,
  annotate = FALSE,
  refline = 0,
  cex = 0.9,
  psize = 1.5,
  xlab = "Pearson correlation coefficient (r)",
  header = "Study"
)

dev.off()

tiff(
  "Figure4_Biomass_MC_StudyLevel.tiff",
  width = 8,
  height = 8,
  units = "in",
  res = 600,
  compression = "lzw"
)

forest(
  res_study,
  transf = tanh,
  slab = study_pooled$Study,
  annotate = FALSE,
  refline = 0,
  cex = 0.9,
  psize = 1.5,
  xlab = "Pearson correlation coefficient (r)",
  header = "Study"
)

dev.off()

png(
  "Figure4_Biomass_MC_StudyLevel.png",
  width = 3500,
  height = 3500,
  res = 600
)

forest(
  res_study,
  transf = tanh,
  refline = 0,
  xlab = "Pearson correlation coefficient (r)",
  slab = study_pooled$Study,
  annotate = FALSE,
  cex = 0.9,
  psize = 1.5,
  header = "Study"
)

dev.off()

##Forest plot — Temperature → Microcystin##
# 1. Order data
temp_data <- temp_data[order(temp_data$StudyID), ]

# 2. Refit model
res_temp <- rma.mv(
  yi, vi,
  random = ~ 1 | StudyID,
  data   = temp_data,
  method = "REML"
)

# 3. Clean labels
temp_data$StudyLabel <- gsub(",", "", temp_data$Study)

# 4. Slab labels
#temp_data$StudyLabel <- gsub(",", "", temp_data$Study)
#slab = temp_data$StudyLabel

# 5. Breaks
study_breaks <- cumsum(table(temp_data$StudyID))

# 6. Plot
##############################################################
# Figure 2
# Temperature → Microcystin Forest Plot
##############################################################

# Clean study names only
temp_data$StudyLabel <- gsub(",", "", temp_data$Study)

##############################################################
# PDF version
##############################################################

pdf(
  "Figure2_Temperature_MC.pdf",
  width = 10,
  height = 16
)

forest(
  res_temp,
  transf = tanh,
  refline = 0,
  xlab = "Pearson correlation coefficient (r)",
  slab = temp_data$StudyLabel,
  annotate = FALSE,
  cex = 0.8,
  psize = 1.2,
  header = "Study"
)

dev.off()

##############################################################
# TIFF version
##############################################################

tiff(
  "Figure2_Temperature_MC.tiff",
  width = 10,
  height = 16,
  units = "in",
  res = 600,
  compression = "lzw"
)

forest(
  res_temp,
  transf = tanh,
  refline = 0,
  xlab = "Pearson correlation coefficient (r)",
  slab = temp_data$StudyLabel,
  annotate = FALSE,
  cex = 0.8,
  psize = 1.2,
  header = "Study"
)

dev.off()
png(
  "Figure2_Temperature_MC.png",
  width = 4000,
  height = 7500,
  res = 600
)

forest(
  res_temp,
  transf = tanh,
  refline = 0,
  xlab = "Pearson correlation coefficient (r)",
  slab = temp_data$StudyLabel,
  annotate = FALSE,
  cex = 0.8,
  psize = 1.2,
  header = "Study"
)

dev.off()

##Regtest and Funnel — Nitrate → Microcystin
res_nit_simple <- rma(
  yi, vi,
  data = nit_data,
  method = "REML"
)

regtest(res_nit_simple)
png("funnel_nitrate_MC.png", width = 1800, height = 1400, res = 600)
funnel(res_nit_base, atransf = tanh, xlab = "Correlation coefficient (Pearson’s r)")
dev.off()

##Regtest and Funnel — Biomass → Microcystin
res_biomass_simple <- rma(
  yi, vi,
  data = biomass_data,
  method = "REML"
)
regtest(res_biomass_simple)
png("funnel_biomass_MC.png", width = 1800, height = 1400, res = 600)
funnel(res_biomass_mc_base, atransf = tanh, xlab = "Correlation coefficient (Pearson’s r)")
dev.off()

##Regtest and Funnel — Temp → Microcystin
res_temp_simple <- rma(
  yi, vi,
  data = temp_data,
  method = "REML"
)
regtest(res_temp_simple)
png("funnel_temp_MC.png", width = 1800, height = 1400, res = 600)
funnel(res_temp, atransf = tanh, xlab = "Correlation coefficient (Pearson’s r)")
dev.off()

getwd()