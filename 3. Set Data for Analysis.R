############################
###Create analysis datasets
############################

#Load the cleaned dataset
library(dplyr)
library(metafor)
df <- read.csv("data/derived-data/Meta_analysis_cleaned_v2.csv")

##Confirm structre
head(df)

##Verify numerical column
str(df[, c("r_Microcystin", "Zr_Microcystin", "r_Biomass", "Zr_Biomass", "n")])

##Confirm the predictor category
unique(df$Predictor)

##Confirm missing rows
df %>% filter(Predictor == "") %>% head(20)

##Clean missing rows
df <- df %>% 
  filter(!(Predictor == "" & is.na(r_Microcystin) & is.na(r_Biomass)))
unique(df$Predictor)
nrow(df)
head(df)

##Create MC dataset for Predictors and MC
# All microcystin effect sizes with Temperature or Nitrate as predictor
micro_data <- df %>%
  filter(!is.na(r_Microcystin),
         !is.na(n),
         Predictor %in% c("Temperature", "Nitrate")) %>%
  mutate(
    yi = Zr_Microcystin,
    vi = 1/(n - 3)
  )
micro_data <- micro_data %>%
  filter(is.finite(vi), vi > 0)

###The above Drops any row without r_Microcystin or n
##Keeps only Temperature and Nitrate predictors
##Defines:
##yi = Fisher’s Z (effect size)
##vi = 1/(n − 3) (sampling variance), which metafor uses

##Confirm how many effect sizes is left
nrow(micro_data)
table(micro_data$Predictor)

#####################################################
##Split the effect sizes into Temperature and Nitrate
#####################################################
temp_data <- micro_data %>%
  filter(Predictor == "Temperature")

nit_data <- micro_data %>%
  filter(Predictor == "Nitrate")
nrow(temp_data)  # ~85
nrow(nit_data)   # ~60

##########Creating the data for Biomass to MC analysis
biomass_data <- df %>%
  filter(!is.na(r_Biomass),
         !is.na(n)) %>%
  mutate(
    yi = Zr_Biomass,
    vi = 1/(n - 3)
  )
biomass_data <- biomass_data %>%
  filter(is.finite(vi), vi > 0)

nrow(biomass_data)

############Confirm output of the above code
head(temp_data[, c("Study", "Predictor", "Temperature", "r_Microcystin", "yi", "vi")])

head(nit_data[, c("Study", "Predictor", "Nitrate.concentration", "r_Microcystin", "yi", "vi")])

head(biomass_data[, c("Study", "r_Biomass", "yi", "vi")])

##guarantees you never again reach rma.mv() with a broken V matrix.
stopifnot(
  all(is.finite(temp_data$vi)),
  all(is.finite(nit_data$vi)),
  all(is.finite(biomass_data$vi))
)

##Freeze ES identifiers as factors
temp_data$StudyID <- factor(temp_data$StudyID)
nit_data$StudyID  <- factor(nit_data$StudyID)
biomass_data$StudyID <- factor(biomass_data$StudyID)

unique(temp_data$Study)
