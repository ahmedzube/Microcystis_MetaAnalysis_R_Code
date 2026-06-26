library(readxl)
library(dplyr)
getwd()
dir()
library(readxl)
library(dplyr)

# 1. Load the raw file
df <- read_excel("data/raw-data/Meta analysis Spreadsheet paper.xlsx")

# 2. Function to fix Unicode minus signs
fix_minus <- function(x) {
  x <- gsub("\u2212", "-", x)   # Unicode minus U+2212
  x <- gsub("–", "-", x)        # En-dash (U+2013)
  return(x)
}

# 3. Clean all character columns
df_clean <- df %>%
  mutate(across(where(is.character), fix_minus))

# 4. Convert key columns to numeric
numeric_cols <- c("r_Microcystin", "Zr_Microcystin",
                  "r_Biomass", "Zr_Biomass", "n")

for (col in numeric_cols) {
  df_clean[[col]] <- suppressWarnings(as.numeric(df_clean[[col]]))
}

# 5. Save cleaned file in derived-data folder
write.csv(df_clean,
          "data/derived-data/Meta_analysis_cleaned_v2.csv",
          row.names = FALSE)

file.exists("data/derived-data/Meta_analysis_cleaned_v2.csv")
file.info("data/derived-data/Meta_analysis_cleaned_v2.csv")$size
df_clean <- read.csv("data/derived-data/Meta_analysis_cleaned_v2.csv")

