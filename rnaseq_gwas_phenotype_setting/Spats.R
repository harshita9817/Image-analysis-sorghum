mock_winzo<- read.csv("raw_data_treatment_m_winsorized.csv")
#--------------------------------------------------------#
#  FULL SPATS ANALYSIS FOR INOCULATED DATASET (WIDE BLUEs)
#--------------------------------------------------------#
#--------------------------------------------------------#
#  FULL SPATS ANALYSIS FOR INOCULATED DATASET (WIDE BLUEs)
#--------------------------------------------------------#

library(SpATS)
library(dplyr)
library(tidyr)

#--------------------------------------------------------#
# 1. Prepare dataset
#--------------------------------------------------------#
df <- inoculated_winzo

# Clean and convert types
df$Row <- as.numeric(df$Row)
df$Column <- as.numeric(df$Column)
df$Genotype <- as.factor(df$PINumber)   # Genotype ID
df$Repeat <- as.factor(df$Repeat)
df$Treatment <- as.factor(df$Treatment) # 'i' for inoculated

#--------------------------------------------------------#
# 2. Define traits to analyze
#--------------------------------------------------------#
traits <- c("percent_pigment",
            "percent_chlorosis",
            "percent_necrosis",
            "percent_all_response")

# Compute spatial grid segments
nseg_row <- round(length(unique(df$Row)) / 2 + 1)
nseg_col <- round(length(unique(df$Column)) / 2 + 1)

# Make output folders
dir.create("SpATS_outputs", showWarnings = FALSE)
dir.create("SpATS_outputs/plots", showWarnings = FALSE)

# Initialize storage
all_BLUEs_long <- data.frame()
all_H2 <- data.frame()

#--------------------------------------------------------#
# 3. Loop through all traits
#--------------------------------------------------------#
for (trait in traits) {
  cat("\n----------------------------------------\n")
  cat("Analyzing:", trait, "...\n")
  
  # Fixed model (BLUEs)
  m_fixed <- SpATS(
    response = trait,
    genotype = "Genotype",
    genotype.as.random = FALSE,
    fixed = ~ Repeat,
    spatial = ~ SAP(Row, Column, nseg = c(nseg_row, nseg_col)),
    data = df
  )
  
  # Extract coefficients
  coeffs <- m_fixed$coeff
  coeffs_df <- data.frame(
    Genotype = names(coeffs),
    value = as.numeric(coeffs),
    stringsAsFactors = FALSE
  )
  
  # Compute BLUEs
  intercept <- coeffs_df$value[coeffs_df$Genotype == "Intercept"]
  
  BLUEs <- coeffs_df %>%
    filter(Genotype != "Intercept") %>%
    transmute(
      Genotype,
      Trait = trait,
      BLUE = value + intercept
    )
  
  all_BLUEs_long <- rbind(all_BLUEs_long, BLUEs)
  
  # Random model (heritability)
  m_random <- SpATS(
    response = trait,
    genotype = "Genotype",
    genotype.as.random = TRUE,
    fixed = ~ Repeat,
    spatial = ~ SAP(Row, Column, nseg = c(nseg_row, nseg_col)),
    data = df
  )
  
  h2 <- getHeritability(m_random)
  all_H2 <- rbind(all_H2, data.frame(Trait = trait, Heritability = h2))
  
  # Save spatial trend plot
  png_filename <- paste0("SpATS_outputs/plots/", trait, "_SpatialTrend.png")
  png(png_filename, width = 1000, height = 800, res = 150)
  plot(m_fixed, which = 1, main = paste("Spatial Trend:", trait))
  dev.off()
  
  cat("Saved plot:", png_filename, "\n")
}

#--------------------------------------------------------#
# 4. Convert BLUEs to wide format (traits = columns)
#--------------------------------------------------------#
all_BLUEs_wide <- all_BLUEs_long %>%
  pivot_wider(names_from = Trait, values_from = BLUE)

# Merge with unique genotype metadata (PINumber & Treatment)
meta_df <- df %>%
  distinct(Genotype, PINumber, Treatment)

all_BLUEs_wide <- meta_df %>%
  left_join(all_BLUEs_wide, by = "Genotype")

#--------------------------------------------------------#
# 5. Save outputs
#--------------------------------------------------------#
write.csv(all_BLUEs_wide, "SpATS_outputs/Inoculated_BLUEs_WIDE.csv", row.names = FALSE)
write.csv(all_H2, "SpATS_outputs/Inoculated_Heritabilities.csv", row.names = FALSE)

cat("\n✅ DONE! Results saved in 'SpATS_outputs/' folder.\n")
cat("Files created:\n")
cat(" - Inoculated_BLUEs_WIDE.csv (one row per genotype)\n")
cat(" - Inoculated_Heritabilities.csv (per trait)\n")
cat(" - plots/ (spatial trend maps)\n")
