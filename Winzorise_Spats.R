# ---- Packages ----
library(readxl)
library(dplyr)
library(SpATS)

# ---- Load data ----
df <- read_excel("SAP2020_merged_v3.1.xls")

# ---- Winsorize helper ----
winsorize <- function(x, lower = 0.05, upper = 0.95) {
  q <- quantile(x, c(lower, upper), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

# ---- Winsorize all numeric traits ----
df_wins <- df %>% mutate(across(where(is.numeric), winsorize))

# ---- Save & reload ----
write.csv(df_wins, "2020_winsorized.csv", row.names = FALSE)
df_wins <- read.csv("2020_winsorized.csv")

# ---- One trait ----
trait <- "DaysToBloom"
df_wins$Row <- as.numeric(df_wins$Row)
df_wins$Column <- as.numeric(df_wins$Column)
df_wins$Genotype <- as.factor(df_wins$SorghumAccession)

# ---- Fixed model ----
m_fixed <- SpATS(trait,
                 genotype = "Genotype",
                 genotype.as.random = FALSE,
                 spatial = ~ SAP(Row, Column, nseg = c(33, 17)),  # keep nseg
                 data = df_wins)

# Extract BLUEs (genotype effect + intercept)
coeffs <- as.data.frame(m_fixed$coeff)
coeffs$Genotype <- rownames(coeffs)
intercept <- coeffs$value[coeffs$Genotype == "Intercept"]

BLUEs <- coeffs %>%
  filter(Genotype != "Intercept") %>%
  transmute(Genotype,
            !!paste0(trait, "_BLUE") := value + intercept)

# ---- Random model (heritability) ----
m_random <- SpATS(trait,
                  genotype = "Genotype",
                  genotype.as.random = TRUE,
                  spatial = ~ SAP(Row, Column, nseg = c(33, 17)), # keep nseg
                  data = df_wins)
#for nseg Row do it = length(unique$Row) what ever number you get divide it by 2 and add 1 to it. Same goes to the column. THose are the number you will have to add to the row and column section. 
h2 <- getHeritability(m_random)
print(h2)

# ---- Plot spatial trend ----
plot(m_fixed, which = 1, main = paste("Spatial Trend:", trait))
