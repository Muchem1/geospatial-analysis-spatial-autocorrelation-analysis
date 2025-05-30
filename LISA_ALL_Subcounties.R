# Load required packages
library(sf)
library(spdep)
library(tidyverse)
library(tmap)


data <- read_csv("D:/Rabies_Analysis/rbs/points01_with_All_ADM2_EN_00.csv") %>%
  group_by(ADM2_EN) %>%
  summarise(Number_Sick = sum(Number_Sick, na.rm = TRUE)) %>%
  ungroup()

# View(data)
# str(data)

ADM2_EN_sf <- st_read("D:/DATA CENTER/ken_adm_iebc_20191031_shp/ken_admbnda_adm2_iebc_20191031.shp") %>%
  st_transform(4326) %>%
  select(ADM2_EN)


merged_sf <- ADM2_EN_sf %>%
  #  Left-join brings in Number_Sick, with NA where no match
  left_join(data, by = "ADM2_EN") %>%                      # :contentReference[oaicite:0]{index=0}
  #  Replace those NAs with zeros
  mutate(Number_Sick = replace_na(Number_Sick, 0))

invalid_idx <- which(!st_is_valid(merged_sf))
if (length(invalid_idx) > 0) {
  message("Found ", length(invalid_idx), " invalid geometries; attempting repair…")
  merged_sf[invalid_idx, ] <- st_make_valid(merged_sf[invalid_idx, ])
  # Double‐check
  stopifnot(all(st_is_valid(merged_sf)))
}

#Create spatial weights matrix (row‐standardized)
nb <- poly2nb(merged_sf, queen = TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

#  Calculate Global Moran's I
moran_results <- moran.test(merged_sf$Number_Sick, 
                            listw = lw, 
                            randomisation = TRUE)

# Print results
print(moran_results)

# Calculate Local Moran’s I
local_mi <- localmoran(merged_sf$Number_Sick, lw, zero.policy = TRUE)
local_mi <- bind_cols(data.frame(local_mi,check.names = F), attr(local_mi, "quadr"))
colnames(local_mi)
# View(local_mi)

# Attach local Moran’s I statistics back to the sf object
merged_sf <- merged_sf %>%
  mutate(
    Local_I       = local_mi[, "Ii"],     # the observed local Moran's I
    Expectation   = local_mi[, "E.Ii"],    # expected value under randomness
    Variance      = local_mi[, "Var.Ii"],  # variance of Ii
    Z_score       = local_mi[, "Z.Ii"],   # standardized Z‐score
    P_value       = local_mi[, "Pr(z != E(Ii))"],# p‐value
    quadr         = local_mi[, "mean"]
  )


# View(merged_sf)
# 
# st_write(
#   merged_sf,
#   dsn         = "D:/Rabies_Analysis/rbs/Final_Moran/subcounty_LISA.shp",
#   delete_layer = TRUE
# )
# Print the first few rows of the Local Moran’s I results
# print(head(merged_sf %>%
#              st_drop_geometry() %>%
#              select(ADM2_EN, Number_Sick, Local_I, Z_score, P_value)))


# # Save Local Moran’s I table
# write_csv(
#   merged_sf %>%
#     st_drop_geometry() %>%
#     select(ADM2_EN, Number_Sick, Z_score, Local_I, P_value, quadr),
#   "D:/Rabies_Analysis/rbs/Final_Moran/local_moran_results_ADM2_EN.csv"
# )
