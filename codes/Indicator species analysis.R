library(readxl)
library(tidyverse)
library(indicspecies)
library(vegan)
library(writexl)

#relative abundance data
rel16 = read_xlsx("Data/rel16.xlsx") #16S data
rel18 = read_xlsx("Data/rel18.xlsx") #18S data


#16S SML vx SSW indicator species analysis--------------------------------------------


meta_cols <- c("type", "Day", "station", "Time", "no", "ord")

clr16_wide <- rel16 |>
  dplyr::select(-c( "phylum", "class", "order", "genus"))|>
  pivot_wider(names_from = otu, values_from = relAb, values_fill = 0) |>
  (\(df) {
    otu_cols <- setdiff(names(df), meta_cols)
    clr_mat <- compositions::clr(df[, otu_cols] + 1) |> as.data.frame()
    bind_cols(df[, meta_cols], clr_mat)
  })()

otu_cols  <- setdiff(names(clr16_wide), meta_cols)



clr_mat  <- clr16_wide[, otu_cols]
clr_meta <- clr16_wide[, meta_cols]

clr16 = clr16_wide|>
  pivot_longer(names_to = "otu", values_to = "relAb", starts_with("ASV"))|>
  inner_join(rel16|>select(-relAb))


perm_ctrl <- how(nperm = 999, 
                 blocks = clr_meta$station) 

indval_type <- multipatt(clr_mat,
                         cluster = clr_meta$type,
                         func = "r.g",
                         control = perm_ctrl)


r16_ind_results = as.data.frame(indval_type[["sign"]])|>
  rownames_to_column(var = "otu")|>
  filter(p.value<0.05)|>
  inner_join(rel16|>select(phylum, genus, otu)|>distinct(.keep_all = T))

#counting number of significant asvs
r16_ind_results|>
  mutate(type = ifelse(index == 1, "ML", "SS"))|>
  inner_join(rel16|>group_by(otu, phylum, type)|>
               summarise(relAb = mean(relAb)))|>
  group_by(type)|>
  summarise(asvs = n(),
            percent_count = (asvs*100)/1127, 
            relabund = sum(relAb))

EnrichedTaxa16 = rel16|>
  group_by(otu, phylum, type)|>
  summarise(relAb = mean(relAb))|>
  inner_join(r16_ind_results)|>
  inner_join(as.data.frame(indval_type[["str"]])|>
               rownames_to_column(var = "otu"))|>
  group_by(index)|>
  arrange(desc(relAb))|>
  #slice(1:30)|>
  mutate(genus = gsub("_insertae_sedis_genus", "", genus))|>
  mutate(genus = gsub("_unclassified", "", genus))|>
  mutate(ASV = gsub("000", "", otu))|>
  mutate(Taxa = paste0(genus, "_", ASV))|>   
  mutate(
    Habitat = factor(index,
                     levels = c(1, 2),
                     labels = c("SML-enriched", "SSW-enriched")),
    Taxa = factor(Taxa, levels = rev(unique(Taxa)))
  )|>pivot_longer(names_to = "Sample Type", values_to = "Stat", cols = c(ML, SS))|>
  filter(!type != `Sample Type`)|>
  mutate(p.value = ifelse(Stat < 0, 1, p.value),
         `Sample Type` = ifelse(`Sample Type` == "ML", "SML", "SSW"))|>
  filter(Stat > 0)|>
  ungroup()|>
  select(Phylum = phylum, Genus = genus, ASV, `Mean Abundance` = relAb, `Association Strength (r.g)` = Stat, `P value` = p.value, Habitat)


SMLEnrichedTaxa16 = EnrichedTaxa16|>
  filter(Habitat == "SML-enriched")

SSWEnrichedTaxa16 = EnrichedTaxa16|>
  filter(Habitat == "SSW-enriched")

write_xlsx(SMLEnrichedTaxa16, "Data/SMLEnrichedTaxa16.xlsx")
write_xlsx(SSWEnrichedTaxa16, "Data/SSWEnrichedTaxa16.xlsx")



#18S SML vx SSW indicator species analysis--------------------------------------------
meta_cols18 <- c("type", "Day", "station", "Time", "no", "ord")
clr18_wide <- rel18 |>
  dplyr::select(-c( "Division", "species", "class", "order", "genus"))|>
  pivot_wider(names_from = otu, values_from = relAb, values_fill = 0) |>
  (\(df) {
    otu_cols <- setdiff(names(df), meta_cols18)
    clr_mat <- compositions::clr(df[, otu_cols] + 1) |> as.data.frame()
    bind_cols(df[, meta_cols], clr_mat)
  })()

otu_cols18  <- setdiff(names(clr18_wide), meta_cols18)



clr_mat18  <- clr18_wide[, otu_cols18]
clr_meta18 <- clr18_wide[, meta_cols18]

clr18 = clr18_wide|>
  pivot_longer(names_to = "otu", values_to = "relAb", starts_with("ASV"))|>
  inner_join(rel18|>select(-relAb))

perm_ctrl18 <- how(nperm = 999, 
                   blocks = clr_meta18$station)  
indval_type18 <- multipatt(clr_mat18,
                           cluster = clr_meta18$type,
                           func = "r.g",
                           control = perm_ctrl18)



r18_ind_results = as.data.frame(indval_type18[["sign"]])|>
  rownames_to_column(var = "otu")|>
  filter(p.value<0.05)|>
  inner_join(rel18|>select(Division, genus, otu)|>distinct(.keep_all = T))

r18_ind_results|>
  mutate(type = ifelse(index == 1, "ML", "SS"))|>
  inner_join(rel18|>group_by(otu, Division, type)|>
               summarise(relAb = mean(relAb)))|>
  group_by(type)|>
  summarise(asvs = n(),
            percent_count = (asvs*100)/481, 
            relabund = sum(relAb))


EnrichedTaxa18 = rel18|>
  group_by(otu, Division, type)|>
  summarise(relAb = mean(relAb))|>
  inner_join(r18_ind_results)|>
  inner_join(as.data.frame(indval_type18[["str"]])|>
               rownames_to_column(var = "otu"))|>
  group_by(index)|>
  arrange(desc(relAb))|>
  #slice(1:30)|>
  mutate(genus = gsub("_insertae_sedis_genus", "", genus))|>
  mutate(genus = gsub("_unclassified", "", genus))|>
  mutate(ASV = gsub("000", "", otu))|>
  mutate(Taxa = paste0(genus, "_", ASV))|>   
  mutate(
    Habitat = factor(index,
                     levels = c(1, 2),
                     labels = c("SML-enriched", "SSW-enriched")),
    Taxa = factor(Taxa, levels = rev(unique(Taxa)))
  )|>pivot_longer(names_to = "Sample Type", values_to = "Stat", cols = c(ML, SS))|>
  filter(!type != `Sample Type`)|>
  mutate(p.value = ifelse(Stat < 0, 1, p.value),
         `Sample Type` = ifelse(`Sample Type` == "ML", "SML", "SSW"))|>
  filter(Stat > 0)|>
  ungroup()|>
  select(Division, Genus = genus, ASV, `Mean Abundance` = relAb, `Association Strength (r.g)` = Stat, `P value` = p.value, Habitat)


SMLEnrichedTaxa18 = EnrichedTaxa18|>
  filter(Habitat == "SML-enriched")

SSWEnrichedTaxa18 = EnrichedTaxa18|>
  filter(Habitat == "SSW-enriched")

write_xlsx(SMLEnrichedTaxa18, "Data/SMLEnrichedTaxa18.xlsx")
write_xlsx(SSWEnrichedTaxa18, "Data/SSWEnrichedTaxa18.xlsx")



#16S indicator species analysis for the different stations--------------------------------------------

perm_ctrls <- how(nperm = 999)  # equivalent to strata
indval_station <- multipatt(clr_mat,
                            cluster = clr_meta$station,
                            func = "r.g",
                            control = perm_ctrls)

r16_station_results = as.data.frame(indval_station[["sign"]])|>
  rownames_to_column(var = "otu")|>
  filter(p.value<0.05)|>
  inner_join(rel16|>select(phylum, genus, otu)|>distinct(.keep_all = T))

StationsTaxa16 = rel16|>
  # group_by(no, otu, phylum)|>
  # summarise(relAb = sum(relAb))|>
  group_by(otu, phylum)|>
  summarise(relAb = mean(relAb))|>
  inner_join(r16_station_results)|>
  inner_join(as.data.frame(indval_station[["str"]])|>
               rownames_to_column(var = "otu"))|>
  ungroup()|>
  select(1:15)|>
  group_by(index)|>
  arrange(desc(relAb))|>
  slice(1:30)|>
  mutate(genus = gsub("_insertae_sedis_genus", "", genus))|>
  mutate(genus = gsub("_unclassified", "", genus))|>
  mutate(ASV = gsub("000", "", otu))|>
  mutate(Taxa = paste0(genus, "_", otu))|>   # replace df with your object name
  pivot_longer(names_to = "Stations", values_to = "Stat", 
               cols = c(`Virginia\nCoast`,`Delaware\nCoast Summer`,`Delaware\nCoast Fall`, `Continental\nSlope`))|>
  mutate(Stations = factor(Stations, levels = c("Virginia\nCoast", "Delaware\nCoast Summer", "Delaware\nCoast Fall", "Continental\nSlope")))|>
  mutate(pstat = stat/Stat)|>
  filter(pstat == 1)|>
  filter(!is.na(p.value))|>
  ungroup()|>
  select(Phylum = phylum, Genus = genus, ASV, `Mean Abundance` = relAb, `Association Strength (r.g)` = Stat, `P value` = p.value, Stations)

write_xlsx(StationsTaxa16, "Data/StationsTaxa16.xlsx")  


#18S indicator species analysis for the different stations--------------------------------------------

indval_station18 <- multipatt(clr_mat18,
                              cluster = clr_meta18$station,
                              func = "r.g",
                              control = perm_ctrls)
summary(indval_station18, indvalcomp = TRUE, alpha = 0.05)

r18_station_results = as.data.frame(indval_station18[["sign"]])|>
  rownames_to_column(var = "otu")|>
  filter(p.value<0.05)|>
  inner_join(rel18|>select(Division, genus, otu)|>distinct(.keep_all = T))



StationsTaxa18 = rel18|>
  group_by(otu, Division)|>
  summarise(relAb = mean(relAb))|>
  inner_join(r18_station_results)|>
  inner_join(as.data.frame(indval_station18[["str"]])|>
               rownames_to_column(var = "otu"))|>
  ungroup()|>
  select(1:15)|>
  group_by(index)|>
  arrange(desc(relAb))|>
  slice(1:30)|>
  mutate(genus = gsub("_unclassified", "", genus))|>
  mutate(ASV = gsub("000", "", otu))|>
  mutate(Taxa = paste0(genus, "_", otu))|>   
  pivot_longer(names_to = "Stations", values_to = "Stat", 
               cols = c(`Virginia\nCoast`,`Delaware\nCoast Summer`,`Delaware\nCoast Fall`, `Continental\nSlope`))|>
  mutate(Stations = factor(Stations, levels = c("Virginia\nCoast", "Delaware\nCoast Summer", "Delaware\nCoast Fall", "Continental\nSlope")))|>
  mutate(pstat = stat/Stat)|>
  filter(pstat == 1)|>
  filter(!is.na(p.value))|>
  ungroup()|>
  select(Division, Genus = genus, ASV, `Mean Abundance` = relAb, `Association Strength (r.g)` = Stat, `P value` = p.value, Stations)

write_xlsx(StationsTaxa18, "Data/StationsTaxa18.xlsx")  


