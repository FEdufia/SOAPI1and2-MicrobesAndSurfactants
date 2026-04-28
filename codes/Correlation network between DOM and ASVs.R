library(readxl)
library(tidyverse)
library(igraph)
library(ggraph)

#relative abundance data
rel16 = read_xlsx("Data/rel16.xlsx") #16S data
rel18 = read_xlsx("Data/rel18.xlsx") #18S data
metdat = read_xlsx("Data/Metdata2.xlsx")

#selecting 16s top phyla, classes and orders---------------
top_phylum16 <- rel16 %>%
  group_by(no,  phylum) %>%
  summarise(agg_rel_abund=sum(relAb)) %>%
  arrange(desc(agg_rel_abund))|>
  group_by( phylum) %>%
  summarise(median=median(agg_rel_abund)) %>%
  arrange((desc(median))) %>% 
  top_n(11, median) %>%
  pull(phylum)

#selecting 18s top phyla, classes and orders---------------
top_Division18 <- rel18 %>%
  group_by(no,  Division) %>%
  summarise(agg_rel_abund=sum(relAb)) %>%
  arrange(desc(agg_rel_abund))|>
  group_by( Division) %>%
  summarise(median=mean(agg_rel_abund)) %>%
  arrange((desc(median))) %>% 
  top_n(6, median) %>%
  pull(Division)


#metadata for ASVs
phydat =  rel16|>
  select(otu, phylum)|>
  distinct()|>
  mutate(phylum = ifelse(phylum %in% top_phylum16, phylum, "Others"))

#16S surfactant significant OTUs ---------------------------------------
meta167 <- c("type", "no", "ord")

clr16_wide2 = rel16|>
  dplyr::select(no, ord, otu, type, relAb)|>
  pivot_wider(names_from = otu, values_from = relAb, values_fill = 0) |>
  (\(df) {
    otu_cols <- setdiff(names(df), meta167)
    clr_mat <- compositions::clr(df[, otu_cols] + 1) |> as.data.frame()
    bind_cols(df[, meta167], clr_mat)
  })()


ASV_DOM = metdat |>
  inner_join(clr16_wide2) |>
  mutate(
    Surfactant = (`An_612 nm` / `DOC ppm`),
    across(c(Comp1, Comp2, Comp3, Comp4),
           ~ . / `DOC ppm`)
  )

# ASVs
asv_columns <- ASV_DOM|> 
  dplyr::select(starts_with("ASV")) |> 
  names()

# FDOM + surfactant variables
dom_columns <- ASV_DOM |> 
  dplyr::select(
    Surfactant,
    starts_with("Comp")
  ) |> 
  names()


significant_results <- data.frame()

for (dom_var in dom_columns) {
  
  for (asv_var in asv_columns) {
    
    dom_values <- ASV_DOM[[dom_var]]
    asv_values <- ASV_DOM[[asv_var]]
    
    complete_pairs <- sum(complete.cases(dom_values, asv_values))
    
    if (complete_pairs > 3) {
      
      correlation_result <- cor.test(
        dom_values,
        asv_values,
        method = "pearson",
        use = "pairwise.complete.obs", exact = F
      )
      
      result_row <- data.frame(
        DOM_variable = dom_var,
        ASV = asv_var,
        CorrelationCoefficient = correlation_result$estimate,
        PValue = correlation_result$p.value
      )
      
      significant_results <- rbind(significant_results, result_row)
    }
  }
}

significant_otus = significant_results|>
  filter(PValue < 0.05)

significant_otus|>
  mutate(Correlation = ifelse(CorrelationCoefficient>0, "Positive", "Negative") )|>
  group_by(DOM_variable)|>
  summarise(n = n())

significant_otus|>
  mutate(Correlation = ifelse(CorrelationCoefficient>0, "Positive", "Negative") )|>
  group_by(DOM_variable, Correlation)|>
  summarise(n = n())

ASV_DOM_Corr16 =  rel16|>
  select(ASV = otu, genus, phylum)|>
  distinct()|>
  inner_join(significant_otus)|>
  mutate(Correlation = ifelse(CorrelationCoefficient>0, "Positive", "Negative"))

write_xlsx(ASV_DOM_Corr16, "Data/16S_ASVsCorrelatingWithDOM.xlsx")

#16S correlation_network -------------------------------------

edges <- significant_otus|>
  mutate(
    weight_signed = CorrelationCoefficient,   
    weight = abs(CorrelationCoefficient)      
  ) |>
  select(
    from = DOM_variable,
    to   = ASV,
    weight,
    weight_signed
  )


# DOM nodes (colored by name)
dom_nodes <- data.frame(
  name = unique(edges$from),
  Type = "DOM",
  phylum = unique(edges$from)  # DOM name used as color
)

# ASV nodes (colored by phylum)
asv_nodes <- data.frame(
  name = unique(edges$to),
  Type = "ASV"
) %>%
  left_join(phydat, by = c("name" = "otu"))

# Combine nodes
nodes <- bind_rows(dom_nodes, asv_nodes)

# Fill NAs (should only be DOM nodes) with their names
nodes$phylum[is.na(nodes$phylum)] <- nodes$name[is.na(nodes$phylum)]

net <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
E(net)$weight <- edges$weight 

dom_palette = c(
  "Surfactant" = "#6B8FA8",
  "Comp1"      = "#A8906E",
  "Comp2"      = "#E8D4B8",
  "Comp3"      = "#B8C8D8",
  "Comp4"      = "#969696"
)

category_colors <- c(
  "Pseudomonadota" = "#1f78b4",          
  "Cyanobacteriota" = "#33a02c",        
  "Bacteroidota" = "#ff7f00",           
  "Actinomycetota" = "#e31a1c",         
  "Verrucomicrobiota" = "#6a3d9a",      
  "Planctomycetota" = "#b15928",        
  "SAR324_clade(Marine_group_B)" = "#a6cee3",  
  "Marinimicrobia_(SAR406_clade)" = "#b2df8a", 
  "Bacillota" = "#fb9a99",              
  "Bdellovibrionota" = "#fdbf6f",       
  "Dadabacteria" ="yellow3",
  "Others" = "#999999"                  
)



# Combine palettes for all nodes
node_palette <- c(dom_palette, category_colors)

gg1 = ggraph(net, layout = "fr") +
  geom_edge_link(aes( width = NULL), color = "grey", alpha = 0.5) +
  geom_node_point(aes(
    color = phylum,                     # phylum or DOM name
    shape = Type,
    size = ifelse(Type == "DOM", 6, 2)  # bigger for DOM
  )) +
  scale_color_manual(values = node_palette) +
  scale_size_identity() +
  theme_void() +
  theme(legend.position = "right")#Legend was edited using inkscape

#### network stats 16s####
library(bipartite)

web_matrix <- edges %>%
  mutate(value = 1) %>%   
  select(from, to, value) %>%
  pivot_wider(
    names_from = to,
    values_from = value,
    values_fill = 0
  ) %>%
  column_to_rownames("from") %>%
  as.matrix()

networklevel(
  web_matrix,
  index="ALLBUTDD", 
  level = "both",
  weighted = FALSE,                
  ISAmethod = "Bluethgen",
  SAmethod = "Bluethgen",
  H2_integer = FALSE,              
  normalise = TRUE
)

shared_mat <- web_matrix %*% t(web_matrix)

significant_otus|>
  inner_join(phydat, by = join_by("ASV" == "otu"))|>
  group_by( phylum)|>
  summarise(n = n())



gg2 = significant_otus |>
  inner_join(phydat, by = join_by("ASV" == "otu")) |>
  mutate(Correlation = ifelse(CorrelationCoefficient > 0, "Positive Correlations", "Negative Correlations"))|>
  group_by(DOM_variable, phylum, Correlation) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(phylum = factor(phylum, levels = c(top_phylum16, "Others"))) |>
  mutate(DOM_variable = factor(DOM_variable, levels = c("Comp4", "Surfactant", "Comp1", "Comp2", "Comp3"), labels=c("Tyrosine-like\nC4", "Surfactant", "Tryptophan-like\nC1", "Humic-like\n C2", "Humic-like\nC3")))|>
  mutate(Correlation = factor(Correlation, levels = c("Positive Correlations", "Negative Correlations")))|>
  ggplot(aes(x = DOM_variable, # reorder(DOM_variable, -n, FUN = sum),
             y = n,
             fill = phylum)) +
  geom_col(color = "black") +
  theme_few(base_size = 14) +
  scale_fill_manual(values = category_colors) +facet_wrap(~Correlation, scales = "fixed", ncol = 1)+
  labs(x = "DOM variable", y = "Number of significant ASVs")+theme(legend.position = "none")




#18S Correlations with surfactants----------------------------------
phydat18 =  rel18|>
  select(otu, Division)|>
  distinct()|>
  mutate(Division = ifelse(Division %in% top_Division18, Division, "Others"))

#16S surfactant significant OTUs ---------------------------------------
meta168 <- c("type", "no", "ord")

clr18_wide2 = rel18|>
  dplyr::select(no, ord, otu, type, relAb)|>
  pivot_wider(names_from = otu, values_from = relAb, values_fill = 0) |>
  (\(df) {
    otu_cols <- setdiff(names(df), meta168)
    clr_mat <- compositions::clr(df[, otu_cols] + 1) |> as.data.frame()
    bind_cols(df[, meta168], clr_mat)
  })() 



ASV_DOM18 = metdat |>
  inner_join(clr18_wide2)|>
  mutate(
    Surfactant = (`An_612 nm` / `DOC ppm`),
    across(c(Comp1, Comp2, Comp3, Comp4),
           ~ . / `DOC ppm`)
  )

# ASVs
asv_columns18 <- ASV_DOM18  |> 
  dplyr::select(starts_with("ASV")) |> 
  names()

# DOM + surfactant variables
dom_columns18 <- ASV_DOM18  |> 
  dplyr::select(
    Surfactant,
    starts_with("Comp")
  ) |> 
  names()

significant_results18 <- data.frame()

for (dom_var in dom_columns18) {
  
  for (asv_var in asv_columns18) {
    
    dom_values <- ASV_DOM18 [[dom_var]]
    asv_values <- ASV_DOM18 [[asv_var]]
    
    complete_pairs <- sum(complete.cases(dom_values, asv_values))
    
    if (complete_pairs > 3) {
      
      correlation_result <- cor.test(
        dom_values,
        asv_values,
        method = "pearson",
        use = "pairwise.complete.obs", exact = F
      )
      
      result_row <- data.frame(
        DOM_variable = dom_var,
        ASV = asv_var,
        CorrelationCoefficient = correlation_result$estimate,
        PValue = correlation_result$p.value
      )
      
      significant_results18 <- rbind(significant_results18, result_row)
    }
  }
}

significant_otus18 = significant_results18|>
  filter(PValue < 0.05)
significant_otus18|>
  mutate(Correlation = ifelse(CorrelationCoefficient>0, "Positive", "Negative") )|>
  group_by(DOM_variable)|>
  summarise(n = n())

significant_otus18|>
  mutate(Correlation = ifelse(CorrelationCoefficient>0, "Positive", "Negative") )|>
  group_by(DOM_variable, Correlation)|>
  summarise(n = n())

ASV_DOM_Corr18 =  rel18|>
  select(ASV = otu, genus, Division)|>
  distinct()|>
  inner_join(significant_otus18)|>
  mutate(Correlation = ifelse(CorrelationCoefficient>0, "Positive", "Negative"))

write_xlsx(ASV_DOM_Corr18, "Data/18S_ASVsCorrelatingWithDOM.xlsx")

#18S correlation_network -------------------------------------

edges18 <- significant_otus18|>
  mutate(
    weight_signed = CorrelationCoefficient,   
    weight = abs(CorrelationCoefficient)      
  ) |>
  select(
    from = DOM_variable,
    to   = ASV,
    weight,
    weight_signed
  )



# DOM nodes (colored by name)
dom_nodes18 <- data.frame(
  name = unique(edges18$from),
  Type = "DOM",
  phylum = unique(edges18$from)  # DOM name used as color
)

# ASV nodes (colored by phylum)
asv_nodes18 <- data.frame(
  name = unique(edges18$to),
  Type = "ASV"
) %>%
  left_join(phydat18, by = c("name" = "otu"))

# Combine nodes
nodes18 <- bind_rows(dom_nodes18, asv_nodes18)

# Fill NAs (should only be DOM nodes) with their names
nodes18$Division[is.na(nodes18$Division)] <- nodes18$name[is.na(nodes18$Division)]

category3_colors <- c(
  "Alveolata" = "#8da0cb",             
  "Chlorophyta" = "#66c2a5",          
  "Stramenopiles" = "#e5c494",         
  "Eukaryota_unclassified" = "#FFC300", 
  "Picozoa" = "#e78ac3",               
  "Rhizaria" = "#ffd92f",              
  "Cryptophyta" = "#a6d854", 
  "Rhodophyta" = "#34495E",
  "Others" = "#999999"                  
)


net18 <- graph_from_data_frame(edges18, vertices = nodes18, directed = FALSE)
E(net18)$weight <- edges18$weight 

node_palette18 <- c(dom_palette, category3_colors)

gg3 = ggraph(net18, layout = "fr") +
  geom_edge_link(aes( width = NULL),color = "grey", alpha = 0.5) +
  # scale_edge_color_gradient2(
  #   high = "blue",
  #   mid = "grey80",
  #   low = "red",
  #   midpoint = 0,
  #   name = "Correlation"
  # ) +
  geom_node_point(aes(
    color = Division,                     # phylum or DOM name
    shape = Type,
    size = ifelse(Type == "DOM", 6, 2)  # bigger for DOM
  )) +
  scale_color_manual(values = node_palette18) +
  scale_size_identity() +
  theme_void() +
  theme(legend.position = "right")

#### network stats 18s####
library(bipartite)

web_matrix18 <- edges18 %>%
  mutate(value = 1) %>%   # presence/absence of significant correlation
  select(from, to, value) %>%
  pivot_wider(
    names_from = to,
    values_from = value,
    values_fill = 0
  ) %>%
  column_to_rownames("from") %>%
  as.matrix()


networklevel(
  web_matrix18,
  index="ALLBUTDD", 
  level = "both",
  weighted = FALSE,                
  ISAmethod = "Bluethgen",
  SAmethod = "Bluethgen",
  H2_integer = FALSE,             
  normalise = TRUE
)

shared_mat18 <- web_matrix18 %*% t(web_matrix18)


gg4 = significant_otus18|>
  inner_join(phydat18, by = join_by("ASV" == "otu"))|>
  mutate(Correlation = ifelse(CorrelationCoefficient > 0, "Positive Correlations", "Negative Correlations"))|>
  group_by(DOM_variable, Division, Correlation)|>
  summarise(n = n(), .groups = "drop") |>
  mutate(Correlation = factor(Correlation, levels = c("Positive Correlations", "Negative Correlations")))|>
  mutate(Division = factor(Division, levels = c(top_Division18, "Others"))) |>
  mutate(DOM_variable = factor(DOM_variable, levels = c("Comp4", "Surfactant", "Comp1", "Comp2", "Comp3"), labels=c("Tyrosine-like\nC4", "Surfactant", "Tryptophan-like\nC1", "Humic-like\n C2", "Humic-like\nC3")))|>
  mutate(Correlation = factor(Correlation, levels = c("Positive Correlations", "Negative Correlations")))|>
  ggplot(aes(x = reorder(DOM_variable, -n, FUN = sum),
             y = n,
             fill = Division)) +
  geom_col(color = "black") +
  theme_few(base_size = 14) +
  scale_fill_manual(values = category3_colors) + facet_wrap(~Correlation, scales = "fixed", ncol = 1)+
  labs(x = "DOM variable", y = "Number of significant ASVs")+theme(legend.position = "none")

#### figure 4 ####
(gg1 + gg2 ) +
  plot_layout(widths = c(2, 1.5)) +
  plot_annotation(
    tag_levels = 'a',
    tag_prefix = '(',
    tag_suffix = ')'
  ) &
  theme(plot.tag = element_text(face = "plain", size = 15))

#### figure 5 ####
(gg3 + gg4 ) +
  plot_layout(widths = c(2, 1.5)) +
  plot_annotation(
    tag_levels = 'a',
    tag_prefix = '(',
    tag_suffix = ')'
  ) &
  theme(plot.tag = element_text(face = "plain", size = 15))
