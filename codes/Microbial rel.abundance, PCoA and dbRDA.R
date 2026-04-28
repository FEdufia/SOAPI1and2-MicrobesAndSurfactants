library(readxl)
library(tidyverse)
library(ggh4x)
library(ggpubr)
library(ggthemes)
library(patchwork)
library(cowplot)
library(vegan)
#relative abundance data
rel16 = read_xlsx("Data/rel16.xlsx") #16S data
rel18 = read_xlsx("Data/rel18.xlsx") #18S data


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



#relative abundance stack bars for 16s--------------
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


#### 16s stack bars for station ####
(f2a = rel16|>
   mutate(Time = factor(Time, levels = c("Sunrise", "Midday", "Sunset", "Night")))|>
   
   #filter(phylum %in% top_phylum16)|>
   mutate(phylum = ifelse((phylum %in% top_phylum16), phylum, "Others"))|>
   mutate(phylum = factor(phylum, levels = c(top_phylum16, "Others")))|>
   group_by(type, Day, station, Time, no, phylum)|>
   summarise(relAb = sum(relAb))|>
   group_by( type, station, phylum)|>
   summarise(relAb = mean(relAb))|>
   filter(relAb > 0)|>
   mutate(type = ifelse(type == "ML", "SML", "SSW"))|>
   mutate(type = factor(type, levels = c("SML", "SSW")))|>
   mutate(station = factor(station, levels = c("Virginia\nCoast","Delaware\nCoast Summer" ,"Continental\nSlope", "Delaware\nCoast Fall")))|>
   ggplot( aes(x = type, y = relAb, fill = phylum))+
   scale_fill_manual(values = category_colors)+
   scale_color_manual(values = category_colors)+
   geom_col(#stat = "identity",
     position = "fill", alpha = 0.8, color = "black")+facet_grid(~station, scales = "free_x")+
   theme_minimal_hgrid(font_size = 14)+labs(y = "Mean Relative Abundance", x = "Sample Type", fill = "Phylum" )+
   scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))
)


#### 18s stack bars for station ####
category3_colors <- c(
  "Alveolata" = "#8da0cb",             # Vibrant Orange-Red
  "Chlorophyta" = "#66c2a5",           # Bright Green
  "Stramenopiles" = "#e5c494",         # Bright Blue
  "Eukaryota_unclassified" = "#FFC300", # Golden Yellow
  "Picozoa" = "#e78ac3",               # Hot Pink
  "Rhizaria" = "#ffd92f",              # Rich Purple
  "Cryptophyta" = "#a6d854", 
  "Rhodophyta" = "#34495E",# Dark Slate Blue
  "Others" = "#999999"                  # Neutral Gray
)

(f2c = rel18|>
    mutate(Time = factor(Time, levels = c("Sunrise", "Midday", "Sunset", "Night")))|>
    #filter(phylum %in% top_phylum16)|>
    mutate(Division = ifelse((Division %in% top_Division18), Division, "Others"))|>
    mutate(Division = factor(Division, levels = c(top_Division18, "Others"
                                                  
    )))|>
    group_by(type, Day, station, Time, no, Division)|>
    summarise(relAb = sum(relAb))|>
    group_by( type, station, Division)|>
    summarise(relAb = mean(relAb))|>
    filter(relAb > 0)|>
    mutate(type = ifelse(type == "ML", "SML", "SSW"))|>
    mutate(type = factor(type, levels = c("SML", "SSW")))|>
    mutate(station = factor(station, levels = c("Virginia\nCoast","Delaware\nCoast Summer","Continental\nSlope", "Delaware\nCoast Fall" )))|>
    ggplot( aes(x = type, y = relAb, fill = Division))+
    scale_fill_manual(values = category3_colors)+
    geom_bar(stat = "identity", position = "fill", color = "black")+facet_grid(~station, scales = "free_x")+
    theme_minimal_hgrid(font_size = 14)+labs(y = "Mean Relative Abundance", x = "Sample Type", fill = "Division" )+scale_y_continuous(labels = scales::percent, expand = c(0.01, 0.01))# +theme(panel.spacing = unit(0.1, "lines"))+scale_x_discrete(expand = c(0, 0))
)




#PERMANOVA 16S-------------------------------------
# Prepare wide CLR matrix
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

# Distance matrix (Euclidean is appropriate for CLR-transformed data)
dist_mat <- dist(clr_mat, method = "euclidean")  # Aitchison distance

# Global PERMANOVA first
adonis2(dist_mat ~ station + type, data = clr_meta, permutations = 999, by = "term")


#clr transformed data-----------------------------------
####pcoa and dbrda data prep for 16s ####

metdat = read_xlsx("Data/Metdata2.xlsx")
meta166 <- c("type", "Day", 
             "Time", "no", "ord")



wrel16 = rel16 |>
  dplyr::select(-c( "phylum", "class", "order", "genus", "station"))|>
  pivot_wider(names_from = otu, values_from = relAb, values_fill = 0) |>
  (\(df) {
    otu_cols <- setdiff(names(df), meta166)
    clr_mat <- compositions::clr(df[, otu_cols] + 1) |> as.data.frame()
    bind_cols(df[, meta166], clr_mat)
  })() 

metdat = read_xlsx("Data/Metdata2.xlsx")

renv = metdat|>
  inner_join(wrel16)


rmenv = renv|>
  select(`DOC ppm`, Surfactant = `An_612 nm`, `Tryptophan-like C1` = Comp1, `Humic-like C2` = Comp2,
         `Humic-like C3` =  Comp3, `Tyrosine-like C4` =  Comp4, starts_with("ASV")) %>%
  mutate(across(-c(`DOC ppm`, starts_with("ASV")), ~ . / `DOC ppm`)) %>%
  select(-`DOC ppm`)


otudat = rmenv|>
  select(starts_with("ASV"))

DOM16 = rmenv|>
  select(
    Surfactant, 
    `Tryptophan-like C1`, `Humic-like C2` , 
    `Humic-like C3` , `Tyrosine-like C4`
  )

otumatrix = vegdist(otudat, method =  "euclidean")


#PCOA for 16S --------------------------------------

pcoa_result <- cmdscale(otumatrix, eig = T, k = 3)  
# Coordinates
pcoa_coords <- data.frame(pcoa_result$points)%>%cbind(., as.data.frame(renv)) |>
  mutate(station = case_when(station == "Virginia Coast" ~ "Virginia Coast",
                             station == "Delaware Coast" ~ "Delaware Coast Summer",
                             station == "Shelf Station" ~ "Continental Slope",
                             station == "Fall" ~ "Delaware Coast Fall"))|>
  mutate(station = factor(station, levels = c("Virginia Coast", "Delaware Coast Summer", "Continental Slope", "Delaware Coast Fall")))|>
  mutate(type = ifelse(type == "ML", "SML", "SSW"))
colnames(pcoa_coords)[1:3] <- c("PCoA1", "PCoA2","PCoA3" )  # Rename axes

# Variance explained
variance_explained <- pcoa_result$eig / sum(pcoa_result$eig)

permanova_label <- "PERMANOVA\nStation: R² = 0.600, p = 0.001\n Sample Type: R² = 0.036, p = 0.012"

(f2b = ggplot(pcoa_coords, aes(x = PCoA1, y = PCoA2, color = station, shape = type)) +
    geom_point(size = 4) +
    #geom_line(aes(group = no), color = "gray50", linetype = "dashed")+
    #theme_minimal() +
    labs(
      color = "Station", shape = "Sample Type",
      x = paste0("PCoA1 (", round(variance_explained[1] * 100, 2), "%)"),
      y = paste0("PCoA2 (", round(variance_explained[2] * 100, 2), "%)")
    )+scale_color_manual(values = c("#4daf4a", "#a65628" ,
                                    "#1919ff","#ffae19")
    )+
    geom_point(aes(color = station), alpha = 0.7, size =1) +
    geom_point(colour = "white", size = 1.5)+
    theme_few(base_size = 14)+
    annotate("text", x = 3, y = -0.5,
             label = permanova_label,
             hjust = 1.1, vjust = 1,
             size = 3.5, fontface = "italic")
)



#PERMANOVA 18S-------------------------------------
# Prepare wide CLR matrix
meta_cols <- c("type", "Day", "station", "Time", "no", "ord")
clr18_wide <- rel18 |>
  dplyr::select(-c( "Division", "species", "class", "order", "genus"))|>
  pivot_wider(names_from = otu, values_from = relAb, values_fill = 0) |>
  (\(df) {
    otu_cols <- setdiff(names(df), meta_cols)
    clr_mat <- compositions::clr(df[, otu_cols] + 1) |> as.data.frame()
    bind_cols(df[, meta_cols], clr_mat)
  })()

otu_cols18  <- setdiff(names(clr18_wide), meta_cols)



clr_mat18  <- clr18_wide[, otu_cols18]
clr_meta18 <- clr18_wide[, meta_cols]

# Distance matrix (Euclidean is appropriate for CLR-transformed data)
dist_mat18 <- dist(clr_mat18, method = "euclidean")  # Aitchison distance

# Global PERMANOVA first
adonis2(dist_mat18 ~ station + type, data = clr_meta18, permutations = 999, by = "term")




#### dataprep for 18S PCOA and dbRDA####

meta188 <- c("type", "Day", 
             "Time", "no", "ord")

wrel18 = rel18 |>
  dplyr::select(-c( "Division", "species", "class", "order", "genus", "station"))|>
  pivot_wider(names_from = otu, values_from = relAb, values_fill = 0) |>
  (\(df) {
    otu_cols <- setdiff(names(df), meta188)
    clr_mat <- compositions::clr(df[, otu_cols] + 1) |> as.data.frame()
    bind_cols(df[, meta188], clr_mat)
  })() 




r18env = metdat|>
  inner_join(wrel18)

r18menv = r18env|>
  select(`DOC ppm`, Surfactant = `An_612 nm`, `Tryptophan-like C1` = Comp1, `Humic-like C2` = Comp2,
         `Humic-like C3` =  Comp3, `Tyrosine-like C4` =  Comp4, starts_with("ASV")) %>%
  mutate(across(-c(`DOC ppm`, starts_with("ASV")), ~ . / `DOC ppm`)) %>%
  select(-`DOC ppm`)


otudat18= r18menv|>
  select(starts_with("ASV"))

DOM18 = r18menv|>
  select(-starts_with("ASV"))

otuenv18 = r18env
otumatrix18 = vegdist(otudat18, method = "euclidean")



# PCoA for 18S ------------------
pcoa_result18 <- cmdscale(otumatrix18, eig = T, k = 3)  
# Coordinates
pcoa_coords18 <- data.frame(pcoa_result18$points)%>%cbind(., as.data.frame(r18env)) |>
  mutate(station = case_when(station == "Virginia Coast" ~ "Virginia Coast",
                             station == "Delaware Coast" ~ "Delaware Coast Summer",
                             station == "Shelf Station" ~ "Continental Slope",
                             station == "Fall" ~ "Delaware Coast Fall"))|>
  mutate(station = factor(station, levels = c("Virginia Coast", "Delaware Coast Summer", "Continental Slope", "Delaware Coast Fall")))|>
  mutate(type = ifelse(type == "ML", "SML", "SSW"))

colnames(pcoa_coords18)[1:3] <- c("PCoA1", "PCoA2","PCoA3" )  # Rename axes

# Variance explained
variance_explained18 <- pcoa_result18$eig / sum(pcoa_result18$eig)

permanova_label18 <- "PERMANOVA\nStation: R² = 0.479, p = 0.001\n Sample Type: R² = 0.032, p = 0.046"

(f2d = ggplot(pcoa_coords18, aes(x = PCoA1, y = PCoA2, color = station, shape = type)) +
    geom_point(size = 4) +
    #geom_line(aes(group = ord), color = "gray50", linetype = "dashed")+
    #theme_minimal() +
    labs(
      color = "Station", shape = "Sample Type",
      x = paste0("PCoA1 (", round(variance_explained18[1] * 100, 2), "%)"),
      y = paste0("PCoA2 (", round(variance_explained18[2] * 100, 2), "%)")
    )+
    scale_color_manual(values = c("#4daf4a", "#a65628" , "#1919ff","#ffae19"))+
    #geom_point(aes(color = station), alpha = 0.7, size = 4) +
    geom_point(colour = "white", size = 1.5)+
    theme_few(base_size = 14)+
    annotate("text", x = 5, y = -2,
             label = permanova_label18,
             hjust = 1.1, vjust = 1,
             size = 3.5, fontface = "italic")
)


####FIGURE 2####


(f2a + f2b + f2c + f2d) +
  plot_annotation(
    tag_levels = 'a',
    tag_prefix = '(',
    tag_suffix = ')'
  ) &
  theme(plot.tag = element_text(face = "plain", size = 15))




#dbRDA for 16S----------------------------------------------------------------

rda1 = dbrda(otudat ~ . , data = DOM16, scale = T, na.action = na.exclude, distance = "euclidean")


# Model summaries
RsquareAdj(rda1) # Rsquared 
vif.cca(rda1) # variance inflation factor (<10 OK)
anova.cca(rda1, permutations = 999) # full model


R.sum <- summary(rda1)
#R.sum$cont   # Prints the "Importance of components" table
R.sum$cont$importance[2, "dbRDA1"]

R.sum$cont$importance[2, "dbRDA2"]


df1  <- data.frame(R.sum$sites[,1:3])%>%cbind(., (renv))|>
  mutate(type = ifelse(type == "ML", "SML", "SSW"))|>
  mutate(station = factor(station, levels = c("Virginia Coast", "Delaware Coast", "Shelf Station", "Fall"), labels = c("Virginia Coast", "Delaware Coast Summer", "Continental Slope", "Delaware Coast Fall")))  # PC1 and PC2

#rownames(df1)= df$Location
df2  <- data.frame(R.sum$species[,1:3])
df3 <- data.frame(R.sum$biplot[,1:3])


#### axis 1 vs axis 2 ############
(f3a = ggplot(df1, aes(x=dbRDA1, y=dbRDA2)) + 
   geom_point( aes( color = station, shape = type
   ), 
   size=3) +
   scale_color_manual(values = c("#4daf4a", "#a65628" ,
                                 "#1919ff","#ffae19")
   )+
   geom_point(aes(shape = type) , colour = "white",  size = 1)+
   geom_segment(data = df3, aes( x=0, xend=dbRDA1, y=0, yend=dbRDA2), 
                color="black", 
                arrow=arrow(length=unit(0.01,"npc")))+
   geom_text(data=df3, 
             aes(x=dbRDA1,y=dbRDA2,label=rownames(df3) 
                 , hjust=0.5*(1-sign(dbRDA1)),vjust=0.5*(1-sign(dbRDA2))
             ),color="black",
             size=4)+
   labs(x = paste0("dbRDA1 (", round(R.sum$cont$importance[2, "dbRDA1"] * 100, 2), "%)"),
        y = paste0("dbRDA2 (", round(R.sum$cont$importance[2, "dbRDA2"] * 100, 2), "%)"), color = "Station", shape = "Sample Type")+
   
   theme_few(base_size = 14)+
   geom_hline(yintercept=0, linetype="dotted") +
   geom_vline(xintercept=0, linetype="dotted")+theme_few(base_size = 15)+rremove("grid")+
   annotate("text", x = 0.5, y = -2,
            label = "dbRDA: df(5,23), R² = 0.32, p = 0.004",
            hjust = 1.1, vjust = 1,
            size = 3.5, fontface = "italic")
)


#dbRDA for 18S ---------------------------------------------------------------------

rda1p = dbrda(otudat18 ~ . , data = DOM18, scale = T, na.action = na.exclude, distance = "euclidean")
rda1p

R.sump <- summary(rda1p)
#R.sum$cont   # Prints the "Importance of components" table
R.sump$cont$importance[2, "dbRDA1"]
# 0.74785
R.sump$cont$importance[2, "dbRDA2"]

# Model summaries
RsquareAdj(rda1p) # adjusted Rsquared 
vif.cca(rda1p) # variance inflation factor (<10 OK)
anova.cca(rda1p, permutations = 999) # full model


R.sum18= summary(rda1p)

df1p  <- data.frame(R.sum18$sites[,1:3])%>%cbind(., (r18env))|>
  mutate(type = ifelse(type == "ML", "SML", "SSW"))|>
  mutate(station = factor(station, levels = c("Virginia Coast", "Delaware Coast", "Shelf Station", "Fall"), 
                          labels = c("Virginia Coast", "Delaware Coast Summer", "Continental Slope", "Delaware Coast Fall")))  # PC1 and PC2
#rownames(df1)= df$Location
df2p  <- data.frame(R.sum18$species[,1:3])
df3p <- data.frame(R.sum18$biplot[,1:3])


#### axis 1 vs axis 2 ############
(f3b = ggplot(df1p, aes(x=dbRDA1, y=dbRDA2)) + 
   #scale_color_manual(values =c("#7FC97F","#377EB8","#FDB462","#E31A1C") )+
   geom_point( aes( color = station, shape = type
   ), #colour="black", #pch=21,
   size=3) +
   scale_color_manual(values = c("#4daf4a", "#a65628" ,
                                 "#1919ff","#ffae19")
   )+
   #geom_point(aes(color = station), alpha = 0.7, size = 4) +
   geom_point(aes(shape = type) , colour = "white",  size = 1)+
   # stat_ellipse(geom = "polygon",
   #              aes(fill = station),
   #              alpha = 0.1, level = 0.80)+
   #scale_fill_manual(values =c("#7FC97F","#377EB8","#FDB462","#E31A1C") )+
   geom_segment(data = df3p, aes( x=0, xend=dbRDA1, y=0, yend=dbRDA2), 
                color="black", 
                arrow=arrow(length=unit(0.01,"npc")))+
   geom_text(data=df3p, 
             aes(x=dbRDA1,y=dbRDA2,label=rownames(df3p) 
                 , hjust=0.5*(1-sign(dbRDA1)),vjust=0.5*(1-sign(dbRDA2))
             ),color="black",
             size=4)+
   labs(x = paste0("dbRDA1 (", round(R.sum18$cont$importance[2, "dbRDA1"] * 100, 2), "%)"),
        y = paste0("dbRDA2 (", round(R.sum18$cont$importance[2, "dbRDA2"] * 100, 2), "%)"), color = "Station", shape = "Sample Type")+
   
   theme_few(base_size = 14)+
   # geom_point(data=df2,
   #            aes(x=dbRDA1,y=dbRDA2,label=rownames(df2)),color="#0970b9", size=4)+
   geom_hline(yintercept=0, linetype="dotted") +
   geom_vline(xintercept=0, linetype="dotted")+theme_few(base_size = 15)+rremove("grid")+
   annotate("text", x = -0, y = -1.4,
            label = "dbRDA: df(5,21), R² = 0.28, p = 0.012",
            hjust = 1.1, vjust = 1,
            size = 3.5, fontface = "italic")
)


#### FIGURE 3 ####
f3a +f3b+ plot_layout(guides = "collect")+
  plot_annotation(
    tag_levels = 'a',
    tag_prefix = '(',
    tag_suffix = ')'
  ) &
  theme(plot.tag = element_text(face = "plain", size = 15))

