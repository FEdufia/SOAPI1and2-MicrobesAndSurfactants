library(readxl)
library(tidyverse)
library(ggh4x)
library(ggpubr)
library(ggthemes)
library(patchwork)
library(cowplot)

metdat = read_xlsx("Data/Metdata2.xlsx")


#### figure 1 ####


# Prepare data once
dom_plot_dat <- metdat |>
  mutate(station = case_when(
    station == "Virginia Coast" ~ "Virginia\nCoast",
    station == "Delaware Coast" ~ "Delaware\nCoast Summer",
    station == "Shelf Station" ~ "Continental\nSlope",
    station == "Fall" ~ "Delaware\nCoast Fall")) |>
  mutate(station = factor(station, levels = c("Virginia\nCoast","Delaware\nCoast Summer",
                                              "Continental\nSlope","Delaware\nCoast Fall"))) |>
  select(station, type, Surfactant = `An_612 nm`, `DOC uM`, `DOC ppm`, Comp1, Comp2, Comp3, Comp4) |>
  mutate(type = ifelse(type == "ML", "SML", "SSW"),
         type = factor(type, levels = c("SML", "SSW"))) |>
  pivot_longer(names_to = "DOM", values_to = "values", -c(station, type, `DOC ppm`)) |>
  filter(DOM != "DOC uM") |>
  mutate(values = values / `DOC ppm`) |>
  mutate(DOM = factor(DOM,
                      levels = c("Surfactant", "Comp1", "Comp4", "Comp2", "Comp3"),
                      labels = c("Surfactant (µM/mg DOC)", "Tryptophan-like C1 (R.U/mg DOC)",
                                 "Tyrosine-like C4 (R.U/mg DOC)", "Humic-like C2 (R.U/mg DOC)",
                                 "Humic-like C3 (R.U/mg DOC)")))

# Shared plot function
make_dom_panel <- function(dom_var, letter, show_x = TRUE, show_legend = TRUE) {
  p <- dom_plot_dat |>
    filter(DOM == dom_var) |>
    ggboxplot(x = "station", y = "values", fill = "type") +
    theme_few(base_size = 15) +
    scale_fill_manual(values = c("skyblue", "grey80")) +
    labs(y = "", fill = "", x = "",
         title = paste0("(", letter, ") ", dom_var)) +
    theme(
      plot.title = element_text(size = 12, face = "plain"),
      legend.position = if (show_legend) "bottom" else "none",
      axis.text.x  = if (show_x) element_text() else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank()
    )
  return(p)
}

# Build individual panels
p_a <- make_dom_panel("Surfactant (µM/mg DOC)",          "a", show_x = TRUE,  show_legend = FALSE)
p_b <- make_dom_panel("Tryptophan-like C1 (R.U/mg DOC)", "b", show_x = FALSE, show_legend = FALSE)
p_c <- make_dom_panel("Tyrosine-like C4 (R.U/mg DOC)",   "c", show_x = FALSE, show_legend = FALSE)
p_d <- make_dom_panel("Humic-like C2 (R.U/mg DOC)",      "d", show_x = TRUE,  show_legend = FALSE)
p_e <- make_dom_panel("Humic-like C3 (R.U/mg DOC)",      "e", show_x = TRUE,  show_legend = FALSE)



legend_plot <- dom_plot_dat |>
  filter(DOM == "Surfactant (µM/mg DOC)") |>
  ggboxplot(x = "station", y = "values", fill = "type") +
  scale_fill_manual(values = c("skyblue", "grey80"), name = "") +
  theme_few() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.5, "cm")
  )

legend_grob <- get_legend(legend_plot)
legend_plot <- wrap_elements(legend_grob)
# Row 1: [a] [b] [c]
# Row 2: [legend] [d] [e]

layout <- "
ABC
LDE
"

p_a  + p_b + p_c + p_d + p_e +legend_plot+ 
  plot_layout(design = layout)


#### FIGURE S1 correlation matrix b/n DOM ####

dom_mat<- metdat|>
  select(Surfactant = `An_612 nm`, `DOC ppm`, 
         `Tryptophan-like C1` = Comp1, `Marine\nHumic-like C2` = Comp2, `Terestrial\nHumic-like C3` = Comp3, `Tyrosine-like C4` = Comp4) |>
  na.omit() %>%
  mutate(across(-`DOC ppm`, ~ . / `DOC ppm`)) %>%
  select(-`DOC ppm`)


library(Hmisc)

cor_res <- rcorr(as.matrix(dom_mat), type = "pearson")


cor_matrix <- cor_res$r
p_matrix  <- cor_res$P

sig_levels <- matrix("",
                     nrow = nrow(p_matrix),
                     ncol = ncol(p_matrix))

sig_levels[p_matrix < 0.05]  <- "*"
sig_levels[p_matrix < 0.01]  <- "**"
sig_levels[p_matrix < 0.001] <- "***"

library(corrplot)
corrplot(cor_matrix,
         method = "color",
         type = "upper",
         order = "hclust",
         #col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 1.2,
         p.mat = p_matrix,
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig")


