# Produces graphics for figures in Derivatization Manuscript

# load libraries ----
library(ggplot2)
library(Rmisc)
library(RColorBrewer)
library(ggh4x)
library(scales)
library(ggprism)
library(dplyr)
library(ggbreak)
library(fillpattern)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(forcats)
library(ggVennDiagram)

# Formatting and colors -----
axis_title_size <- 11
axis_text_size <- 11

# colors
gentle_rainbow <- c("#FF595E", "#FFCA3A", "#8AC926", 
                    "#118AB2", "#0B4982", "#6A4C93")
color_set <- rev(gentle_rainbow)
color_pair <- c(color_set[1], "#706b70")
show_col(color_pair)
color_bars <- c(color_set[1], color_set[4], "#706b70")
show_col(color_bars)
chrom_colors <- c("#6A4C93", "#9480b4", "#c0b5d4", 
                  "#706b70", "#999299", "#ada7ad")

chrom_colors_purple <- colorRampPalette(c("#6A4C93","#c0b5d4"))
chrom_colors_grey <- colorRampPalette(c("#494949","#CACACA"))
chrom_colors_long <- c(chrom_colors_purple(5),
                       chrom_colors_grey(5))
show_col(color_set)
lty_pal <- c("solid", "dotted", "dashed",
             "solid", "dotted", "dashed")



# Main text figures
# Figure 1: Ionization Efficiencies -----
# Panel B ---
ionizationData_1B <- read.csv("Figure1/Figure1B_alkyl_IE.csv",
                           header = TRUE)
# calculate pmoles (2µL injections)
ionizationData_1B <- ionizationData_1B %>%
  mutate(Amount = Concentration * 2)

# sum peak areas (2nd, and 3rd charge state + sodium adducts)
ionizationData_1B <- ionizationData_1B %>%
  mutate(Area = Area..2 + Area..2Na + Area..3)

# remove blanks
ionizationData_1B <- subset(ionizationData_1B, Sample != "water")
ionizationData_1B <- subset(ionizationData_1B, !(is.na(Concentration)))

# perform linear regression
unique_samples <- unique(ionizationData_1B$Sample)
models <- list()
for (s in unique_samples) {
  subset_data <- ionizationData_1B[ionizationData_1B$Sample == s, ]
  models[[s]] <- lm(Area ~ Amount, data = subset_data)
}
slopes <- sapply(models, function(m) coef(m)["Amount"])
sorted_samples <- names(sort(slopes, decreasing = TRUE))
sorted_samples <- sub("\\.Amount$", "", sorted_samples)
ionizationData_1B$Sample <- factor(ionizationData_1B$Sample, levels = sorted_samples)
ionizationData_1B <- ionizationData_1B[order(ionizationData_1B$Sample), ]
sorted_slopes <- sort(slopes, decreasing = TRUE)

# get list of R squared values
r_squared_vals <- list()
for (s in unique_samples) {
  r_squared_vals[[s]] <- summary(models[[s]])$adj.r.squared
}

# summarize for figure generation
data1B_summary <- summarySE(ionizationData_1B, 
                        measurevar = "Area",
                        groupvars = c("Amount", "Sample"))

makeIEplot_1B <- ggplot(data = data1B_summary,
                      aes(x = Amount,
                          y = Area,
                          color = Sample,
                          fill = Sample)) +
  geom_errorbar(aes(ymin = Area-se, ymax = Area+se), 
                width = 2.5) +
  geom_point(size = 0.5) +
  geom_smooth(method='lm', formula=y~x, se=FALSE, linewidth = 0.5) +
  scale_fill_manual(values = color_set)+
  scale_color_manual(values = color_set)+
  scale_y_continuous(expand = c(0,0), labels = function(x) x / 1e5,
                     limits = c(0, 500000))+
  scale_x_continuous(guide = "prism_minor",
                     minor_breaks = seq (0, 70, 4),
                     expand = c(0,0),
                     limits = c(0,70))+
  guides(y = "prism_minor")+
  xlab("Amount (pmoles)") +
  ylab(expression("Peak Area *"*10^5*"")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.minor.ticks.length = rel(0.5),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        # legend.position = "inside",
        # legend.position.inside = c(0.15, 0.7),
        )

makeIEplot_1B

ggsave("plots/Figure1B.pdf",
       width = 2.16, height = 2, unit = "in")
# Panel C ---

chromatogramData_1C <- read.csv("Figure1/Figure1C_chromatogram.csv",
                                header = TRUE)

makeChromatogram_1C <- ggplot(data = chromatogramData_1C,
                        aes(x = X.Minutes.,
                            y = Y.Counts. / 10000)) +
  facet_wrap2(~ fct_rev(Analyte),
              axes = "x",
              remove_labels = "x",
             ncol = 1) +
  geom_line(aes(color = fct_rev(Analyte))) +
  geom_ribbon(alpha = 0.6, aes(ymin = 0, ymax = Y.Counts. / 10000, 
                  fill = fct_rev(Analyte))) + 
  scale_color_manual(values = color_pair) +
  scale_fill_manual(values = color_pair) +
  xlab("Time (minutes)") +
  ylab(expression("Intensity *"*10^6*"")) +
  scale_x_continuous(breaks = seq(9,15, by=1), limits=c(9,15),
                     expand = c(0,0.2)) +
  scale_y_continuous(expand = c(0,0.05)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y=unit(0.5, "lines"),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.text = element_text(color = "black",
                                   size = axis_text_size),
        # legend.key = element_blank(),
        legend.key.size = unit(0.1, 'in'),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.6, 0.25),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

makeChromatogram_1C

ggsave("plots/Figure1C.pdf", dpi = 300,
       width = 2.16, height = 2, unit = "in")

# Panel D ---
ligationLOD_data1D <- read.csv("Figure1/Figure1D_Ligation_LOD.csv",
                               header = TRUE)
data1D_summary <- summarySE(ligationLOD_data1D, 
                            measurevar = "Peak.Area",
                            groupvars = c("Amount.Injected", "Analyte"))
# linear models to be used to calculate slope of cal curves
# this can be used as a readout for sensitivity
mod_model <- lm(Peak.Area ~ Amount.Injected, data = subset(ligationLOD_data1D, 
                                            Analyte == "Ligation Product"))
std_model <-lm(Peak.Area ~ Amount.Injected, data = subset(ligationLOD_data1D, 
                                            Analyte == "7mer"))
mod_m <- mod_model$coefficients[[2]]
std_m <- std_model$coefficients[[2]]
# fold change in slope
mod_m / std_m

makeIEplot_1D <- ggplot(data = data1D_summary,
                        aes(x = Amount.Injected,
                            y = Peak.Area,
                            color = Analyte,
                            fill = Analyte)) +
  geom_errorbar(aes(ymin = Peak.Area-se, ymax = Peak.Area+se), 
                width = 2, position=position_dodge(0)) +
  geom_point(size = 0.5) +
  geom_smooth(method='lm', formula=y~x, se=FALSE, linewidth = 0.5) +
  scale_fill_manual(values = rev(color_pair)) +
  scale_color_manual(values = rev(color_pair)) +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 900000, 300000),
                     minor_breaks = seq(0, 900000, 100000),
                     limits = c(0, 900000),
                     guide = guide_axis(minor.ticks = TRUE),
                     labels = function(x) x / 1e5) + 
  scale_x_continuous(minor_breaks = seq (0, 40, 2),
                     expand = c(0,0),
                     guide = guide_axis(minor.ticks = TRUE),
                     limits = c(0,42))+
  xlab("Amount (pmoles)") +
  ylab(expression("Peak Area *"*10^5*"")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.key = element_blank(),
        legend.position = "none")

makeIEplot_1D

ggsave("plots/Figure1D.pdf", dpi = 300,
       width = 2.16, height = 2, unit = "in")



# Figure 2: Ligation to standards ----

# 2B Bar graph comparing peak areas
mix_ligation_2B <- read.csv("Figure2/Figure2B_Ligation_Areas.csv")
data2B_summary <- summarySE(mix_ligation_2B, 
                             measurevar = "Peak.Area",
                             groupvars = c("Type", "Oligo.Length"))

data2B_summary$Oligo.Length<- as.factor(data2B_summary$Oligo.Length)

# reorder dataframe
data2B_summary <- (data2B_summary %>%
  mutate(Type = factor(Type, levels = c("Modified", "Unmodified", "Standard"))) %>%
           arrange(Type))
pd <- position_dodge( width = 0.8)

make_mod_plot_2B <- ggplot(data = data2B_summary,
                        aes(fill = Type,
                            x = Oligo.Length,
                            y = Peak.Area)) +
  geom_bar(position = pd, stat= "identity",
           width = 0.6,
           size = 1.5,
           aes(color = Type)) +
  geom_errorbar(aes(ymin = Peak.Area-se, ymax = Peak.Area+se), 
                width = 0.4, position = pd) +
  scale_color_manual(values = color_bars) +
  scale_fill_manual(values = fill_alpha(color_bars, 0.4)) +
  xlab("RNA acceptor length (nucleotides)") +
  scale_y_continuous(expand = c(0,0), labels = function(x) x / 1e6,
                     limits = c(0, 15000000))+
  ylab(expression("Peak Area *"*10^6*"")) +
theme(plot.margin = margin(12,12,12,12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(color = "black",
                                 size = axis_text_size),
      axis.text.y = element_text(color = "black",
                                 size = axis_text_size),
      axis.title.x = element_text(size = axis_title_size),
      axis.title.y = element_text(size = axis_title_size),
      legend.text = element_text(color = "black",
                                 size = axis_text_size),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.5, 0.85),
      strip.background = element_blank(),
      strip.text = element_blank())

make_mod_plot_2B
ggsave("plots/Figure2B.pdf", dpi = 300,
       width = 4, height = 3, unit = "in")

# 2C example chromatograms pre and post ligation
chromatogramData_2C <- read.csv("Figure2/Figure2C_chromatogram.csv",
                                header = TRUE)

chromatogramData_2C <- subset(chromatogramData_2C, Y.Counts != 0)

chromatogramData_2C_std <- subset(chromatogramData_2C, Type == "Standard" 
                              & X.Minutes < 18)
chromatogramData_2C_mod <- subset(chromatogramData_2C, Type == "Modified" 
                                  & X.Minutes > 18)
chromatogramData_2C <- rbind(chromatogramData_2C_std, chromatogramData_2C_mod)

chromatogramData_2C$Oligo.Length <- as.factor(chromatogramData_2C$Oligo.Length)
chromatogramData_2C$Type <- as.factor(chromatogramData_2C$Type)

chromatogramData_2C <- chromatogramData_2C %>%
  arrange(X.Minutes, Type, Oligo.Length)




# change shading inside geom_ribbon

# options(bitmapType = "cairo")

makeChromatogram_2C <- ggplot(data = chromatogramData_2C,
                              aes(x = X.Minutes,
                                  y = Y.Counts / 10000)) +
  facet_wrap2(~ Type,
              axes = "x",
              remove_labels = "x",
              ncol = 1) +
  geom_line(size = 1, alpha = 0.9,
            aes(color = interaction(factor(Oligo.Length), 
                                    factor(Type)),
                linetype = factor(Oligo.Length))) +
  geom_ribbon(show.legend = FALSE,
              alpha = 0.6,
              aes(
                  ymin = 0,
                  ymax = Y.Counts / 10000,
                  fill = interaction(factor(Oligo.Length), 
                                     factor(Type))
                )) +
  scale_color_manual(guide = "none",
                     values = chrom_colors) +
  scale_fill_manual(guide = "none",
                    values = chrom_colors) +
  xlab("Time (minutes)") +
  ylab(expression("Intensity *"*10^6*"")) +
  scale_linetype_manual(name = "Oligo Length", 
                        values = lty_pal,
                        guide = guide_legend(
                          override.aes = list(color = "black",
                                              fill = "black"))) +
  scale_x_continuous(breaks = seq(13,25, by=1), limits=c(13,25),
                      expand = c(0,0.2)) +
  scale_x_break(c(17, 21))+
  scale_y_continuous(expand = c(0,0.05)) +
  # labs(linetype = "Oligo Length") +
  theme(plot.margin = margin(8,150,8,8),
        panel.spacing.y=unit(0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(angle = 0,
                                   color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.length.x.top = unit(0, "pt"),
        axis.title.y = element_text(angle = 90,
                                    size = axis_title_size),
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = axis_title_size),
        legend.background = element_blank(),
        legend.position = "none",
        # legend.position.inside = c(0.2, 0.2),
        strip.background = element_blank(),
        strip.text = element_blank())
          

makeChromatogram_2C

ggsave("plots/Figure2C.pdf", dpi = 300,
       width = 3, height = 2.25, unit = "in")

# panel D spectra
spectra_2D <- read.csv("Figure2/Figure2D_ligation_spectra.CSV")
lig_prod <- subset(spectra_2D, Identity == "Ligation Product")
rna_acceptor <- subset(spectra_2D, Identity == "RNA")

make_spectrum <- function(exp_mass, conj_data, col, num_breaks) {
  # determine bounds
  x_min <- as.numeric(exp_mass) - 1
  x_max <- as.numeric(exp_mass) + 2
  # x_min <- as.numeric(728)
  # x_max <- as.numeric(756)
  # subset data so it runs faster
  data <- subset(conj_data, X.Thomsons. > as.numeric(x_min) &
                   X.Thomsons. < as.numeric(x_max))
  make_spectrum <- ggplot(data = conj_data,
                          aes(x = X.Thomsons.,
                              y = Y.Counts./10000)) +
    geom_point(size = 0.1, 
               color = col,
               fill = col) +
    geom_line(size = 0.2, 
              color = col) +
    scale_x_continuous(limits = c(as.numeric(x_min), as.numeric(x_max)),
                       breaks = seq(floor(x_min), ceiling(x_max), by = 1),
                       minor_breaks = seq(floor(x_min), ceiling(x_max), 
                                          by = 1.0 / num_breaks),
                       guide = guide_axis(minor.ticks = TRUE)) +
    scale_y_continuous(limits = c(0, 4.2)) +
    xlab("m/z") +
    ylab(expression("Intensity *"*10^4*"")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.minor.ticks.length = rel(0.5),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black",
                                     angle = 45,
                                     hjust = 1,
                                     vjust = 1,
                                     size = axis_text_size),
          axis.text.y = element_text(color = "black",
                                     size = axis_text_size),
          # axis.title.x = element_blank(),
          axis.title.x = element_text(face = "bold.italic",
                                      size = axis_title_size),
          # axis.title.y = element_blank(),
          axis.title.y = element_text(size = axis_title_size),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.title = element_text(size = axis_title_size),
          legend.text = element_text(color = "black",
                                     size = axis_text_size),
          plot.margin=grid::unit(c(2,2,2,2), "mm"))
  # return ggplot object
  return (make_spectrum)
}

make_spectrum(1488.245, lig_prod, color_set[1], 4)
ggsave("plots/Figure2D_right.pdf", dpi = 300,
       width = 1.65, height = 2, unit = "in")
make_spectrum(1354.847, rna_acceptor, "#000000", 3)
ggsave("plots/Figure2D_left.pdf", dpi = 300,
       width = 1.65, height = 2, unit = "in")

# Figure 3: T1 Workflow optimization ----
areas_31mer <- read.csv("Figure3_areas.csv")
# summarize for figure generation
areas_31mer_summary <- summarySE(areas_31mer, 
                            measurevar = "peak.area",
                            groupvars = c("oligo", "donor"))
make_barplot_3B <- ggplot(data = areas_31mer_summary,
                          aes(fill = donor,
                              x = oligo,
                              y = peak.area)) +
  geom_bar(position = pd, stat= "identity",
           width = 0.6,
           size = 1,
           aes(color = donor)) +
  geom_errorbar(aes(ymin = peak.area-se, ymax = peak.area+se), 
                width = 0.4, position = pd) +
  scale_color_manual(values = color_bars[c(1,3)]) +
  scale_fill_manual(values = fill_alpha(
                    (color_bars[c(1,3)]), 0.4)) +
  xlab("Digest Fragment Sequence") +
  ylab(expression("Peak Area *"*10^6*"")) +
  scale_y_continuous(expand = c(0,0),
                     labels = function(x) x / 1e6) +
  theme(plot.margin = margin(12,12,12,12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size,
                                   angle = 20,
                                   hjust = 1),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size,),
        axis.title.y = element_text(size = axis_title_size),
        legend.text = element_text(color = "black",
                                   size = axis_text_size),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.6),
        strip.background = element_blank(),
        strip.text = element_blank())

make_barplot_3B
title <- "Figure3B"
ggsave(paste0("plots/", title, ".pdf"), dpi = 300,
       width = 3.25, height = 2.5, unit = "in")

# Figure 4: Application to tRNA digest ----

# 4B peak area changes

tRNA_peakarea_4B <- read.csv("Figure4/tRNA_peakareas.csv")
tRNA_peakarea_4B$Identity <- factor(
  tRNA_peakarea_4B$Identity,
  levels = c("RNA", "Ligation Product"))
tRNA_peakarea_4B$Fragment <- factor(
  tRNA_peakarea_4B$Fragment,
  levels = unique(tRNA_peakarea_4B$Fragment[order(-tRNA_peakarea_4B$Order)])
)

tRNA_foldchange <- pivot_wider(tRNA_peakarea_4B[,-4], 
                               names_from = Identity, 
                               values_from = Peak.Area)
tRNA_foldchange$foldchange <- tRNA_foldchange$`Ligation Product` / tRNA_foldchange$RNA

make_barplot_4B <- ggplot(data = tRNA_peakarea_4B,
                          aes(fill = Identity,
                              x = Fragment,
                              y = Peak.Area)) +
  geom_bar(position = pd, stat= "identity",
           width = 0.6,
           size = 1,
           aes(color = Identity)) +
  scale_color_manual(values = rev(color_bars[c(1,3)])) +
  scale_fill_manual(values = fill_alpha
                    (rev(color_bars[c(1,3)]), 0.4)) +
  xlab("tRNA fragment sequence") +
  ylab(expression("Peak Area *"*10^5*"")) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0),
                     position = "right",
                     labels = function(x) x / 1e5) +
  theme(plot.margin = margin(12,12,12,12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(color = "black",
                                   angle = -60,
                                   vjust = 1,
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.text = element_text(color = "black",
                                   size = axis_text_size),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        strip.background = element_blank(),
        strip.text = element_blank())

make_barplot_4B
title <- "Figure4B"
ggsave(paste0("plots/", title, ".pdf"), dpi = 300,
       width = 3, height = 3.8, unit = "in")

# 4C example chromatograms for two of the oligonucleotides
tRNA_data <- read.csv("Figure4/Figure4C_chromatogram_SM.CSV")
tRNA_lig_data <-  read.csv("Figure4/Figure4C_chromatogram_lig_prod.CSV")

tRNA_data_4C <- rbind(tRNA_data, tRNA_lig_data)
tRNA_data_4C <- tRNA_data_4C[order(tRNA_data_4C$Sequence), ]
tRNA_data_4C <- subset(tRNA_data_4C, Sequence == "CUCAG" |
                      Sequence == "m5U_CG")
# To make chromatogram in graphical abstract use 
# tRNA_data_4C <- subset(tRNA_data_4C,
#                       Sequence == "m5U_CG")

makeChromatogram_4C <- ggplot(data = tRNA_data_4C,
                              aes(x = X.Minutes.,
                                  y = Y.EIC. / 1000)) +
  facet_wrap2(~ Identity,
              axes = "x",
              remove_labels = "x",
              ncol = 1) +
  geom_line(size = 0.5, alpha = 0.9,
            aes(color = interaction(factor(Sequence), 
                                    factor(Identity))))+
                # linetype = factor(Sequence))) +
  geom_ribbon(show.legend = FALSE,
              alpha = 0.6,
              aes(
                ymin = 0,
                ymax = Y.EIC. / 1000,
                fill = interaction(factor(Sequence), 
                                        factor(Identity))
              )) +
  scale_color_manual(guide = "none",
                     values = chrom_colors_long[c(1,3,6,8)]) +
  scale_fill_manual(guide = "none",
                    values = chrom_colors_long[c(1,3,6,8)]) +
  xlab("Time (minutes)") +
  ylab(expression("Intensity *"*10^4*"")) +
  scale_linetype_manual(name = "Oligo Length", 
                        values = lty_pal,
                        guide = "none") +
  scale_x_continuous(breaks = seq(5,30, by=5), limits=c(4,31),
                     expand = c(0,0.2)) +
  # scale_x_break(c(17, 21))+
  scale_y_continuous(expand = c(0,0.05)) +
  # labs(linetype = "Oligo Length") +
  theme(plot.margin = margin(8,8,8,8),
        panel.spacing.y=unit(0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(angle = 0,
                                   color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.length.x.top = unit(0, "pt"),
        axis.title.y = element_text(angle = 90,
                                    size = axis_title_size),
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = axis_title_size),
        legend.background = element_blank(),
        legend.position = "none",
        # legend.position.inside = c(0.8, 0.2),
        strip.background = element_blank(),
        strip.text = element_blank())

makeChromatogram_4C

ggsave("plots/Figure4C.pdf", dpi = 300,
       width = 4, height = 2.25, unit = "in")

# panel D example spectra for TPCG before and after ligation

tRNAlig_spec_4D <- read.csv("Figure4/Figure4D_TPCG_TAG_spectrum.CSV")
tRNA_spec_4D <- read.csv("Figure4/Figure4D_TPCG_spectrum.CSV")

make_tRNA_spectrum <- function(exp_mass, conj_data, col) {
  # determine bounds
  x_min <- as.numeric(exp_mass) - 1
  x_max <- as.numeric(exp_mass) + 2
  # x_min <- as.numeric(728)
  # x_max <- as.numeric(756)
  # subset data so it runs faster
  data <- subset(conj_data, X.Thomsons. > as.numeric(x_min) &
                   X.Thomsons. < as.numeric(x_max))
  make_spectrum <- ggplot(data = conj_data,
                               aes(x = X.Thomsons.,
                                   y = Y.Counts./1000)) +
    geom_point(size = 0.1, 
               color = col,
               fill = col) +
    geom_line(size = 0.2, 
              color = col) +
    scale_x_continuous(limits = c(as.numeric(x_min), as.numeric(x_max))) +
    scale_y_continuous(limits = c(0, 4.2)) +
    xlab("m/z") +
    ylab(expression("Intensity *"*10^3*"")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.minor.ticks.length = rel(0.5),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black",
                                     angle = 20,
                                     hjust = 1,
                                     vjust = 1,
                                     size = axis_text_size),
          axis.text.y = element_text(color = "black",
                                     size = axis_text_size),
          # axis.title.x = element_blank(),
          axis.title.x = element_text(face = "bold.italic",
                                        size = axis_title_size),
          # axis.title.y = element_blank(),
          axis.title.y = element_text(size = axis_title_size),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.title = element_text(size = axis_title_size),
          legend.text = element_text(color = "black",
                                     size = axis_text_size),
          plot.margin=grid::unit(c(2,2,2,2), "mm"))
  # return ggplot object
  return (make_spectrum)
}

make_tRNA_spectrum(1033.546,tRNAlig_spec_4D, col = color_set[1])
title <- "TPCGcctag_decyl_spectrum"
ggsave(paste0("plots/", title, ".pdf"), dpi = 300,
       width = 2, height = 1.8, unit = "in")
make_tRNA_spectrum(606.092,tRNA_spec_4D, col = "#000000")
title <- "TPCG_spectrum"
ggsave(paste0("plots/", title, ".pdf"), dpi = 300,
       width = 2, height = 1.8, unit = "in")





# Supplemental Figures -----

# S Fig 1: High res mass spectrum of oligo derivatives -----
# makes spectrum, centers x axis on exp_mass
make_spectrum <- function(exp_mass, conj_data) {
  # determine bounds
  x_min <- as.numeric(exp_mass) - 1
  x_max <- as.numeric(exp_mass) + 2
  # x_min <- as.numeric(728)
  # x_max <- as.numeric(756)
  # subset data so it runs faster
  data <- subset(conj_data, X.Thomsons. > as.numeric(x_min) &
                   X.Thomsons. < as.numeric(x_max))
  make_conj_spectrum <- ggplot(data = conj_data,
                               aes(x = X.Thomsons.,
                                   y = Y.Counts./10000)) +
    geom_point(size = 0.1) +
    geom_line(linewidth = 0.2) +
    scale_x_continuous(limits = c(as.numeric(x_min), as.numeric(x_max))) +
    xlab("m/z") +
    ylab(expression("Intensity *"*10^4*"")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.minor.ticks.length = rel(0.5),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black",
                                     size = axis_text_size),
          axis.text.y = element_text(color = "black",
                                     size = axis_text_size),
          axis.title.x = element_blank(),
          # axis.title.x = element_text(face = "bold.italic",
          #                               size = axis_title_size),
          axis.title.y = element_blank(),
          # axis.title.y = element_text(size = axis_title_size),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.title = element_text(size = axis_title_size),
          legend.text = element_text(color = "black",
                                     size = axis_text_size),
          plot.margin=grid::unit(c(0,0,0,0), "mm"))
  # return ggplot object
  return (make_conj_spectrum)
}

files <- list.files(path="supplemental/conjugate_spectra/", full.names=TRUE, recursive=FALSE)

for (file in files) {
  conj_data <- read.csv(file, skip = 1)
  exp_mass <- conj_data$X.Thomsons.[which.max(conj_data$Y.Counts.)]
  # print(exp_mass)
  plot <- make_spectrum(exp_mass, conj_data)
  plot
  title <- tools::file_path_sans_ext(basename(file))
  title <- paste0(title,"_ultrazoom")
  ggsave(paste0("plots/Supplement/", title, ".pdf"), dpi = 300,
         width = 1.4, height = 0.8, unit = "in")
}
# S Fig 2: Ion Tagged Oligonucleotide IEs -----
# summarize for figure generation
ITO_data <- read.csv("supplemental/FigureS2_IonTag_IEs.csv", header = TRUE)

dataS2_summary <- summarySE(ITO_data, 
                            measurevar = "Area",
                            groupvars = c("Concentration..uM.", "Oligo", "Order"))
dataS2_summary <- dataS2_summary %>%
  mutate(Oligo = factor(Oligo, 
                        levels = c("decyl", "octyl", "methyl", "unmodified"))) 
S2_colors <- color_set[c(1,3,4,5)]

makeIEplot_S2 <- ggplot(data = dataS2_summary,
                        aes(x = Concentration..uM.,
                            y = Area,
                            color = Oligo,
                            fill = Oligo)) +
  geom_errorbar(aes(ymin = Area-se, ymax = Area+se), 
                width = 0.5) +
  geom_point(size = 0.5) +
  geom_smooth(method='lm', formula=y~x, se=FALSE, linewidth = 0.5) +
  scale_fill_manual(values = S2_colors)+
  scale_color_manual(values = S2_colors)+
  scale_y_continuous(expand = c(0,0), labels = function(x) x / 1e6,
                     limits = c(0, 3500000))+
  scale_x_continuous(guide = "prism_minor",
                     breaks = seq (0, 14, 4),
                     minor_breaks = seq (0, 14, 1),
                     expand = c(0,0),
                     limits = c(0,14))+
  guides(y = "prism_minor")+
  xlab("Concentration (µM)") +
  ylab(expression("Peak Area *"*10^5*"")) +
  labs(color = "3′ chemistry",
       fill = "3′ chemistry") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.minor.ticks.length = rel(0.5),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.7),
  )

makeIEplot_S2

ggsave("plots/FigureS2.pdf",
       width = 3.5, height = 3, unit = "in")

# S Fig 3: T4 RNL1 Ligation optimization ----- 
SM_data <- read.csv("supplemental/FigureS3_kinetics_SM.csv")
names(SM_data) <- sub('X', '', names(SM_data))
SM_data <- pivot_longer(SM_data, !Time, 
                        names_to = "T4.RNL1", 
                        values_to = "Peak.Area")
SM_data <- na.omit(SM_data)
Prod_data <- read.csv("supplemental/FigureS3_kinetics_prod.csv")
names(Prod_data) <- sub('X', '', names(Prod_data))
Prod_data <- pivot_longer(Prod_data, !Time, 
                          names_to = "T4.RNL1", 
                          values_to = "Peak.Area")
Prod_data <- na.omit(Prod_data)

make_SM_opt <- ggplot(data = SM_data,
                      aes(x = Time,
                          y = Peak.Area / 1000000, 
                          color = T4.RNL1)) +
  geom_point(size = 1) +
  geom_line(size = 1, alpha = 0.5) +
  xlab("Time (minutes)") +
  ylab(expression("EICPeak Volume *"*10^6*"")) +
  labs(color = "Concentration of \n T4 RNA Ligase \n (U / µL)") +
  scale_color_manual(values = colors) +
  theme(strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.minor.ticks.length = rel(0.5),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.8, 0.6),
        legend.title = element_text(size = axis_title_size),
        legend.text = element_text(color = "black",
                                   size = axis_text_size))

make_SM_opt
ggsave("plots/FigureS3_SM.pdf", dpi = 300,
       width = 3.2, height = 2.5, unit = "in")

make_Prod_opt <- ggplot(data = Prod_data,
                        aes(x = Time,
                            y = Peak.Area / 1000000, 
                            color = T4.RNL1)) +
  geom_point(size = 1) +
  geom_line(size = 1, alpha = 0.5) +
  xlab("Time (minutes)") +
  ylab(expression("EICPeak Volume *"*10^6*"")) +
  scale_color_manual(values = colors) +
  theme(strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.minor.ticks.length = rel(0.5),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = axis_title_size),
        legend.text = element_text(color = "black",
                                   size = axis_text_size))

make_Prod_opt
ggsave("plots/FigureS3_Prod.pdf", dpi = 300,
       width = 3.2, height = 2.5, unit = "in")
# S Fig 6: RNA Acceptors remaining post ligation ----

chromatogramData_S6 <- read.csv("supplemental/FigureS6_acceptor_EICs.csv",
                                skip = 1,
                                header = TRUE)

chromatogramData_S6 <- chromatogramData_S6 %>% arrange(Oligo)

makeChromatogram_S6 <- ggplot(data = chromatogramData_S6,
                              aes(x = X.Minutes.,
                                  y = Y.Counts. / 10000)) +
  facet_wrap2(~ fct_rev(Oligo),
              axes = "x",
              remove_labels = "x",
              ncol = 1) +
  geom_line(aes(color = fct_rev(Oligo))) +
  geom_ribbon(alpha = 0.6, aes(ymin = 0, ymax = Y.Counts. / 1000, 
                               fill = fct_rev(Oligo))) + 
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  xlab("Time (minutes)") +
  ylab(expression("Intensity *"*10^3*"")) +
  # scale_x_continuous(breaks = seq(9,15, by=1), limits=c(9,15),
  #                    expand = c(0,0.2)) +
  scale_y_continuous(expand = c(0,0.05)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y=unit(0.5, "lines"),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.text = element_text(color = "black",
                                   size = axis_text_size),
        # legend.key = element_blank(),
        legend.key.size = unit(0.1, 'in'),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.9),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

makeChromatogram_S6

ggsave("plots/FigureS6.pdf", dpi = 300,
       width = 6.5, height = 4, unit = "in")
# S Fig 7: Biotinylation of Ribonuclease T1----
BiotinT1_data <- read.csv("supplemental/FigureS7_T1_spectra.CSV")
BiotinT1_data <- subset(BiotinT1_data, X.Thomsons. > 2500 &
                          X.Thomsons. < 3000)

max_stand <- max(subset(BiotinT1_data, Sample == "standard")$Y.Counts.)
max_30 <- max(subset(BiotinT1_data, Sample == "30 min")$Y.Counts.)
max_90 <- max(subset(BiotinT1_data, Sample == "90 min")$Y.Counts.)

BiotinT1_data <- BiotinT1_data |>
  mutate(Y.Counts.Normalized = case_when(
    Sample == "standard" ~ Y.Counts. / max_stand,
    Sample == "30 min"   ~ Y.Counts. / max_30,
    Sample == "90 min"   ~ Y.Counts. / max_90,
    .default = NA_real_
  ))

BiotinT1_data$Sample = factor(BiotinT1_data$Sample, 
                              levels=c('90 min','30 min','standard'))

x_min <- 2750.0
x_max <- 3000.0
make_spectrum <- ggplot(data = BiotinT1_data,
                        aes(x = X.Thomsons.,
                            y = Y.Counts.Normalized, 
                            color = Sample)) +
  facet_wrap2(~ fct_rev(Sample),
              axes = "x",
              remove_labels = "x",
              ncol = 1) +
  # geom_point(size = 0.1) +
  geom_line(size = 0.2) +
  scale_x_continuous(limits = c(as.numeric(x_min), as.numeric(x_max))) +
  xlab("m/z") +
  ylab("Relative Intensity") +
  scale_color_manual(values = colors) +
  theme(strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.minor.ticks.length = rel(0.5),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(face = "bold.italic",
                                    size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = axis_title_size),
        legend.text = element_text(color = "black",
                                   size = axis_text_size))
make_spectrum

ggsave("plots/FigureS7A.pdf", dpi = 300,
       width = 6.5, height = 4, unit = "in")

# load in pytheas final reports and clean data 
extract_digest_time <- function(str) {
  hour <- grepl("hr", str)
  # make num
  digest_time_str <- strsplit(str, "-")[[1]][2]
  digest_time_num <- as.numeric(gsub("([0-9]+).*$", "\\1", digest_time_str))
  # return NA for controls
  if (is.na(digest_time_num)) {
    return(NULL)
  }
  else if (hour) {
    return(as.numeric(60*digest_time_num))
  }
  # convert to mins
  return(digest_time_num)
}
combined_data_list <- list()
for (file in list.files(path="supplemental/pytheas_reports/", 
                        pattern="*.csv", 
                        full.names=TRUE, 
                        recursive=FALSE)) {
  print(file)
  pytheas_data <- data.frame(read.csv(file))
  # extract sample name from data file string
  sample_name <- strsplit(strsplit(file, "/")[[1]][4], "_")[[1]][1]
  pytheas_data$sample <- sample_name
  pytheas_data$dig_time <- extract_digest_time(sample_name)
  combined_data_list <- append(combined_data_list, list(pytheas_data)) 
}

pytheas_df <- bind_rows(combined_data_list)
trimws(pytheas_df$sequence)
pytheas_df$sequence <- as.character(pytheas_df$sequence)
pytheas_df <- subset(pytheas_df,
                     length >= 4 &
                       Rank == 1.0 &
                       molecule.ID != "decoy")

n_distinct(pytheas_df$sample)

## Filter step to remove samples we aren't interested in
pytheas_df <- subset(pytheas_df, 
                     sample != "positive" &
                       sample != "23S-15min" &
                       sample != "post-standard" &
                       sample != "pre-standard")

# arrange in order of score
pytheas_df <- pytheas_df %>% arrange(desc(Score..Sp.))
# remove duplicate sequences in the same sample (highest score match is kept)
pytheas_df_unique <- pytheas_df %>%
  distinct(sequence, sample, .keep_all = TRUE)

# Subsetting data to find unique counts 
unique_counts <- pytheas_df_unique %>%
  group_by(pytheas_df_unique$sequence) 

# using all hits
length(unique(pytheas_df_unique$sequence))
hits_seq_sample <- pytheas_df_unique %>%
  group_by(sequence, sample) %>%
  summarise(average_length = mean(length), 
            dig_time = dig_time,
            score = mean(Score..Sp.), n = n())
# discard any oligo matches that appear in multiple samples
unique_hits_seq <- hits_seq_sample %>%
  group_by(sequence) %>%
  filter(n() == 1) %>%
  ungroup()

shared_hits_seq <- hits_seq_sample %>%
  group_by(sequence) %>%
  filter(n_distinct(sample) == 4) %>%
  summarise(
    average_length = mean(average_length),
    dig_time = first(dig_time),   # or mean(dig_time) if it varies
    score = mean(score),
    n = sum(n),
    sample = "all"
  ) %>%
  ungroup()

hits <- rbind(unique_hits_seq, shared_hits_seq)
hits <- subset(hits, score > 0.004)

unique_hits_seq$missed_cleavages <- (str_count(unique_hits_seq$sequence, "G")) -1

# Filter for sequences found in only one sample
unique_sequences <- unique_counts %>%
  filter(sample_count == 1)

# Subset original dataframe
df_unique <- df %>%
  filter(sequence %in% unique_sequences$sequence)

# visualization
pytheas_df <- subset(pytheas_df,
                     Score..Sp. > 0.004 &
                       molecule.ID != "decoy" &
                       length >= 4)

oligos_2mins <- as.list(subset(pytheas_df, 
                                 dig_time == 2)$sequence)
oligos_2mins <- oligos_2mins[!duplicated(oligos_2mins)]

oligos_5mins <- as.list(subset(pytheas_df, 
                                 dig_time == 5)$sequence)
oligos_5mins <- oligos_5mins[!duplicated(oligos_5mins)]

oligos_60mins <- as.list(subset(pytheas_df, 
                                  dig_time == 60)$sequence)
oligos_60mins <- oligos_60mins[!duplicated(oligos_60mins)]

oligos_120mins <- as.list(subset(pytheas_df, 
                                  dig_time == 120)$sequence)
oligos_120mins <- oligos_120mins[!duplicated(oligos_120mins)]

oligo_list <- list(oligos_2mins,
                    oligos_5mins,
                    oligos_60mins,
                   oligos_120mins)

ggVennDiagram(oligo_list, label_alpha = 0,
              stroke_size = 0.5, set_name_size = 0.5,
              category.names = c("2 min",
                                 "5 min",
                                 "1 hr",
                                 "2 hr")) +
  scale_fill_gradientn(colors = c("white", color_set[3], chrom_colors[3], color_set[6]),
                       values = c(0, 0.15, 0.2, 1))

ggsave("plots/pytheas-venn-diagram.pdf", dpi = 300,
       width = 4, height = 3 , unit = "in")

# summarize length data by digestion time
oligo_lengths <- summarySE(data = subset(pytheas_df_unique, 
                                         molecule.ID != "decoy" &
                                           Rank == 1.0 & Score..Sp. > 0.05 ), 
                           measurevar = "length",
                           groupvars = "dig_time")

#ANOVA
hits$sample <- as.factor(hits$sample)
summary(hits)
res_aov <- aov(average_length ~ sample,
               data = hits)
summary(res_aov)

# visualize plots to assess normality of data
par(mfrow = c(1, 2)) # combine plots
# histogram
hist(res_aov$residuals)
# QQ-plot
library(car)
qqPlot(res_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# statistical test for normality
shapiro.test(res_aov$residuals)
# this test produces a very small p value suggesting that
# our data are NOT normally distributed, but since we have >30 data points
# from each group, we can proceed given the central limit theorem (?)

# statistical test
leveneTest(average_length ~ as.factor(sample),
           data = hits)

# Post hoc tests
library(multcomp)
# Tukey HSD test:
post_test <- glht(res_aov,
                  linfct = mcp(sample = "Tukey")
)

summary(post_test)

# instead using a ranked statistical test (non-parametric)
# there is no assumption of normality
kruskal.test(average_length ~ sample, 
             data = hits)

pairwise.wilcox.test(hits$average_length,
                     as.factor(hits$sample),
                     p.adjust.method = "BH")

hits <- hits %>%
  mutate(sample = factor(sample,
                           levels = c("23S-2min", 
                                      "23S-5min", 
                                      "23S-1hr",
                                      "23S-2hr",
                                      "all"),
                           labels = c("2 min only", 
                                      "5 min only", 
                                      "1 hr only", 
                                      "2 hr only",
                                      "all samples")
  ))

ggplot(subset(hits, sample != "all samples"),
       aes(x = as.factor(sample), 
           y = as.numeric(average_length))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(
    position = position_jitter(height = 0, width = 0.2),
    shape = 1,
    alpha = 0.6
  ) +
  scale_y_continuous(
    breaks = seq(5, 16, by = 2),
    minor_breaks = seq(5, 16, by = 1)
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 11),
        axis.text.y = element_text(color = "black",
                                   size = 11),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 11),
        legend.text = element_text(color = "black",
                                   size = 11),
        legend.title = element_text(color = "black",
                                    size = 11),
        legend.background = element_blank()) +
  # legend.position = "NONE",
  labs(x = "Sample", y = "Length", )

ggsave("plots/pytheas-box-wisker.pdf", dpi = 300,
       width = 6, height = 2.5, unit = "in")


# S Fig 8: T4 PNK Optimization -----

data  <- read.csv("supplemental/FigureS8_T4PNK_opt.csv")

colors <- rev(color_set[2:5])

custom_order <- c("3′ cP", "3′ P", "3′ OH")
data$X3..Chemistry <- factor(data$X3..Chemistry, levels = custom_order)

make_T4PNK_plot_S3 <- ggplot(data = data,
                           aes(fill = X3..Chemistry,
                               x = Oligonucleotide,
                               y = Peak.Area)) +
  facet_wrap2(facets = vars(T4.PNK),
              ncol = 1,
              axes = "x") +
  geom_bar(position = pd, stat= "identity",
           width = 0.6,
           size = 1.5,
           aes(color = X3..Chemistry)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = fill_alpha(colors, 0.4)) +
  xlab("T1 digestion fragment") +
  scale_y_continuous(expand = c(0,0), labels = function(x) x / 1e5)+
  ylab(expression("Extracted Ion Chromatogram Peak Area *"*10^5*"")) +
  # labs(color = "3′ Chemistry") +
  theme(plot.margin = margin(12,12,12,12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.title = element_blank(),
        legend.text = element_text(color = "black",
                                   size = axis_text_size),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.9),
        strip.background = element_blank(),
        strip.text = element_blank())

make_T4PNK_plot_S8
ggsave("plots/FigureS8.pdf", dpi = 300,
       width = 6.5, height = 7, unit = "in")

# S Fig 10: 200 ng tRNA digest

tRNA_peakarea_S10 <- read.csv("Figure4/tRNA_peakareas_200ng.csv")
tRNA_peakarea_S10$Identity <- factor(
  tRNA_peakarea_S10$Identity,
  levels = c("RNA", "Ligation Product"))

tRNA_peakarea_S10$amount <- 200.0
tRNA_peakarea_4B$amount <- 500.0
all_tRNA_peak_areas <- rbind(tRNA_peakarea_4B, tRNA_peakarea_S10)

all_tRNA_peak_areas$Fragment <- factor(
  all_tRNA_peak_areas$Fragment,
  levels = unique(all_tRNA_peak_areas$Fragment[order(all_tRNA_peak_areas$Order)])
)

make_barplot_S10 <- ggplot(data = all_tRNA_peak_areas,
                           aes(fill = interaction(Identity:as.factor(amount)),
                               x = Fragment,
                               y = Peak.Area)) +
  geom_bar(position = pd, stat= "identity",
           width = 0.6,
           size = 1,
           aes(color = interaction(Identity:as.factor(amount)))) +
  facet_wrap2(~ fct_rev(as.factor(amount)),
              axes = "x",
              nrow = 1) +
  scale_color_manual(values = rev(chrom_colors[c(1,2,4,6)])) +
  scale_fill_manual(values = fill_alpha
                    (rev(chrom_colors[c(1,2,4,6)]), 0.4)) +
  xlab("tRNA fragment sequence") +
  ylab(expression("Peak Area *"*10^5*"")) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 500000),
                     labels = function(x) x / 1e5) +
  theme(plot.margin = margin(12,12,12,12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.y = element_text(color = "black",
                                   size = axis_text_size),
        axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust =01,
                                   size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        legend.text = element_text(color = "black",
                                   size = axis_text_size),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8),
        strip.background = element_blank(),
        strip.text = element_blank())


make_barplot_S10
title <- "FigureS10"
ggsave(paste0("plots/Supplement", title, ".pdf"), dpi = 300,
       width = 6.5, height = 4.5, unit = "in")
