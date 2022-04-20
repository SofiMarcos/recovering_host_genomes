# ALLELE FREQUENCIES #

library(tidyverse)
library(janitor)
library(ggpattern)

# Load Cobb data  -----------------
internal_c <- read.table("3_population_genetics_inference/data/freq/internal-cobb-freq.txt",
                         header = TRUE, sep = "") %>%
  clean_names()
internal_c$method <- "internal"

external_c <- read.table("3_population_genetics_inference/data/freq/external-cobb-freq.txt",
                         header = TRUE, sep = "") %>%
  clean_names()
external_c$method <- "external"

combined_c <- read.table("3_population_genetics_inference/data/freq/combined-cobb-freq.txt",
                         header = TRUE, sep = "") %>%
  clean_names()
combined_c$method <- "combined"

diverse_c <- read.table("diverse-cobb-freq.txt", header = TRUE, sep = "") %>%
  clean_names()
diverse_c$method <- "diverse"

cobb <- rbind(internal_c, external_c, combined_c, diverse_c)
rm(internal_c, external_c, combined_c, diverse_c)



# Load Ross data  -----------------
setwd("3_population_genetics_inference/data/freq")

internal_c <- read.table("3_population_genetics_inference/data/freq/internal-ross-freq.txt",
                         header = TRUE, sep = "") %>%
  clean_names()
internal_c$method <- "internal"

external_c <- read.table("3_population_genetics_inference/data/freq/external-ross-freq.txt",
                         header = TRUE, sep = "") %>%
  clean_names()
external_c$method <- "external"

combined_c <- read.table("3_population_genetics_inference/data/freq/combined-ross-freq.txt",
                         header = TRUE, sep = "") %>%
  clean_names()
combined_c$method <- "combined"

diverse_c <- read.table("3_population_genetics_inference/data/freq/diverse-ross-freq.txt",
                        header = TRUE, sep = "") %>%
  clean_names()
diverse_c$method <- "diverse"

ross <- rbind(internal_c, external_c, combined_c, diverse_c)
rm(internal_c, external_c, combined_c, diverse_c)


# Allele frequencies -----------------
#
# cobb
cobb %>%
  mutate(method = factor(
    method, levels = c("internal","external","combined","diverse"))
  ) %>%
  ggplot(aes(x = maf, group = method, color = method)) + ylim(0,3) +
  geom_density(adjust = 1.5, alpha = 0.4, size = 1) +
  geom_vline(xintercept = 0.05, size = 1) +
  # text
  ylab("maf") + xlab("Density") + ggtitle(label = "Cobb") +
  # color
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#9281e2")) +
  theme_bw() +
  theme(
    # panel
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend
    legend.position = "none",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    # text
    axis.title.x = element_text(size = 9,face = "bold"),
    axis.title.y = element_text(size = 9,face = "bold"),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold")
  )

# ross
ross %>%
  mutate(method = factor(
    method, levels = c("internal","external","combined","diverse")
  )) %>%
  ggplot(aes(x = maf, group = method, color = method)) +
  geom_density(adjust = 1.5, alpha = .4, size = 1) +
  geom_vline(xintercept = 0.05, size = 1) +
  # text
  ylab("maf") + xlab("Density") + ylim(0,3) + ggtitle(label = "Ross") +
  # color
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#9281e2")) +
  theme_bw() +
  theme(
    # panel
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend
    legend.position = "none",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    # text
    axis.title.x = element_text(size = 9,face = "bold"),
    axis.title.y = element_text(size = 9,face = "bold"),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold")
  )


# Allele frequencies higher than 0.5 -----------------
# ross
ross %>%
  mutate(
    maf = if_else(
      condition = maf > 0.5, true = 1 - maf , false = maf)) %>%
  filter(maf>0)

ross_clean %>%
  mutate(method = factor(
    method, levels = c("internal","external","combined","diverse")
  )) %>%
  ggplot(aes(x = maf, color = method)) + ylim(0,3.5) +
  geom_density(adjust = 1.5, alpha = 0.4, size = 1, outline.type = "full") +
  # text
  ylab("Density") + xlab("maf") + ggtitle(label = "Ross") +
  # color
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#9281e2")) +
  theme_bw() +
  theme(
    # panel
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend
    legend.position = "none",
    legend.title = element_text(size = 9, face = "bold"),
    # text
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold")
  )

# cobb
cobb_clean <- cobb %>%
  mutate(
    maf = if_else(
      condition = maf > 0.5, true = 1 - maf , false = maf)) %>%
  filter(maf > 0)

cobb_clean %>%
  mutate(method = factor(
    method, levels = c("internal","external","combined","diverse")
  )) %>%
  ggplot(aes(x = maf, color = method)) + ylim(0,3.5) +
  geom_density(adjust = 1.5, alpha = 0.4, size = 1, outline.type = "full") +
  # text
  ylab("Density") + xlab("maf") + ggtitle(label = "Cobb") +
  # color
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#9281e2")) +
  theme_bw() +
  theme(
    # panel
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend
    legend.position = "none",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    # text
    axis.title.x = element_text(size = 9,face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold")
  )


# maf frequencies  -----------------
mafs <- read.csv("3_population_genetics_inference/data/freq/Cobb-RossMAFs.csv",
                 header = TRUE, sep = ";") %>%
  clean_names()


mafs %>%
  pivot_longer(internal:diverse, names_to = "panel", values_to = "percentage") %>%
  mutate(panel = factor(
    panel, levels = c("internal","external","combined","diverse")
    )) %>%
  ggplot(aes(x = panel, y = percentage, fill = panel, pattern = breed)) +
  # text
  xlab("") + ylab("mafs >0.05 (%)") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  # pattern
  geom_bar_pattern(position = "dodge", stat = "identity", alpha = 0.6,
                   color = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c("none","stripe")) +
  # color
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#9281e2")) +
  theme_bw() +
  theme(
    # panel
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    # text
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9,face = "bold"),
    axis.text = element_text(size = 9),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    plot.title = element_text(size = 11, face = "bold")
  ) +
  # guides
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))
