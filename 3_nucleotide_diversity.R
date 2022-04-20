# NUCLEOTIDE DIVERSITY #

library(tidyverse)
library(ggpubr)
library(rstatix)
library(janitor)

# Load data  -----------------
# internal
cobbi <- read.csv("3_population_genetics_inference/data/nucleo_div/internal-cobb.windowed.pi",
                  header = TRUE, sep = "")
cobbi$panel <- "internal"
cobbi$breed <- "cobb"

rossi <- read.csv("3_population_genetics_inference/data/nucleo_div/internal-ross.windowed.pi",
                  header = TRUE, sep = "")
rossi$panel <- "internal"
rossi$breed <- "ross"

# external
cobbe <- read.csv("3_population_genetics_inference/data/nucleo_div/external-cobb.windowed.pi",
                  header = TRUE, sep = "")
cobbe$panel <- "external"
cobbe$breed <- "cobb"

rosse <- read.csv("3_population_genetics_inference/data/nucleo_div/external-ross.windowed.pi",
                  header = TRUE, sep = "")
rosse$panel <- "external"
rosse$breed <- "ross"

# combined
cobbc <- read.csv("3_population_genetics_inference/data/nucleo_div/combined-cobb.windowed.pi",
                  header = TRUE, sep = "")
cobbc$panel <- "combined"
cobbc$breed <- "cobb"

rossc <- read.csv("3_population_genetics_inference/data/nucleo_div/combined-ross.windowed.pi",
                  header = TRUE, sep = "")
rossc$panel <- "combined"
rossc$breed <- "ross"

# diverse
cobbd <- read.csv("3_population_genetics_inference/data/nucleo_div/diverse-cobb.windowed.pi",
                  header = TRUE, sep = "")
cobbd$panel <- "diverse"
cobbd$breed <- "cobb"

rossd <- read.csv("3_population_genetics_inference/data/nucleo_div/diverse-ross.windowed.pi",
                  header = TRUE, sep = "")
rossd$panel <- "diverse"
rossd$breed <- "ross"

# HD
cobbv <- read.csv("3_population_genetics_inference/data/nucleo_div/cobb.windowed.pi",
                  header = TRUE, sep = "")
cobbv$panel <- "validation"
cobbv$breed <- "cobb"
rossv <- read.csv("3_population_genetics_inference/data/nucleo_div/ross.windowed.pi",
                  header = TRUE, sep = "")
rossv$panel <- "validation"
rossv$breed <- "ross"


# stats  -----------------
# combine tables to compare statistically
all <- rbind(cobbi, cobbe, cobbc, cobbd, rossi, rosse, rossc, rossd)
all$PI <- all$PI*100

all_good <- all %>%
  mutate(panel = factor(
    panel, levels = c("validation", "internal",
                      "external","combined",
                      "diverse")
  ))

stat.test <- all_good %>%
  group_by(breed) %>%
  t_test(PI ~ panel,  alternative = "two.sided",mu = 0,conf.level = 0.95) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance(cutpoints = c(0, 0.001, 0.05, 1),
                   symbols = c("**", "*", "ns")
                   )
stat.test

stat.test <- stat.test %>% add_xy_position(data = all_good, x = "panel")
stat.test <- stat.test[-c(2,3,5,8,9,11),]


# plot  -----------------
# plot all the tables: panels and validation samples
all <- rbind(cobbi, cobbe, cobbc, cobbd, cobbv,
             rossi, rosse, rossc, rossd, rossv)
all$PI <- all$PI*100

all_good <- all %>%
  mutate(panel = factor(
    panel, levels = c("validation", "internal",
                      "external","combined",
                      "diverse")
    ))

bxp <- ggboxplot(
  all_good,
  # panel
  x = "panel",
  y = "PI",
  facet.by = "breed",
  xlab = "",
  # color
  fill = "panel",
  palette = c("grey", "#E69F00", "#56B4E9", "#009E73","#9281e2"),
  alpha = 0.6
  )

bxp + stat_pvalue_manual(stat.test)

bxp + stat_pvalue_manual(stat.test) +
  rremove("x.text") + rremove("x.ticks") + rremove("legend")


# Validation samples   -----------------
nucdiv <- read.csv("3_population_genetics_inference/data/nucleo_div/HD-samples-nucdiv.csv",
                   header = TRUE, sep = ";") %>%
  clean_names()

nucdiv$ind2 <- nucdiv$individuals
nucdiv$ind2 <- sub("ca03_14", "0.75", nucdiv$ind2)
nucdiv$ind2 <- sub("ca19_18", "2.83", nucdiv$ind2)
nucdiv$ind2 <- sub("ca24_01", "0.5", nucdiv$ind2)
nucdiv$ind2 <- sub("cc19_18", "1.25", nucdiv$ind2)
nucdiv$ind2 <- sub("ca06_05", "2.29", nucdiv$ind2)
nucdiv$ind2 <- sub("ca13_03", "0.97", nucdiv$ind2)
nucdiv$ind2 <- sub("ca21_09", "0.28", nucdiv$ind2)
nucdiv$ind2 <- sub("cb20_03", "3.73", nucdiv$ind2)
nucdiv$ind2 <- sub("cc09_13", "1.64", nucdiv$ind2)
nucdiv$ind2 <- sub("cb07_06", "1.86", nucdiv$ind2)
nucdiv$ind2 <- as.character(nucdiv$ind2)

nucdiv %>%
  pivot_longer(HD:diverse, names_to = "method", values_to = "concordance") %>%
  mutate(Method = factor(Method, levels = c("HD","nternal",
                                            "external","combined","diverse")
                         )) %>%
  ggplot(aes(x = ind2, y = Concordance, color = Method)) +
  geom_point(size = 2) +
  facet_grid(~breed) +
  # text
  ylab("NucDiv") + xlab("Depth") +
  # color
  scale_color_manual(values = c("black", "#E69F00",
                                "#56B4E9", "#009E73","#9281e2")) +
  theme_bw() +
  theme(legend.position = "bottom")
