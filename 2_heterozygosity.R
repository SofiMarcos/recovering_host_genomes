# HETEROZYGOSITY #

library(tidyverse)
library(janitor)

# Load data  -----------------
# internal
internal_c <- read.csv("3_population_genetics_inference/data/het/internal-cobb.het",
                       header = TRUE, sep = "")
internal_r <- read.csv("3_population_genetics_inference/data/het/internal-ross.het",
                       header = TRUE, sep = "")
# external
external_c <- read.csv("3_population_genetics_inference/data/het/external-cobb.het",
                       header = TRUE, sep = "")
external_r <- read.csv("3_population_genetics_inference/data/het/external-ross.het",
                       header = TRUE, sep = "")
# combined
combined_c <- read.csv("3_population_genetics_inference/data/het/combined-cobb.het",
                       header = TRUE, sep = "")
combined_r <- read.csv("3_population_genetics_inference/data/het/combined-ross.het",
                       header = TRUE, sep = "")
# diverse
diverse_c <- read.csv("3_population_genetics_inference/data/het/diverse-cobb.het",
                      header = TRUE, sep = "")
diverse_r <- read.csv("3_population_genetics_inference/data/het/diverse-ross.het",
                      header = TRUE, sep = "")
# HD
HD_c <- read.csv("3_population_genetics_inference/data/het/HD-cobb.het",
                 header = TRUE, sep = "")
HD_r <- read.csv("3_population_genetics_inference/data/het/HD-ross.het",
                 header = TRUE, sep = "")

# Prepare data  -----------------
internal_c$O.HET. <- 1 - (internal_c$O.HOM./internal_c$N.NM.)
internal_c$E.HET. <- 1 - (internal_c$E.HOM./internal_c$N.NM.)
internal_c$Panel <- 'I'
internal_c$Breed <- 'Cobb'

internal_r$O.HET. <- 1 - (internal_r$O.HOM./internal_r$N.NM.)
internal_r$E.HET. <- 1 - (internal_r$E.HOM./internal_r$N.NM.)
internal_r$Panel <- 'I'
internal_r$Breed <- 'Ross'

external_c$O.HET. <- 1 - (external_c$O.HOM./external_c$N.NM.)
external_c$E.HET. <- 1 - (external_c$E.HOM./external_c$N.NM.)
external_c$Panel <- 'E'
external_c$Breed <- 'Cobb'

external_r$O.HET. <- 1 - (external_r$O.HOM./external_r$N.NM.)
external_r$E.HET. <- 1 - (external_r$E.HOM./external_r$N.NM.)
external_r$Panel <- 'E'
external_r$Breed <- 'Ross'

combined_c$O.HET. <- 1 - (combined_c$O.HOM./combined_c$N.NM.)
combined_c$E.HET. <- 1 - (combined_c$E.HOM./combined_c$N.NM.)
combined_c$Panel <- 'C'
combined_c$Breed <- 'Cobb'

combined_r$O.HET. <- 1 - (combined_r$O.HOM./combined_r$N.NM.)
combined_r$E.HET. <- 1 - (combined_r$E.HOM./combined_r$N.NM.)
combined_r$Panel <- 'C'
combined_r$Breed <- 'Ross'

diverse_c$O.HET. <- 1 - (diverse_c$O.HOM./diverse_c$N.NM.)
diverse_c$E.HET. <- 1 - (diverse_c$E.HOM./diverse_c$N.NM.)
diverse_c$Panel <- 'D'
diverse_c$Breed <- 'Cobb'

diverse_r$O.HET. <- 1 - (diverse_r$O.HOM./diverse_r$N.NM.)
diverse_r$E.HET. <- 1 - (diverse_r$E.HOM./diverse_r$N.NM.)
diverse_r$Panel <- 'D'
diverse_r$Breed <- 'Ross'

HD_c$O.HET. <- 1 - (HD_c$O.HOM./HD_c$N.NM.)
HD_c$E.HET. <- 1 - (HD_c$E.HOM./HD_c$N.NM.)
HD_c$Panel <- 'V'
HD_c$Breed <- 'Cobb'

HD_r$O.HET. <- 1 - (HD_r$O.HOM./HD_r$N.NM.)
HD_r$E.HET. <- 1 - (HD_r$E.HOM./HD_r$N.NM.)
HD_r$Panel <- 'V'
HD_r$Breed <- 'Ross'

HD <- rbind(HD_c, HD_r)

# stats  -----------------
# combine tables to compare
total <- rbind(internal_c, internal_r, external_c, external_r,
               combined_c, combined_r, diverse_c, diverse_r)

total <- as.data.frame(total)
total$O.HOM. <- NULL
total$E.HOM. <- NULL
total$F <- NULL
total$NM <- NULL
total$FID <- NULL
names(total)[1] <- "IND"

het <- total %>%
  mutate(Panel = factor(Panel, levels = c("V", "I","E","C","D")))

stat.test <- het %>%
  group_by(Breed) %>%
  t_test(O.HET. ~ Panel, alternative = "two.sided",
         mu = 0, paired = TRUE, conf.level = 0.95
  ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance(cutpoints = c(0, 0.001, 0.05, 1),
                   symbols = c("**", "*", "ns"))

stat.test <- stat.test %>% add_xy_position(data = het, x = "Panel")
stat.test <- stat.test[-c(2,3,5,8,9,11),]

# Plot  -----------------
# # combine all tables
total <- rbind(internal_c, internal_r, external_c, external_r,
               combined_c, combined_r, diverse_c, diverse_r, HD)

total <- as.data.frame(total)
total$O.HOM. <- NULL
total$E.HOM. <- NULL
total$F <- NULL
total$NM <- NULL
total$FID <- NULL
names(total)[1] <- "IND"

het <- total %>%
  mutate(Panel = factor(Panel, levels = c("V", "I","E","C","D")))

brp <- ggboxplot(
  het, x = "Panel", y = "O.HET.", facet.by = "Breed", fill = "Panel",
  palette = c("grey", "#E69F00", "#56B4E9", "#009E73","#9281e2"), alpha = 0.6,
  xlab = "")

brp + stat_pvalue_manual(stat.test)
brp + stat_pvalue_manual(stat.test) +
  rremove("x.text") + rremove("x.ticks")


# alternative plot - separate by breed
ggplot(total, aes(x = Panel, y = 100*O.HET., fill = Breed)) +
  geom_boxplot(size = 0.5, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(size = 1,
              position = position_jitterdodge(jitter.width = 0.2,
                                              dodge.width = 0
                                              )
              ) +
  facet_wrap(~Breed) + theme_bw() + ylab("O.Het") + xlab("Target population") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text.x = element_blank())


# Validation samples  -----------------
HD_c <- read.csv("3_population_genetics_inference/data/het/HD-cobb.het",
                 header = TRUE, sep = "")
HD_r <- read.csv("3_population_genetics_inference/data/het/HD-ross.het",
                 header = TRUE, sep = "")

HD_c$O.HET. <- 1 - (HD_c$O.HOM./HD_c$N.NM.)
HD_c$E.HET. <- 1 - (HD_c$E.HOM./HD_c$N.NM.)
HD_c$Panel <- 'True'
HD_c$Breed <- 'Cobb'

HD_r$O.HET. <- 1 - (HD_r$O.HOM./HD_r$N.NM.)
HD_r$E.HET. <- 1 - (HD_r$E.HOM./HD_r$N.NM.)
HD_r$Panel <- 'True'
HD_r$Breed <- 'Ross'

HD <- rbind(HD_c, HD_r)
HD$O.HOM. <- NULL
HD$E.HOM. <- NULL
HD$F <- NULL
HD$NM <- NULL
HD$FID <- NULL
names(HD)[1] <- "IND"

HD <- HD[-c(5, 6), ] # filter the ones that were excluded for the analysis

HD$IND <- sub("C1a", "F1a", HD$IND)
names <- HD$IND

total <- rbind(internal_c, internal_r,
               external_c, external_r,
               combined_c, combined_r,
               diverse_c, diverse_r
               )
total <- as.data.frame(total)
total$O.HOM. <- NULL
total$E.HOM. <- NULL
total$F <- NULL
total$NM <- NULL
total$FID <- NULL
names(total)[1] <- "IND"

validation_samples <- total[total$IND %in% names, ]

val <- rbind(HD, validation_samples)
val$IND2 <- val$IND
val$IND2 <- sub("CA03_14F1a", "0.75", val$IND2)
val$IND2 <- sub("CA19_18F1a", "2.83", val$IND2)
val$IND2 <- sub("CA24_01F1a", "0.5", val$IND2)
val$IND2 <- sub("CC19_18F1a", "1.25", val$IND2)
val$IND2 <- sub("CA06_05F1a", "2.29", val$IND2)
val$IND2 <- sub("CA13_03F1a", "0.97", val$IND2)
val$IND2 <- sub("CA21_09F1a", "0.28", val$IND2)
val$IND2 <- sub("CB20_03F1a", "3.73", val$IND2)
val$IND2 <- sub("CC09_13F1a", "1.64", val$IND2)
val$IND2 <- sub("CB07_06F1a", "1.86", val$IND2)
val$IND2 <- as.character(val$IND2)


# all together
val %>%
  mutate(Panel = factor(Panel, levels = c(
    "True", "Internal","External","Combined","Diverse")
    )) %>%
  ggplot(
    aes(x = IND2, y = 100*O.HET., color = Panel, shape = Breed)) +
  geom_point(size = 2) +
  # text
  xlab("depth") +
  # color
  scale_color_manual(values = c("black", "#E69F00", "#56B4E9",
                                "#009E73","#9281e2")
                     ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "top"
    ) +
  theme_bw()

# One box per breed
val %>%
  mutate(Panel = factor(Panel, levels = c(
    "True","Internal","External","Combined","Diverse")
    )) %>%
  ggplot(aes(x = IND2, y = 100*O.HET., color= Panel)) +
  geom_point(size = 3) +
  facet_wrap(~Breed) +
  # color
  scale_color_manual(values = c("black", "#E69F00", "#56B4E9",
                                "#009E73","#9281e2")) +

  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom"
        ) +
  theme_bw()
