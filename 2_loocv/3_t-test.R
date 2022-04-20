# T-TESTS #

library(ggpubr)
library(rstatix)
library(tidyverse)


# Concordance  -----------------
concordance <- read.csv("2_loocv/Leave-one-out-concordance.csv",
                        header = TRUE, sep = ";")
con <- concordance %>%
  pivot_longer(Internal:Diverse, names_to = "Panel", values_to = "Concordance") %>%
  mutate(Chr = factor(Chr, levels = c("GGA1","GGA7","GGA20"))) %>%
  mutate(Panel = factor(Panel, levels = c("Internal","External","Combined","Diverse")))

con$Panel <- as.factor(con$Panel)

con %>%
  group_by(Chr) %>%
  var.test(Concordance ~ Panel, ratio = 1, alternative = "two.sided", conf.level = 0.95) %>%
  add_significance()

stat.test <- con %>%
  group_by(Chr) %>%
  t_test(Concordance ~ Panel,
         alternative = "two.sided",
         mu = 0, paired = TRUE,
         var.equal = TRUE,
         conf.level = 0.95
         ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance(cutpoints = c(0, 0.001, 0.05, 1),
                    symbols = c("**", "*", "ns")
  )
stat.test

c <- con %>%
  mutate(Chr = factor(Chr, levels = c("GGA1","GGA7","GGA20"))) %>%
  mutate(Panel = factor(Panel, levels = c("Internal","External","Combined","Diverse"))) %>%
  ggplot(aes(x = Panel, y = Concordance, fill = Panel)) +
  geom_boxplot(size = 0.5, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(size = 1, position = position_jitterdodge(jitter.width = 1, dodge.width = 1)) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#9281e2")) +
  facet_grid(~Chr) +
  theme_bw() + ylab("Concordance")
c

bxp <- ggboxplot(
  con, x = "Panel", y = "Concordance", facet.by = "Chr", fill = "Panel",
  palette = c("#E69F00", "#56B4E9", "#009E73","#9281e2"), alpha = 0.6, add = "jitter",
  xlab = "")


# Make facet and add p-values
stat.test <- stat.test %>% add_xy_position(data = con, x = "Panel")
stat.test <- stat.test[-c(2,3,5,8,9,11,14,15,17),]
bxp + stat_pvalue_manual(stat.test)

final <- bxp + rremove("x.text") + rremove("x.ticks") + stat_pvalue_manual(stat.test)
final


# precision -----------------
precision <- read.csv("2_loocv/Leave-one-out-hetprecision.csv",
                      header = TRUE, sep = ";")
pre <- precision %>%
  pivot_longer(Internal:Diverse, names_to = "Panel", values_to = "Het.precision") %>%
  mutate(Chr = factor(Chr, levels = c("GGA1","GGA7","GGA20"))) %>%
  mutate(Panel = factor(Panel, levels = c("Internal","External","Combined","Diverse")))

stat.test <- pre %>%
  group_by(Chr) %>%
  t_test(Het.precision ~ Panel,
         alternative = "two.sided",
         mu = 0, paired = TRUE,
         var.equal = TRUE,
         conf.level = 0.95
         ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance(cutpoints = c(0, 0.001, 0.05, 1),
                   symbols = c("**", "*", "ns"))
stat.test

bnp <- ggboxplot(
  pre, x = "Panel", y = "Het.precision", facet.by = "Chr", fill = "Panel",
  palette = c("#E69F00", "#56B4E9", "#009E73","#9281e2"), alpha = 0.6, add = "jitter",
  xlab = "")


# Make facet and add p-values
stat.test <- stat.test %>% add_xy_position(data = pre, x = "Panel")
stat.test <- stat.test[-c(2,3,5,8,9,11,14,15,17),]
bnp + stat_pvalue_manual(stat.test)

prefinal <- bnp + stat_pvalue_manual(stat.test) +
  rremove("x.text") + rremove("x.ticks") + rremove("legend")
prefinal

library(cowplot)
e <- plot_grid(final,prefinal, labels = c("A", "B"), ncol = 1, nrow = 2, align = "h", axis = "bt")


# Minor alleles -----------------
maf <- read.csv("2_loocv/MAFs_precision.csv",header = TRUE, sep = ";")
pmaf <- maf %>%
  pivot_longer(Internal:Diverse, names_to = "Panel", values_to = "Het.precision") %>%
  mutate(Panel = factor(Panel, levels = c("Internal","External","Combined","Diverse"))) %>%
  mutate(MAF = factor(MAF, levels = c("0-0.05","0.05-0.1","0.1-0.3",">0.3")))

con$Panel <- as.factor(con$Panel)

stat.test <- pmaf %>%
  group_by(MAF) %>%
  t_test(Het.precision ~ Panel,
         alternative = "two.sided",
         mu = 0, paired = TRUE,
         var.equal = TRUE,
         conf.level = 0.95
         ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance(cutpoints = c(0, 0.001, 0.05, 1),
                   symbols = c("**", "*", "ns"))
stat.test

bxp <- ggboxplot(
  pmaf, x = "Panel", y = "Het.precision", fill = "Panel",
  palette = c("#E69F00", "#56B4E9", "#009E73","#9281e2"), alpha = 0.6, add = "jitter",
  xlab = "")

# Make facet and add p-values
stat.test <- stat.test %>% add_xy_position(data = pmaf, x = "Panel")
stat.test <- stat.test[-c(2,3,5,8,9,11,14,15,17,20,21,23),]
bxp + stat_pvalue_manual(stat.test)

test <- facet(bxp + theme_bw(), facet.by = "MAF", ncol = 4) +
  stat_pvalue_manual(stat.test)
maffinal <- ggpar(test,legend = "top") +
  rremove("x.text") + rremove("x.ticks") + rremove("grid") + rremove("legend")
maffinal

# Minor alleles number of variants -----------------
MAF <- read.csv("2_loocv/MAFs_final.csv", header = TRUE, sep = ";")

test <- MAF %>%
  pivot_longer(Internal:Diverse, names_to = "Panel", values_to = "Values") %>%
  mutate(Panel = factor(Panel, levels = c("Internal","External","Combined","Diverse"))) %>%
  mutate(MAF = factor(MAF, levels = c("0-0.05","0.05-0.1","0.1-0.3",">0.3")))

maf <- ggplot(test, aes(x = Panel, y = Values/100000,
                        fill = Panel, group = Imputed, alpha = Imputed)
              ) +
  geom_bar(position = "dodge", stat = "identity",color = "black") +
  facet_wrap(~MAF, nrow = 1, ncol = 4) +
  # color
  scale_fill_manual(values = c("#E69F00","#56B4E9","#009E73","#9281e2")) +
  scale_alpha_manual(values = c(0.5,1)) +
  theme_bw() + ylab("Number of SNPs (K)") + xlab("MAF bins") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  guides(group = guide_legend(nrow = 2, byrow = TRUE))

maf


library(cowplot)
f <- plot_grid(maffinal,maf, labels = c("A", "B"), ncol = 1, nrow = 2, align = "hv", axis = "bt")
