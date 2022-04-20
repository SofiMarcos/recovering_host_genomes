# FST #

library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(qqman)
library(ggsignif)


# Load data  -----------------
# internal
internal <- read.csv("3_population_genetics_inference/data/fst/internal.windowed.weir.fst",
                     header = TRUE, sep = "")

internal <- filter(internal, CHROM != "NC_006127.5")  # delete sex chromosome
internal$ZFST <- (internal$MEAN_FST - mean(internal$MEAN_FST)) / sd(internal$MEAN_FST)
internalq <- subset(internal, ZFST > quantile(internal$ZFST, 0.999))
internalm <- subset(internal, ZFST > quantile(internal$ZFST, 0.99))
internalq$pos <- paste(internalq$CHROM,"_", internalq$BIN_START)
internalm$pos <- paste(internalm$CHROM,"_", internalm$BIN_START)

# external
external <- read.csv("3_population_genetics_inference/data/fst/external.windowed.weir.fst",
                     header = TRUE, sep = "")

external <- filter(external, CHROM != "NC_006127.5")  # delete sex chromosome
external$ZFST <- (external$MEAN_FST - mean(external$MEAN_FST)) / sd(external$MEAN_FST)
externalq <- subset(external, ZFST > quantile(external$ZFST, 0.999))
externalm <- subset(external, ZFST > quantile(internal$ZFST, 0.99))
externalq$pos <- paste(externalq$CHROM,"_", externalq$BIN_START)
externalm$pos <- paste(externalm$CHROM,"_", externalm$BIN_START)

# combined
combined <- read.csv("3_population_genetics_inference/data/fst/combined.windowed.weir.fst",
                     header = TRUE, sep = "")

combined <- filter(combined, CHROM != "NC_006127.5")  # delete sex chromosome
combined$ZFST <- (combined$MEAN_FST - mean(combined$MEAN_FST)) / sd(combined$MEAN_FST)
combinedq <- subset(combined, ZFST > quantile(combined$ZFST, 0.999))
combinedm <- subset(combined, ZFST > quantile(internal$ZFST, 0.99))
combinedq$pos <- paste(combinedq$CHROM,"_", combinedq$BIN_START)
combinedm$pos <- paste(combinedm$CHROM,"_", combinedm$BIN_START)

# diverse
diverse <- read.csv("3_population_genetics_inference/data/fst/diverse.windowed.weir.fst",
                    header = TRUE, sep = "")

diverse <- filter(diverse, CHROM != "NC_006127.5")  # delete sex chromosome
diverse$ZFST <- (diverse$MEAN_FST - mean(diverse$MEAN_FST)) / sd(diverse$MEAN_FST)
diverseq <- subset(diverse, ZFST > quantile(diverse$ZFST, 0.999))
diversem <- subset(diverse, ZFST > quantile(internal$ZFST, 0.99))
diverseq$pos <- paste(diverseq$CHROM,"_", diverseq$BIN_START)
diversem$pos <- paste(diversem$CHROM,"_", diversem$BIN_START)

## common windows and number of variants
rownames(internalm) <- internalm$pos
rownames(externalm) <- externalm$pos
rownames(combinedm) <- combinedm$pos
rownames(diversem) <- diversem$pos

AlData <- list(internalm,externalm,combinedm,diversem)
common_names = Reduce(intersect, lapply(AlData, row.names))

intercommon <- internalm[internalm$pos %in% common_names,]
extercommon <- externalm[externalm$pos %in% common_names,]
combicommon <- combinedm[combinedm$pos %in% common_names,]
divercommon <- diversem[diversem$pos %in% common_names,]


# mean number of SNPs in the common windows
var.test(intercommon$N_VARIANTS, extercommon$N_VARIANTS, ratio = 1,
         alternative = "two.sided",conf.level = 0.95)
var.test(extercommon$N_VARIANTS, combicommon$N_VARIANTS, ratio = 1,
         alternative = "two.sided",conf.level = 0.95)
var.test(combicommon$N_VARIANTS, divercommon$N_VARIANTS, ratio = 1,
         alternative = "two.sided",conf.level = 0.95)

t.test(intercommon$N_VARIANTS, extercommon$N_VARIANTS, alternative = "two.sided",
       mu = 0, paired = TRUE, var.equal = TRUE,conf.level = 0.95)
t.test(extercommon$N_VARIANTS, combicommon$N_VARIANTS, alternative = "two.sided",
       mu = 0, paired = TRUE, var.equal = TRUE,conf.level = 0.95)
t.test(combicommon$N_VARIANTS, divercommon$N_VARIANTS, alternative = "two.sided",
       mu = 0, paired = TRUE, var.equal = TRUE,conf.level = 0.95)

common_windows <- as.data.frame(cbind(
  intercommon$N_VARIANTS, extercommon$N_VARIANTS,
  combicommon$N_VARIANTS, divercommon$N_VARIANTS
  ))
colnames(common_windows) <- c("internal", "external", "combined", "diverse")

m <- common_windows %>%
  pivot_longer(internal:diverse, names_to = "Method", values_to = "Variants") %>%
  mutate(Method = factor(Method,
                         levels = c("internal","external","combined","diverse")
                         )) %>%
  # plot
  ggplot(aes(x = Method, y = Variants, fill = Method)) +
  geom_violin( alpha = 0.6) +
  # text
  ylab("Number of variants per window") +
  # color
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#9281e2")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
        ,axis.text = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        strip.text.x = element_text(size = 9)
        )
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)
         )
m


m <- common_windows %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(internal:diverse, names_to = "Method", values_to = "Variants") %>%
  mutate(Method = factor(Method,
                         levels = c("internal","external","combined","diverse")
                         )) %>%
  ggplot(aes(x = Method, y = Variants, fill = Method)) +
  geom_segment( aes(x = Method, xend = Method, y = 0, yend = Variants)) +
  geom_point(size = 5, alpha = 0.7, shape = 21, stroke = 2) +
  # text
  ylab("Number of variants per window") +
  # color
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#9281e2")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 11,face = "bold"),
        strip.text.x = element_text(size = 9)
  )
m

# Venn diagrams-----------------
library("grid")
library("ggVennDiagram")

all <- list(internalq$pos, externalq$pos, combinedq$pos, diverseq$pos)
names(all) <- c("Internal","External","Combined", "Diverse")

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(all, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(x,
             category.names = c("" , "" , "", ""),
             # Circles
             lwd = 2,lty = 'blank',
             fill = c("#E69F00", "#56B4E9", "#009E73","#9281e2"),
             # Numbers
             cex = 1.5,
             # Set names
             cat.cex = 1.5,
             cat.fontfamily = "Arial",
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.055, 0.055, 0.1, 0.1),
             fontfamily = "Arial"
             )


all <- list(internalm$pos, externalm$pos, combinedm$pos, diversem$pos)
names(all) <- c("Internal","External","Combined", "Diverse")

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(all, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(x,
             category.names = c("" , "" , "", ""),
             # Circles
             lwd = 3,lty = 'blank',
             fill = c("#E69F00", "#56B4E9", "#009E73","#9281e2"),
             # Numbers
             cex = 1.5,
             # Set names
             cat.cex = 1.5,
             cat.fontfamily = "Arial",
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.055, 0.055, 0.1, 0.1),
             fontfamily = "Arial"
             )

library(cowplot)
plot_grid(d,q,m,
          labels = c("E","", "F"),
          ncol = 3,
          nrow = 1,
          align = "h",
          axis = "bt"
          )
