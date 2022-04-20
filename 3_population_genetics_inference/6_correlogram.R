# CORRELOGRAM #

library(ggcorrplot)
library(tidyverse)

# Load data  -----------------
# Ross kinship
kinr <- read.csv("3_population_genetics_inference/data/ibs/CorrelationKinship-ross.csv",
                 sep = ";")
rownames(kinr) <- kinr[,1]
kinr[,1] <- NULL

# Cobb kinship
kinc <- read.csv("3_population_genetics_inference/data/ibs/CorrelationKinship-cobb.csv",
                 sep = ";")
rownames(kinc) <- kinc[,1]
kinc[,1] <- NULL

# Ross ibs
ibsr <- read.csv("3_population_genetics_inference/data/ibs/CorrelationIBS-ross.csv",
                 sep = ";")
rownames(ibsr) <- ibsr[,1]
ibsr[,1] <- NULL

# Cobb ibs
ibsc <- read.csv("3_population_genetics_inference/data/ibs/CorrelationIBS-cobb.csv",
                 sep = ";")
rownames(ibsc) <- ibsc[,1]
ibsc[,1] <- NULL


# Plot  -----------------
library(RColorBrewer)
kinplotr <- ggcorrplot(kinr,
                       hc.order = FALSE,
                       type = "lower",
                       outline.color = "white",
                       lab = TRUE
           ) +
  # text
  ylab("") + xlab("") + ggtitle("Ross") +
  # color
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"
      ) +
  theme_bw()

kinplotr
kinplotrfinal <- kinplotr + scale_fill_gradient2(limit = c(0.7,1),
                                                 low = "#3f7f93",
                                                 high =  "#c25539",
                                                 mid = "white",
                                                 midpoint = 0.88
                                                 )
kinplotrfinal


kinplotc <- ggcorrplot(kinc,
                       hc.order = FALSE,
                       type = "lower",
                       outline.color = "white",
                       lab = TRUE
                       ) +
  # text
  ylab("") + xlab("") + ggtitle("Cobb") +
  # color
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
) +
  theme_bw()

kinplotc
kinplotcfinal <- kinplotc + scale_fill_gradient2(limit = c(0.7,1),
                                                 low = "#3f7f93",
                                                 high =  "#c25539",
                                                 mid = "white",
                                                 midpoint = 0.88
                                                 )
kinplotcfinal

library(cowplot)
plot_grid(kinplotcfinal, kinplotrfinal,
          labels = c("C",""),
          ncol = 2, nrow = 1,
          align = "hv",
          axis = "bw"
          )

pwc <- ggcorrplot(ibsc,
                  hc.order = FALSE,
                  type = "lower",
                  outline.color = "white",
                  lab = TRUE
                  ) +
  # text
  ylab("") + xlab("") + ggtitle("Cobb") +
  # color
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
        ) +
  theme_bw()

pwcfinal <- pwc + scale_fill_gradient2(limit = c(0.7,1),
                                       low = "#3f7f93",
                                       high =  "#c25539",
                                       mid = "white",
                                       midpoint = 0.88
                                       )
pwcfinal

pwr <- ggcorrplot(ibsr,
                  hc.order = FALSE,
                  type = "lower",
                  outline.color = "white",
                  lab = TRUE
                  ) +
  # text
  ylab("") + xlab("") + ggtitle("Ross") +
  # color
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"
        ) +
  theme_bw()

pwrfinal <- pwr + scale_fill_gradient2(limit = c(0.7,1),
                                       low = "#3f7f93",
                                       high =  "#c25539",
                                       mid = "white",
                                       midpoint = 0.88
                                       )
pwrfinal

library(cowplot)
plot_grid(kinplotcfinal, kinplotrfinal,pwcfinal, pwrfinal,
          labels = c("C","","D",""),
          ncol = 2, nrow = 2,
          align = "hv",
          axis = "bw"
          )



# VALIDATION SAMPLES
# --------------------------------
# kinship
kin <- read.csv("3_population_genetics_inference/data/ibs/CorrHD-kin.csv"
                , sep = ";")
rownames(kin) <- kin[,1]
kin[,1] <- NULL

# ibs
ibs <- read.csv("3_population_genetics_inference/data/ibs/CorrHD-ibs.csv"
                , sep = ";")
rownames(ibs) <- ibs[,1]
ibs[,1] <- NULL

kin1 <- ggcorrplot(kin,
                   hc.order = FALSE,
                   type = "lower",
                   outline.color = "white",
                   lab = TRUE
                  ) +
  theme_bw() +
  ylab("") + xlab("") + ggtitle("Kinship") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"

        )
kin3 <- kin1 + scale_fill_gradient2(limit = c(0.7,1),
                                    low = "#6D9EC1",
                                    high =  "c25539",
                                    mid = "white",
                                    midpoint = 0.88
                                    )
kin3

pw <- ggcorrplot(ibs,
                 hc.order = FALSE,
                 type = "lower",
                 outline.color = "white",
                 lab = TRUE
                 ) +
  # text
  ylab("") + xlab("") + ggtitle("Pairwise distance") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"
        ) +
  theme_bw()

pw2 <- pw + scale_fill_gradient2(limit = c(0.7,1),
                                 low = "#6D9EC1",
                                 high =  "c25539",
                                 mid = "white",
                                 midpoint = 0.88
                                 )
pw2

library(cowplot)
plot_grid(pw2, kin3,
          labels = c("C","D"),
          ncol = 2,
          nrow = 1,
          align = "hv",
          axis = "bw"
          )
