# IBS DISTANCE HEATMAP #

library(stringr)
library(gplots)
library(RColorBrewer)
library(janitor)

# Load data - population - all -----------------
# internal
internal <- read.csv("3_population_genetics_inference/data/ibs/internal-maf.mibs",
                     header = FALSE, sep = "")
internal_id <- read.csv("3_population_genetics_inference/data/ibs/internal-maf.mibs.id",
                       header = FALSE, sep = "")

rownames(internal) <- internal_id$V1
colnames(internal) <- internal_id$V1

internal <- 1 - internal
internal[internal == 0] <- NA
internal <- as.matrix(internal)
internal2 <- as.dist(internal)

internal_id$V2 <- c("Cobb","Ross","Cobb","Cobb","Cobb","Cobb",
                    "Ross","Ross","Cobb","Ross","Ross","Ross",
                    "Ross","Cobb","Cobb","Ross","Ross","Ross",
                    "Ross","Cobb","Cobb","Ross","Cobb","Ross",
                    "Ross","Ross","Ross","Cobb","Cobb","Cobb",
                    "Ross","Ross","Ross","Ross","Cobb","Ross",
                    "Ross","Cobb","Cobb","Cobb","Cobb","Cobb",
                    "Ross","Ross","Ross","Cobb","Ross","Ross",
                    "Cobb","Cobb","Cobb","Cobb","Ross","Ross",
                    "Ross","Cobb","Cobb","Ross","Ross","Ross",
                    "Ross","Ross","Ross","Ross","Cobb","Cobb",
                    "Cobb","Cobb","Cobb","Cobb","Cobb","Cobb",
                    "Ross","Ross","Cobb","Ross","Ross","Ross",
                    "Ross","Ross","Ross","Ross","Ross","Ross",
                    "Cobb","Cobb","Cobb","Cobb","Cobb","Cobb",
                    "Cobb","Cobb","Ross","Cobb","Cobb","Ross",
                    "Cobb","Ross","Ross","Ross"
)

# plot
palette.breaks <- c(0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32)
col1 <- c("#F8766D","#00BFC4")

heatmap.2(internal, col = brewer.pal(8, "BuPu"),  # change panel
          trace = "none",
          na.color = "grey",
          symm = TRUE,
          # dendrogram
          dendrogram = "column",
          # key title
          density.info = "none",
          key.title = "none",
          key.xlab = "Pairwise distance",
          keysize = 1.5,
          # row and columns labels size
          cexRow = 0.75,cexCol = 0.5,
          labRow = FALSE, labCol = FALSE,
          # graphic margins
          margins = c(5,7),
          # control orientation
          revC = TRUE,na.rm = TRUE,
          breaks = palette.breaks,
          ColSideColors = col1[as.factor(internal_id$V2)]  # change panel_id
)


# external
external <- read.csv("3_population_genetics_inference/data/ibs/external-maf.mibs",
                     header = FALSE, sep = "")
external_id <- read.csv("3_population_genetics_inference/data/ibs/external-maf.mibs.id",
                        header = FALSE, sep = "")
rownames(external) <- external_id$V1
colnames(external) <- external_id$V1

external <- 1 - external
external[external == 0] <- NA
external <- as.matrix(external)
external2 <- as.dist(external)

external_id$V2 <- c("Cobb","Ross","Cobb","Cobb","Cobb","Cobb",
                    "Ross","Ross","Cobb","Ross","Ross","Ross",
                    "Ross","Cobb","Cobb","Ross","Ross","Ross",
                    "Ross","Cobb","Cobb","Ross","Cobb","Ross",
                    "Ross","Ross","Ross","Cobb","Cobb","Cobb",
                    "Ross","Ross","Ross","Ross","Cobb","Ross",
                    "Ross","Cobb","Cobb","Cobb","Cobb","Cobb",
                    "Ross","Ross","Ross","Cobb","Ross","Ross",
                    "Cobb","Cobb","Cobb","Cobb","Ross","Ross",
                    "Ross","Cobb","Cobb","Ross","Ross","Ross",
                    "Ross","Ross","Ross","Ross","Cobb","Cobb",
                    "Cobb","Cobb","Cobb","Cobb","Cobb","Cobb",
                    "Ross","Ross","Cobb","Ross","Ross","Ross",
                    "Ross","Ross","Ross","Ross","Ross","Ross",
                    "Cobb","Cobb","Cobb","Cobb","Cobb","Cobb",
                    "Cobb","Cobb","Ross","Cobb","Cobb","Ross",
                    "Cobb","Ross","Ross","Ross"
)

# plot
palette.breaks <- c(0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32)
col1 <- c("#F8766D","#00BFC4")

heatmap.2(external, col = brewer.pal(8, "BuPu"),  # change panel
          trace = "none",
          na.color = "grey",
          symm = TRUE,
          # dendrogram
          dendrogram = "column",
          # key title
          density.info = "none",
          key.title = "none",
          key.xlab = "Pairwise distance",
          keysize = 1.5,
          # row and columns labels size
          cexRow = 0.75,cexCol = 0.5,
          labRow = FALSE, labCol = FALSE,
          # graphic margins
          margins = c(5,7),
          # control orientation
          revC = TRUE,na.rm = TRUE,
          breaks = palette.breaks,
          ColSideColors = col1[as.factor(external_id$V2)]  # change panel_id
)

# combined
combined <- read.csv("3_population_genetics_inference/data/ibs/combined-maf.mibs",
                     header = FALSE, sep = "")
combined_id <- read.csv("3_population_genetics_inference/data/ibs/combined-maf.mibs.id",
                       header = FALSE, sep = "")
rownames(combined) <- combined_id$V1
colnames(combined) <- combined_id$V1

combined <- 1 - combined
combined[combined == 0] <- NA
combined <- as.matrix(combined)
combined2 <- as.dist(combined)

combined_id$V2 <- c("Cobb","Ross","Cobb","Cobb","Cobb","Cobb",
                   "Ross","Ross","Cobb","Ross","Ross","Ross",
                   "Ross","Cobb","Cobb","Ross","Ross","Ross",
                   "Ross","Cobb","Cobb","Ross","Cobb","Ross",
                   "Ross","Ross","Ross","Cobb","Cobb","Cobb",
                   "Ross","Ross","Ross","Ross","Cobb","Ross",
                   "Ross","Cobb","Cobb","Cobb","Cobb","Cobb",
                   "Ross","Ross","Ross","Cobb","Ross","Ross",
                   "Cobb","Cobb","Cobb","Cobb","Ross","Ross",
                   "Ross","Cobb","Cobb","Ross","Ross","Ross",
                   "Ross","Ross","Ross","Ross","Cobb","Cobb",
                   "Cobb","Cobb","Cobb","Cobb","Cobb","Cobb",
                   "Ross","Ross","Cobb","Ross","Ross","Ross",
                   "Ross","Ross","Ross","Ross","Ross","Ross",
                   "Cobb","Cobb","Cobb","Cobb","Cobb","Cobb",
                   "Cobb","Cobb","Ross","Cobb","Cobb","Ross",
                   "Cobb","Ross","Ross","Ross"
)

# plot
palette.breaks <- c(0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32)
col1 <- c("#F8766D","#00BFC4")

heatmap.2(combined, col = brewer.pal(8, "BuPu"),  # change panel
          trace = "none",
          na.color = "grey",
          symm = TRUE,
          # dendrogram
          dendrogram = "column",
          # key title
          density.info = "none",
          key.title = "none",
          key.xlab = "Pairwise distance",
          keysize = 1.5,
          # row and columns labels size
          cexRow = 0.75,cexCol = 0.5,
          labRow = FALSE, labCol = FALSE,
          # graphic margins
          margins = c(5,7),
          # control orientation
          revC = TRUE,na.rm = TRUE,
          breaks = palette.breaks,
          ColSideColors = col1[as.factor(combined_id$V2)]  # change panel_id
)


# diverse
diverse <- read.csv("3_population_genetics_inference/data/ibs/diverse-maf.mibs",
                    header = FALSE, sep = "")
diverse_id <- read.csv("3_population_genetics_inference/data/ibs/diverse-maf.mibs.id",
                       header = FALSE, sep = "")
rownames(diverse) <- diverse_id$V1
colnames(diverse) <- diverse_id$V1

diverse <- 1 - diverse
diverse[diverse == 0] <- NA
diverse <- as.matrix(diverse)
diverse2 <- as.dist(diverse)

diverse_id$V2 <- c("Cobb","Ross","Cobb","Cobb","Cobb","Cobb",
                   "Ross","Ross","Cobb","Ross","Ross","Ross",
                   "Ross","Cobb","Cobb","Ross","Ross","Ross",
                   "Ross","Cobb","Cobb","Ross","Cobb","Ross",
                   "Ross","Ross","Ross","Cobb","Cobb","Cobb",
                   "Ross","Ross","Ross","Ross","Cobb","Ross",
                   "Ross","Cobb","Cobb","Cobb","Cobb","Cobb",
                   "Ross","Ross","Ross","Cobb","Ross","Ross",
                   "Cobb","Cobb","Cobb","Cobb","Ross","Ross",
                   "Ross","Cobb","Cobb","Ross","Ross","Ross",
                   "Ross","Ross","Ross","Ross","Cobb","Cobb",
                   "Cobb","Cobb","Cobb","Cobb","Cobb","Cobb",
                   "Ross","Ross","Cobb","Ross","Ross","Ross",
                   "Ross","Ross","Ross","Ross","Ross","Ross",
                   "Cobb","Cobb","Cobb","Cobb","Cobb","Cobb",
                   "Cobb","Cobb","Ross","Cobb","Cobb","Ross",
                   "Cobb","Ross","Ross","Ross"
)

# plot
palette.breaks <- c(0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32)
col1 <- c("#F8766D","#00BFC4")

heatmap.2(diverse, col = brewer.pal(8, "BuPu"),  # change panel
          trace = "none",
          na.color = "grey",
          symm = TRUE,
          # dendrogram
          dendrogram = "column",
          # key title
          density.info = "none",
          key.title = "none",
          key.xlab = "Pairwise distance",
          keysize = 1.5,
          # row and columns labels size
          cexRow = 0.75,cexCol = 0.5,
          labRow = FALSE, labCol = FALSE,
          # graphic margins
          margins = c(5,7),
          # control orientation
          revC = TRUE,na.rm = TRUE,
          breaks = palette.breaks,
          ColSideColors = col1[as.factor(diverse_id$V2)]  # change panel_id
)

# HD   -----------------
HD <- read.csv("3_population_genetics_inference/data/ibs/HD-maf.mibs",
               header = FALSE, sep = "")
HD_id <- read.csv("3_population_genetics_inference/data/ibs/HD-maf.mibs.id",
                  header = FALSE, sep = "")
rownames(HD) <- HD_id$V1
colnames(HD) <- HD_id$V1

HD <- 1 - HD
HD[HD == 0] <- NA
HD <- as.matrix(HD)
HD2 <- as.dist(HD)
HD_id$V2 <- c("Cobb", "Cobb", "Cobb", "Ross", "Ross",
              "Ross", "Cobb", "Ross", "Ross", "Ross")

# # Plot 12 validation samples
palette.breaks <- c(0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32)
col1 <- c("#F8766D","#00BFC4")

heatmap.2(combined, col = brewer.pal(6, "BuPu"),
          trace = "none",
          na.color = "grey",
          symm = TRUE,
          # dendrogram
          dendrogram = "column",
          # key title
          density.info = "none",
          key.title = "Histogram",
          key.xlab = "Pairwise distance",
          keysize = 1.5,
          # row and columns labels size
          #labRow = FALSE, labCol = FALSE,
          # graphic margins
          margins = c(1,1),
          # control orientation
          revC = TRUE,na.rm = TRUE,
          breaks = palette.breaks,
          ColSideColors = col1[as.factor(HD_id$V2)])


# Load data - population - by breed -----------------
# Ross  -----------------
# internal
internal <- read.csv("3_population_genetics_inference/data/ibs/internal-ross-maf.mibs",
                     header = FALSE, sep = "")
internal_id <- read.csv("3_population_genetics_inference/data/ibs/internal-ross-maf.mibs.id",
                       header = FALSE, sep = "")

rownames(internal) <- internal_id$V1
colnames(internal) <- internal_id$V1

internal <- 1 - internal
internal[internal == 0] <- NA
internal <- as.matrix(internal)
internal2 <- as.dist(internal)


# external
external <- read.csv("3_population_genetics_inference/data/ibs/external-ross-maf.mibs",
                     header = FALSE, sep = "")
external_id <- read.csv("3_population_genetics_inference/data/ibs/external-ross-maf.mibs.id",
                       header = FALSE, sep = "")

rownames(external) <- external_id$V1
colnames(external) <- external_id$V1

external <- 1 - external
external[external == 0] <- NA
external <- as.matrix(external)
external2 <- as.dist(external)


# combined
combined <- read.csv("3_population_genetics_inference/data/ibs/combined-ross-maf.mibs",
                     header = FALSE, sep = "")
combined_id <- read.csv("3_population_genetics_inference/data/ibs/combined-ross-maf.mibs.id",
                       header = FALSE, sep = "")
rownames(combined) <- combined_id$V1
colnames(combined) <- combined_id$V1

combined <- 1 - combined
combined[combined == 0] <- NA
combined <- as.matrix(combined)
combined2 <- as.dist(combined)


# diverse
diverse <- read.csv("3_population_genetics_inference/data/ibs/diverse-ross-maf.mibs",
                    header = FALSE, sep = "")
diverse_id <- read.csv("3_population_genetics_inference/data/ibs/diverse-ross-maf.mibs.id",
                      header = FALSE, sep = "")
rownames(diverse) <- diverse_id$V1
colnames(diverse) <- diverse_id$V1

diverse <- 1 - diverse
diverse[diverse == 0] <- NA
diverse <- as.matrix(diverse)
diverse2 <- as.dist(diverse)

# Cobb  -----------------
# internal
internal <- read.csv("3_population_genetics_inference/data/ibs/internal-cobb-maf.mibs",
                     header = FALSE, sep = "")
internal_id <- read.csv("3_population_genetics_inference/data/ibs/internal-cobb-maf.mibs.id",
                        header = FALSE, sep = "")

rownames(internal) <- internal_id$V1
colnames(internal) <- internal_id$V1

internal <- 1 - internal
internal[internal == 0] <- NA
internal <- as.matrix(internal)
internal2 <- as.dist(internal)


# external
external <- read.csv("3_population_genetics_inference/data/ibs/external-cobb-maf.mibs",
                     header = FALSE, sep = "")
external_id <- read.csv("3_population_genetics_inference/data/ibs/external-cobb-maf.mibs.id",
                        header = FALSE, sep = "")

rownames(external) <- external_id$V1
colnames(external) <- external_id$V1

external <- 1 - external
external[external == 0] <- NA
external <- as.matrix(external)
external2 <- as.dist(external)


# combined
combined <- read.csv("3_population_genetics_inference/data/ibs/combined-cobb-maf.mibs",
                     header = FALSE, sep = "")
combined_id <- read.csv("3_population_genetics_inference/data/ibs/combined-cobb-maf.mibs.id",
                        header = FALSE, sep = "")
rownames(combined) <- combined_id$V1
colnames(combined) <- combined_id$V1

combined <- 1 - combined
combined[combined == 0] <- NA
combined <- as.matrix(combined)
combined2 <- as.dist(combined)


# diverse
diverse <- read.csv("3_population_genetics_inference/data/ibs/diverse-cobb-maf.mibs",
                    header = FALSE, sep = "")
diverse_id <- read.csv("3_population_genetics_inference/data/ibs/diverse-cobb-maf.mibs.id",
                       header = FALSE, sep = "")
rownames(diverse) <- diverse_id$V1
colnames(diverse) <- diverse_id$V1

diverse <- 1 - diverse
diverse[diverse == 0] <- NA
diverse <- as.matrix(diverse)
diverse2 <- as.dist(diverse)

# plot
palette.breaks <- c(0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32)
col1 <- c("#F8766D","#00BFC4")

heatmap.2(diverse, col = brewer.pal(8, "BuPu"),  # change panel
          trace = "none",
          na.color = "grey",
          symm = TRUE,
          # dendrogram
          dendrogram = "column",
          # key title
          density.info = "none",
          key.title = "none",
          key.xlab = "Pairwise distance",
          keysize = 1.5,
          # row and columns labels size
          cexRow = 0.75,cexCol = 0.5,
          labRow = FALSE, labCol = FALSE,
          # graphic margins
          margins = c(5,7),
          # control orientation
          revC = TRUE,na.rm = TRUE,
          breaks = palette.breaks,
          ColSideColors = col1[as.factor(diverse_id$V2)]  # change panel_id
          )


# Mantel test  -----------------
library(ade4)
# validation samples
mantel.rtest(HD2, internal2, nrepet = 9999)
mantel.rtest(HD2, external2, nrepet = 9999)
mantel.rtest(HD2, combined2, nrepet = 9999)
mantel.rtest(HD2, diverse2, nrepet = 9999)


# population
mantel.rtest(internal2, external2, nrepet = 9999)
mantel.rtest(external2, combined2, nrepet = 9999)
mantel.rtest(combined2, diverse2, nrepet = 9999)

mantel.rtest(internal2, combined2, nrepet = 9999)
mantel.rtest(internal2, diverse2, nrepet = 9999)
mantel.rtest(external2, diverse2, nrepet = 9999)
