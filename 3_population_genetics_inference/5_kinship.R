# KINSHIP HEATMAP #

library(stringr)
library(gplots)
library(RColorBrewer)

# Load data  -----------------
# Ross
# internal
internal <- read.csv("3_population_genetics_inference/data/kin/internal-maf.king",
                     header = FALSE, sep = "")
internal_id <- read.csv("3_population_genetics_inference/data/kin/internal-maf.king.id",
                        header = TRUE, sep = "")

rownames(internal) <- internal_id$IID
colnames(internal) <- internal_id$IID

internal[internal == 0.5000000] <- NA
internal[internal < 0] <- 0
internal <- as.matrix(internal)
internal2 <- as.dist(internal)

# external
external <- read.csv("3_population_genetics_inference/data/kin/external-maf.king",
                    header = FALSE, sep = "")
external_id <- read.csv("3_population_genetics_inference/data/kin/external-maf.king.id",
                       header = TRUE, sep = "")

rownames(external) <- external_id$IID
colnames(external) <- external_id$IID

external[external == 0.5000000] <- NA
external[external < 0] <- 0
external <- as.matrix(external)
external2 <- as.dist(external)

# combined
combined <- read.csv("3_population_genetics_inference/data/kin/combined-maf.king",
                     header = FALSE, sep = "")
combined_id <- read.csv("3_population_genetics_inference/data/kin/combined-maf.king.id",
                        header = TRUE, sep = "")

rownames(combined) <- combined_id$IID
colnames(combined) <- combined_id$IID

combined[combined == 0.5000000] <- NA
combined[combined < 0] <- 0
combined <- as.matrix(combined)
combined2 <- as.dist(combined)

# diverse
diverse <- read.csv("3_population_genetics_inference/data/kin/diverse-maf.king",
                   header = FALSE, sep = "")
diverse_id <- read.csv("3_population_genetics_inference/data/kin/diverse-maf.king.id",
                      header = TRUE, sep = "")
rownames(diverse) <- diverse_id$IID
colnames(diverse) <- diverse_id$IID

diverse[diverse == 0.5000000] <- NA
diverse[diverse < 0] <- 0
diverse <- as.matrix(diverse)
diverse2 <- as.dist(diverse)

# ross samples
#dist<- read.csv("ross.dist", header = FALSE, sep = "")
#samples_id<- read.csv("ross.dist.id", header = FALSE, sep = "")


# Plot  -----------------
palette.breaks <- c(0.015,0.041,0.067,0.093,0.119,0.145,0.171)
diverse_id$V2 <- c("Cobb","Ross","Cobb","Cobb","Cobb","Cobb","Ross",
                   "Ross","Cobb","Ross","Ross","Ross","Ross","Cobb",
                   "Cobb","Ross","Ross","Ross","Ross","Cobb","Cobb",
                   "Ross","Cobb","Ross","Ross","Ross","Ross","Cobb",
                   "Cobb","Cobb","Ross","Ross","Ross","Ross","Cobb",
                   "Ross","Ross","Cobb","Cobb","Cobb","Cobb","Cobb",
                   "Ross","Ross","Ross","Cobb","Ross","Ross","Cobb",
                   "Cobb","Cobb","Cobb","Ross","Ross","Ross","Cobb",
                   "Cobb","Ross","Ross","Ross","Ross","Ross","Ross",
                   "Ross","Cobb","Cobb","Cobb","Cobb","Cobb","Cobb",
                   "Cobb","Cobb","Ross","Ross","Cobb","Ross","Ross",
                   "Ross","Ross","Ross","Ross","Ross","Ross","Ross",
                   "Cobb","Cobb","Cobb","Cobb","Cobb","Cobb","Cobb",
                   "Cobb","Ross","Cobb","Cobb","Ross","Cobb","Ross",
                   "Ross","Ross"
                   )
col1 <- c("#F8766D","#00BFC4")

heatmap.2(diverse, col = brewer.pal(6, "OrRd"),
          trace = "none",
          na.color = "grey",
          symm = TRUE,
          # dendrogram
          dendrogram = "none",
          # key title
          density.info = "none",
          key.title = "Histogram",
          key.xlab = "Kinship",
          keysize = 1.5,
          # row and columns labels size
          #labRow = FALSE, labCol = FALSE,
          # graphic margins
          margins = c(1,1),
          # control orientation
          revC = TRUE, na.rm = TRUE,
          breaks = palette.breaks,
          ColSideColors = col1[as.factor(diverse_id$V2)]
          )


# Validation samples  -----------------
# internal
internal <- read.csv("3_population_genetics_inference/data/kin/internal-HD-maf.king",
                     header = FALSE, sep = "")
internal_id <- read.csv("3_population_genetics_inference/data/kin/internal-HD-maf.king.id",
                        header = TRUE, sep = "")
rownames(internal) <- internal_id$IID
colnames(internal) <- internal_id$IID

internal[internal == 0.5000000] <- NA
internal[internal < 0] <- 0
internal <- as.matrix(internal)
internal2 <- as.dist(internal)


# external
external <- read.csv("3_population_genetics_inference/data/kin/external-HD-maf.king",
                     header = FALSE, sep = "")
external_id <- read.csv("3_population_genetics_inference/data/kin/external-HD-maf.king.id",
                        header = TRUE, sep = "")
rownames(external) <- external_id$IID
colnames(external) <- external_id$IID

external[external == 0.5000000] <- NA
external[external < 0] <- 0
external <- as.matrix(external)
external2 <- as.dist(external)


# combined
combined <- read.csv("3_population_genetics_inference/data/kin/combined-HD-maf.king",
                     header = FALSE, sep = "")
combined_id <- read.csv("3_population_genetics_inference/data/kin/combined-HD-maf.king.id",
                        header = TRUE, sep = "")
rownames(combined) <- combined_id$IID
colnames(combined) <- combined_id$IID

combined[combined == 0.5000000] <- NA
combined[combined < 0] <- 0
combined <- as.matrix(combined)
combined2 <- as.dist(combined)


# diverse
diverse <- read.csv("3_population_genetics_inference/data/kin/diverse-HD-maf.king",
                   header = FALSE, sep = "")
diverse_id <- read.csv("3_population_genetics_inference/data/kin/diverse-HD-maf.king.id",
                      header = TRUE, sep = "")
rownames(diverse) <- diverse_id$IID
colnames(diverse) <- diverse_id$IID

diverse[diverse == 0.5000000] <- NA
diverse[diverse < 0] <- 0
diverse <- as.matrix(diverse)
diverse2 <- as.dist(diverse)

# HD
HD <- read.csv("3_population_genetics_inference/data/kin/HD-maf.king",
               header = FALSE, sep = "")
HD_id <- read.csv("3_population_genetics_inference/data/kin/HD-maf.king.id",
                 header = TRUE, sep = "")
rownames(HD) <- HD_id$IID
colnames(HD) <- HD_id$IID
HD[HD == 0.5000000] <- NA
HD[HD < 0] <- 0
HD <- as.matrix(HD)
HD2 <- as.dist(HD)

# plot
palette.breaks <- c(0, 0.013, 0.026, 0.039, 0.053,
                    0.066, 0.079, 0.092, 0.105,0.118
                    )
HD_id$V2 <- c("Cobb", "Cobb", "Cobb", "Ross", "Ross",
              "Ross", "Cobb", "Ross", "Ross", "Ross"
              )
col1 <- c("#F8766D","#00BFC4")

# detailed heatmap
heatmap.2(diverse, col = brewer.pal(6, "OrRd"),
          trace = "none",
          na.color = "grey",
          symm = TRUE,
          # dendrogram
          dendrogram = "none",
          # key title
          density.info = "none",
          key.title = "Histogram",
          key.xlab = "Kinship",
          keysize = 1.5,
          # row and columns labels size
          #labRow = FALSE, labCol = FALSE,
          # graphic margins
          margins = c(1,1),
          # control orientation
          revC = TRUE,na.rm = TRUE,
          breaks = palette.breaks,
          ColSideColors = col1[as.factor(HD_id$V2)])

# Mantel test  -----------------
library(ade4)

# population
mantel.rtest(internal2, external2, nrepet = 9999)
mantel.rtest(external2, combined2, nrepet = 9999)
mantel.rtest(combined2, diverse2, nrepet = 9999)

mantel.rtest(internal2, combined2, nrepet = 9999)
mantel.rtest(internal2, diverse2, nrepet = 9999)
mantel.rtest(external2, diverse2, nrepet = 9999)

# validation samples
mantel.rtest(HD2, internal2, nrepet = 9999)
mantel.rtest(HD2, external2, nrepet = 9999)
mantel.rtest(HD2, combined2, nrepet = 9999)
mantel.rtest(HD2, diverse2, nrepet = 9999)
