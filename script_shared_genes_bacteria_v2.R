# Check both number of genes and shared genes by strain
# Change working directory, pan_clusters2_table.matrix and metadata.tsv
#  The latter file should have 5 columns, though the 4 latter are empty: Number, isolation, whyDeleted, COLLECTION_DATE, Proteins_Prokka
# AJPerez, 2019 (updated June 2021)
library(magrittr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)

ab <- "ab"
colores <- "Isolation_source"
#colores <- "CRISPR.Cas"
setwd(paste0("./", ab))

# Palette
c10 <- c("red","green","darkorange","darkturquoise","darkred","grey",
         "magenta","darkgreen","blue","yellow3","black","orange","yellow","red","green3")
s4 <- c(1,2,4,3)

# Metadata
mdata <- read.csv(paste0("metadata_", ab, ".tsv"), sep="\t")
isolation <- mdata[[colores]]
isolation2 <- as.numeric(as.factor(isolation))
legend.cols <- as.numeric(as.factor(levels(isolation)))
whydeleted <- mdata$Removable

# Ngenes vs Shared
p <- ggplot(mdata, aes(Ngenes, Nshared, label=ID))+ 
  theme_gray() +
  stat_smooth(method = "lm", lty = 2, col = "black") +
  geom_point(aes(col = mdata[[colores]], shape = as.factor(Removable)), size=3, stroke=2.5) +
  scale_colour_manual(values = c10) +
  scale_shape_manual(values = s4) +
  xlab("Number of genes") +
  ylab("Average number of shared genes") +
  theme(plot.title = element_text(hjust = -0.05), text = element_text(size=16), legend.text = element_text(size = 13),
        legend.position = "top", legend.box = "horizontal", legend.title = element_blank(), 
        legend.spacing.x = unit(0.1, 'cm')) +
  guides(color = guide_legend(nrow = 2, keyheight = 1.3), shape =  guide_legend(nrow = 2, keyheight = 1.3))
  #scale_x_continuous(breaks = seq(0, max(data_length$ngenes)+1000, by = 500), limits =c(0, max(data_length$ngenes)+1000))
p

# CRISPR-Cas
ggplot(mdata, aes(x = CRISPR.Cas, y = Ngenes)) +
  geom_violin(fill=rainbow(2048),alpha = .3) + 
  geom_boxplot(fill = rainbow(5), alpha = .3, outlier.shape = 20, outlier.colour = "red", outlier.size = 4, width=0.5)
ggplot(mdata, aes(x = CRISPR.Cas, y = Nshared)) +
  geom_violin(fill=rainbow(2048),alpha = .3) + 
  geom_boxplot(fill = rainbow(5), alpha = .3, outlier.shape = 20, outlier.colour = "red", outlier.size = 4, width=0.5)

ggplot(mdata, aes(x = Isolation_source, y = Ngenes)) +
  geom_violin(fill=rainbow(7168),alpha = .3) + 
  geom_boxplot(fill = rainbow(14), alpha = .3, outlier.shape = 20, outlier.colour = "red", outlier.size = 4, width=0.5)
ggplot(mdata, aes(x = Isolation_source, y = Nshared)) +
  geom_violin(fill=rainbow(7168),alpha = .3) + 
  geom_boxplot(fill = rainbow(14), alpha = .3, outlier.shape = 20, outlier.colour = "red", outlier.size = 4, width=0.5)


# Lengths
length2 <- data.frame(n=mdata$Ngenes, total=mdata$Prokka_proteins)
length2 <- stack(length2)

redline1 <- quantile(length2$values, probs = 0.25, na.rm = T) - 1.5 * IQR(length2$values, na.rm = T)
h1 <- ggplot(data=length2, aes(values, fill=ind)) + 
  geom_histogram(binwidth = 1, color="grey20", aes(fill=ind), position = "dodge") +
  xlab("Number of genes") +
  ylab("Number of strains") +
  #scale_x_continuous(limits =c(min(length2$values), max(length2$values))) +
  theme(plot.title = element_text(hjust = -0.1)) +
  geom_density(aes(y=1*..count..), size = 1, alpha = .5) +
  geom_vline(xintercept = redline1, color="red", linetype="dashed", size=1.25) +
  theme_grey() +
  theme(legend.title = element_blank(), legend.position=c(0.82, 0.87), text = element_text(size=16))
h1

redline2 <- quantile(mdata$Nshared, probs = 0.25, na.rm = T) - 1.5 * IQR(mdata$Nshared, na.rm = T)
h2 <- ggplot(mdata, aes(Nshared, label=ID)) + 
  geom_histogram(binwidth = 1, color="grey40", fill="grey") +
  xlab("Average number of shared genes") +
  ylab("Number of strains") +
  #scale_x_continuous(breaks = seq(0, 5000, by = 50), limits =c(min(data_length$shared), max(data_length$shared))) +
  geom_density(aes(y=1*..count..), size = 1.2) +
  geom_vline(xintercept = redline2, color="red", linetype="dashed", size=1.25) +
  theme_grey() +
  theme(plot.title = element_text(hjust = -0.1), text = element_text(size=16))
h2

# Test for bimodality
diptest::dip.test(length2$values)
diptest::dip.test(mdata$Nshared)

# Final figure
ab <- ggarrange(h1, h2, labels = c("a)", "b)"), font.label = list(size = 16, color = "black", face = "bold"), nrow = 2)
pdf("ngenes.pdf", width=16, height=8, paper='special')
ggarrange(ab, p, labels = c("", "c)"), font.label = list(size = 16, color = "black", face = "bold"),  widths = c(1/4, 3/4))
dev.off()

# Strains to remove by either low total genes or low shared genes
g1 <- which(length2$values < redline1)
g2 <- which(data_length$shared < redline2)
removables <- sort(union(g1, g2))
removables <- paste0("pa", formatC(removables, width = 5, format ="d", flag = "0"))
write.table(removables, file = "removables.id", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

