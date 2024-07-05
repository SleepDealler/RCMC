#!/usr/bin/env Rscript

# Ustawienie ścieżki do własnego katalogu na pakiety
.libPaths("./lib")
# Załaduj biblioteki
library(WGCNA, lib.loc = "./lib")
library(ggplot2, lib.loc = "./lib")
library(pheatmap, lib.loc = "./lib")
library(RCy3, lib.loc = "./lib")
library(argparse, lib.loc = "./lib")
library(grid, lib.loc = "./lib")
library(gridExtra, lib.loc = "./lib")

# Definiowanie parsera argumentów
parser <- ArgumentParser(description = "Opis skryptu")

# Dodaj wymagane argumenty
parser$add_argument("-i", "--input", type = "character", help = "Input file", required = TRUE)
parser$add_argument("-m", "--metadata", type = "character", help = "Metadata file", required = TRUE)
parser$add_argument("-p", "--phylo", type = "character", help = "Phylogenetic traits file", required = TRUE)

# Dodaj opcjonalne argumenty
parser$add_argument("-d", "--deepsplit", type = "integer", help = "Deep split", default = 2)
parser$add_argument("-mc", "--minClustersize", type = "integer", help = "Minimum cluster size", default = 20)
parser$add_argument("-ch", "--cutHeight", type = "double", help = "Cut height", default = 0.99)
parser$add_argument("-mdt", "--MeDissThres", type = "double", help = "Module eigengene dissimilarity threshold", default = 0.25)

# Parsowanie argumentów
args <- parser$parse_args()

# Przypisanie plików do zmiennych
counts_t_filtered <- read.csv(args$input, row.names = 1)
metadata <- read.csv(args$metadata)
phylo <- read.csv(args$phylo)

# Sprawdzenie, czy identyfikatory wierszy w counts_t_filtered zgadzają się z wartościami w kolumnie 'Library Name' w metadata
if (!all(rownames(counts_t_filtered) %in% metadata$Library.Name)) {
  stop("Error: Some row identifiers in counts_t_filtered do not match any values in the 'Library Name' column of metadata.")
}

if (!all(metadata$Library.Name %in% rownames(counts_t_filtered))) {
  stop("Error: Some 'Library Name' values in metadata do not match any row identifiers in counts_t_filtered.")
}

# Dynamiczne ładowanie kolumn z pliku phylo.csv
selected_columns <- colnames(phylo)
sample_metadata <- metadata[, selected_columns, drop = FALSE]


# Kontynuacja skryptu
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Wykonaj WGCNA z bardziej szczegółowymi parametrami
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(counts_t_filtered, powerVector = powers, verbose = 5)

# Wykres dopasowania topologii skali
pdf(file = "./plots/scale_free_topology_fit.pdf")
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (R^2)",
     type = "n", main = paste("Scale independence"))
grid()
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
grid()
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# Wybierz moc do konstrukcji sieci
if (!is.na(sft$powerEstimate)) {
  softPower <- sft$powerEstimate
} else {
  softPower <- max(powers)
}

# Stwórz macierze adjacency i TOM na przefiltrowanych danych
adjacency <- adjacency(counts_t_filtered, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM

# Hierarchiczne klastrowanie genów
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Dynamiczne cięcie drzewa z wartościami z parsera
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = args$deepsplit, pamRespectsDendro = TRUE,
                             minClusterSize = args$minClustersize,  
                             cutHeight = args$cutHeight)

# Konwersja etykiet modułów na kolory
dynamicColors <- labels2colors(dynamicMods)

# Wykres dendrogramu z kolorami modułów
pdf(file = "./plots/dendrogram_with_module_colors.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# Obliczanie eigengenów modułów
MEList <- moduleEigengenes(counts_t_filtered, colors = dynamicColors)
MEs <- MEList$eigengenes

# Obliczanie dysymilarności eigengenów modułów
MEDiss <- 1 - cor(MEs)

# Klastrowanie eigengenów modułów
METree <- hclust(as.dist(MEDiss), method = "average")

# Merging close modules z wartościami z parsera
MEDissThres <- args$MeDissThres
abline(h = MEDissThres, col = "red")

merge <- mergeCloseModules(counts_t_filtered, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# Rysowanie dendrogramu z kolorami modułów po scaleniu
pdf(file = "./plots/merged_module_dendrogram.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# Utworzenie katalogu wyjściowego, jeśli nie istnieje
if (!dir.exists("output_for_cytoscape")) {
  dir.create("output_for_cytoscape")
}

# Eksport listy genów starych modułów
for (i in 1:length(merge$oldMEs)) {
  modules = c(substring(names(merge$oldMEs)[i], 3))
  genes = colnames(counts_t_filtered)
  inModule = is.finite(match(dynamicColors, modules))
  modGenes = genes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modGenes, modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/orign_CytoscapeInput-edges-", paste(modules, collapse = "-"), ".txt", sep = ""),
                                 nodeFile = paste("output_for_cytoscape/orign_CytoscapeInput-nodes-", paste(modules, collapse = "-"), ".txt", sep = ""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule])
}

# Eksport listy genów nowych modułów
for (i in 1:length(merge$newMEs)) {
  modules = c(substring(names(merge$newMEs)[i], 3))
  genes = colnames(counts_t_filtered)
  inModule = is.finite(match(mergedColors, modules))
  modGenes = genes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modGenes, modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/merge_CytoscapeInput-edges-", paste(modules, collapse = "-"), ".txt", sep = ""),
                                 nodeFile = paste("output_for_cytoscape/merge_CytoscapeInput-nodes-", paste(modules, collapse = "-"), ".txt", sep = ""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = mergedColors[inModule])
}

# Generowanie heatmapy bez tytułu
p <- pheatmap(MEs, cluster_col=TRUE, cluster_row=TRUE, show_rownames=TRUE, show_colnames=TRUE, fontsize=10, silent=TRUE)

# Tworzenie tytułu jako oddzielny obiekt
title <- textGrob("Heatmap of Old Module Eigen-genes and Samples", gp=gpar(fontsize=30, fontface="bold"))

# Łączenie tytułu i heatmapy
combined_plot <- grid.arrange(title, p$gtable, ncol=1, heights=c(0.1, 1))

# Zapis do pliku PDF
pdf(file="./plots/oldMEs.pdf", height=18, width=18)
grid.draw(combined_plot)
dev.off()

# Dynamiczne ładowanie kolumn z pliku phylo.csv
selected_columns <- as.character(phylo[1, ])
sample_metadata <- metadata[, selected_columns]
colnames(sample_metadata) <- c("SampleID", "Zone")  # Przykład, dostosuj według faktycznych kolumn

# Heatmap of new module eigen-genes and sample trait (e.g., Zone)
col_ann <- sample_metadata[, c("SampleID", "Zone")]
rownames(col_ann) <- col_ann[, 1]
col_ann$SampleID <- NULL

# Upewnij się, że poziomy w Zone są zgodne z rzeczywistymi wartościami
col_ann$Zone <- factor(col_ann$Zone, levels = c("wild type"))

# Dostosowanie ann_color do rzeczywistych wartości w Zone
ann_color <- list("Zone" = c("wild type" = "yellow"))

# Przykład sprawdzenia poziomów
print(levels(col_ann$Zone))
print(names(ann_color$Zone))

data <- data.frame(mergedMEs)
rownames(data) <- rownames(counts_t_filtered)  # Ustawianie nazw wierszy dla heatmapy
data <- data[order(match(rownames(data), rownames(col_ann))), ]


# Generowanie heatmapy bez tytułu
p <- pheatmap(data, cluster_col=TRUE, cluster_row=FALSE, show_rownames=TRUE,
              show_colnames=TRUE, fontsize=10,
              annotation_row=col_ann, annotation_colors=ann_color, silent=TRUE)

# Tworzenie tytułu jako oddzielny obiekt
title <- textGrob("Heatmap of New Module Eigen-genes and Samples", gp=gpar(fontsize=30, fontface="bold"))

# Łączenie tytułu i heatmapy
combined_plot <- grid.arrange(title, p$gtable, ncol=1, heights=c(0.1, 1))

# Zapis do pliku PDF
pdf(file="./plots/newMEs.pdf", height=18, width=18)
grid.draw(combined_plot)
dev.off()

# Define numbers of genes and samples
nGenes = ncol(counts_t_filtered)
nSamples = nrow(counts_t_filtered)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(counts_t_filtered, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

# names (colors) of the modules
modNames = substring(names(MEs), 3)

# If you have any specific traits, you can define them here
# For this example, we'll use a random trait
set.seed(12)
random_trait <- rnorm(nSamples)
Verru = as.data.frame(random_trait)
names(Verru) = "RT"

MET = orderMEs(cbind(MEs, Verru))

geneModuleMembership = as.data.frame(cor(counts_t_filtered, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneTraitSignificance = as.data.frame(cor(counts_t_filtered, Verru, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Verru), sep = "")
names(GSPvalue) = paste("p.GS.", names(Verru), sep = "")

# Select the module with the highest average gene significance
averageGS <- colMeans(abs(geneTraitSignificance))
selectedModule <- modNames[which.max(averageGS)]

# Plot the dendrogram and save to PDF with increased title font size
pdf(file = "./plots/eigengene_dendrogram.pdf", height = 8, width = 12)
par(cex = 1.0)
# Plot the dendrogram without the title
plotEigengeneNetworks(MET, "", marDendro = c(0, 4, 2, 0), plotHeatmaps = FALSE)
# Add the title with increased font size
title(main = "Eigengene dendrogram", cex.main = 2.0)
dev.off()


# Plot the heatmap matrix and save to PDF with adjusted margins
pdf(file = "./plots/eigengene_adjacency_heatmap.pdf", height = 8, width = 12)
par(mar = c(10, 10, 4, 2))  # Adjust margins: c(bottom, left, top, right)
plotEigengeneNetworks(MET, "", marHeatmap = c(3, 3, 3, 3), plotDendrograms = FALSE, xLabelsAngle = 90)
title(main = "Eigengene adjacency heatmap", cex.main = 1.5)  # Adjust title font size
dev.off()

# Select module for detailed analysis
module = selectedModule

# Rename to moduleColors
moduleColors = mergedColors
column = match(module, modNames)
moduleGenes = moduleColors == module

# Plot scatterplot of module membership vs. gene significance and save to PDF
pdf(file = "./plots/module_membership_vs_gene_significance.pdf", height = 7, width = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for RandomTrait",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")
dev.off()

# Utworzenie katalogu wyjściowego, jeśli nie istnieje
if (!dir.exists("results")) {
  dir.create("results")
}

# Dynamicznie wygenerowane ścieżki plików dla wybranego koloru modułu
selected_edge_file <- paste0("output_for_cytoscape/merge_CytoscapeInput-edges-", selectedModule, ".txt")
selected_node_file <- paste0("output_for_cytoscape/merge_CytoscapeInput-nodes-", selectedModule, ".txt")

# Wczytaj pliki dla krawędzi i węzłów dla wybranego modułu
edge <- read.delim(selected_edge_file)
colnames(edge) <- c("source", "target", "weight", "direction", "fromAltName", "toAltName")

node <- read.delim(selected_node_file)
colnames(node) <- c("id", "altName", "node_attributes")

# Zapisywanie wyników do plików edge i node
write.table(edge, file = paste0("results/edges_", selectedModule, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(node, file = paste0("results/nodes_", selectedModule, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
