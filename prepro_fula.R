library(biomformat)
library(stringr)
library(ape)

setwd("/scratch/fmorandi/fula_test/Fula-microbiome")

datapath = "./02_pipeline_output/exports/"
outpath = "./03_processed_data/"

##### PREP METADATA #####

# Base metadata
meta = read.table(paste0(outpath, "metadata.tsv"), 
                  row.names=1, fill=T, header = T, sep="\t")
# Read counts
meta_dada2 = read.table(paste0(datapath,"stats.tsv"), 
                        row.names=1, fill=T, header = T, sep="\t")
# Make single metadata table, tidy up
meta = merge(meta, meta_dada2[, c(1, 7)], by="row.names")
vars = colnames(meta)
vars[(length(vars)-1):length(vars)] = c("RawPairs", "FilteredPairs")
colnames(meta) = vars
rownames(meta) = meta$Row.names
meta = meta[, -c(1)]
# Nicer order
meta = with(meta, meta[order(AdultChild, TimePoint, Index), ])
# Empty means NA for these columns
meta[meta$Gender == "", "Gender"] = NA
meta[meta$Ingred == "", "Ingred"] = NA
meta[meta$FirstMeal == "", "FirstMeal"] = NA
meta[meta$Parasites == "", "Parasites"] = NA
# Type conversions
meta$Location = as.factor(meta$Location)
meta$SubLocation = as.factor(meta$SubLocation)
meta$AdultChild = as.factor(meta$AdultChild)
meta$TimePoint = as.factor(meta$TimePoint)
meta$Gender = as.factor(meta$Gender)
meta$FirstMeal = as.factor(meta$FirstMeal)
meta$EntColi = as.logical(sapply(meta$Parasites, grepl, pattern="Ent\\. Coli"))
for (col in c("DisResp", "DisDiab", "DisHBP", "DisLBP", "DisCardio", "DisCancer", "DisConstip", "DisDiarr")) {
  meta[col] = as.logical(meta[[col]])
}
# Grouping variables
meta$Group1 = interaction(meta$Location, meta$AdultChild)
meta$Group2 = as.character(meta$Group1)
meta[meta$AdultChild == "Child", "Group2"] = as.character(interaction(
  as.character(meta[meta$AdultChild == "Child", "Group2"]),
  as.character(meta[meta$AdultChild == "Child", "TimePoint"])))
meta$Group2 = as.factor(meta$Group2)
meta$Group3 = interaction(meta$Group1, meta$TimePoint)
meta$Group4 = interaction(meta$AdultChild, meta$TimePoint)
meta$TimePairs = as.factor(interaction(meta$Index, meta$AdultChild))

rm(meta_dada2, vars)

##### PREP ABUNDANCE AND TAXA: PHYLUM #####

biom = read_biom(paste0(datapath, "table-phylum.biom"))
ab_phy = data.frame(t(as.matrix(biom_data(biom))))
ab_phy = ab_phy[order(rownames(ab_phy)), ]
ab_phy = ab_phy[,order(colnames(ab_phy))]
taxa_strings = as.data.frame(do.call(rbind, biom$rows))
taxa_strings = unlist(taxa_strings$id)
taxa_strings = sort(taxa_strings)
taxa_phy = str_split_fixed(taxa_strings, pattern = ";", n = 2)
taxa_phy = data.frame(taxa_phy)
colnames(taxa_phy) = c("Domain", "Phylum")
taxa_phy$String = apply(taxa_phy, 1, paste, collapse=".")
rownames(taxa_phy) = paste0("phy", rownames(taxa_phy))
colnames(ab_phy) = rownames(taxa_phy)

##### PREP ABUNDANCE AND TAXA: GENUS #####

biom = read_biom(paste0(datapath, "table-genus.biom"))
ab_gen = data.frame(t(as.matrix(biom_data(biom))))
ab_gen = ab_gen[order(rownames(ab_gen)), ]
ab_gen = ab_gen[,order(colnames(ab_gen))]
taxa_strings = as.data.frame(do.call(rbind, biom$rows))
taxa_strings = unlist(taxa_strings$id)
taxa_strings = sort(taxa_strings)
taxa_gen = str_split_fixed(taxa_strings, pattern = ";", n = 6)
taxa_gen = data.frame(taxa_gen)
colnames(taxa_gen) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
taxa_gen$String = apply(taxa_gen, 1, paste, collapse=".")
rownames(taxa_gen) = paste0("gen", rownames(taxa_gen))
colnames(ab_gen) = rownames(taxa_gen)

##### PREP ABUNDANCE AND TAXA: SPECIES #####

biom = read_biom(paste0(datapath, "table-species.biom"))
ab_spe = data.frame(t(as.matrix(biom_data(biom))))
ab_spe = ab_spe[order(rownames(ab_spe)), ]
ab_spe = ab_spe[,order(colnames(ab_spe))]
taxa_strings = as.data.frame(do.call(rbind, biom$rows))
taxa_strings = unlist(taxa_strings$id)
taxa_strings = sort(taxa_strings)
taxa_spe = str_split_fixed(taxa_strings, pattern = ";", n = 7)
taxa_spe = data.frame(taxa_spe)
colnames(taxa_spe) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_spe$String = apply(taxa_spe, 1, paste, collapse=".")
rownames(taxa_spe) = paste0("spe", rownames(taxa_spe))
colnames(ab_spe) = rownames(taxa_spe)

##### PREP ABUNDANCE AND TAXA: ASV #####

biom = read_biom(paste0(datapath, "table-asv.biom"))
ab_asv = data.frame(t(as.matrix(biom_data(biom))))
ab_asv = ab_asv[order(rownames(ab_asv)), ]
asv_ids = as.data.frame(do.call(rbind, biom$rows))
asv_ids = unlist(asv_ids$id)
taxa_asv = as.data.frame(do.call(rbind, lapply(biom$rows, function(x) x$metadata$taxonomy)))
taxa_asv = cbind(asv_ids, taxa_asv)
colnames(taxa_asv) = c("Id", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_asv[taxa_asv == ""] = "__"
taxa_asv$String = apply(taxa_asv[,2:8], 1, paste, collapse=".")
taxa_asv = taxa_asv[order(taxa_asv$String), ]
rownames(taxa_asv) = paste0("asv", 1:nrow(taxa_asv))
colnames(ab_asv) = rownames(taxa_asv)[match(str_sub(colnames(ab_asv), -32, -1), taxa_asv$Id)]
ab_asv = ab_asv[, rownames(taxa_asv)]

##### PREP PHYLOGENETIC TREE #####

tree = ape::read.tree(paste0(datapath, "tree.nwk"))
tree$tip.label = rownames(taxa_asv)[match(tree$tip.label, taxa_asv$Id)]

##### MARK SAMPLES PASSING QC AND SAVE PROTEO % #####

meta[rownames(ab_asv), "PassesQC"] = rowSums(ab_asv) > 14000

ab_phy_norm = 100 * ab_phy / rowSums(ab_phy)
proteo_phy = rownames(taxa_phy)[which(taxa_phy$String == "d__Bacteria.p__Proteobacteria")]
meta[rownames(ab_phy_norm), "Proteobacteria"] = ab_phy_norm[proteo_phy]

##### SAVE DATA #####

ab = list(phy = ab_phy, gen = ab_gen, spe = ab_spe, asv = ab_asv)
taxa = list(phy = taxa_phy, gen = taxa_gen, spe = taxa_spe, asv = taxa_asv)
save(ab, taxa, tree, meta, file = paste0(outpath, "fula_data.Rdata"))
